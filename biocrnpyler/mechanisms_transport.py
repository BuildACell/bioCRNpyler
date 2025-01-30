
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex
from .propensities import ProportionalHillNegative, GeneralPropensity
from .parameter import Parameter, ParameterEntry

class Simple_Diffusion(Mechanism):
    """A mechanism to model the diffusion of a substrate through a membrane channel.
    Does not require energy and follows diffusion rules.
    Reaction schema: substrate <-> product
    """
    def __init__(self, name= "simple_diffusion", 
                 mechanism_type="diffusion", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
    
    def update_species(self, substrate, product, **keywords):
        return [substrate, product]
    def update_reactions(self, substrate, product, component=None, part_id=None,
                         k_diff=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_diff is None):
            raise ValueError("Must pass in a Component or values for k_diff.")
        if k_diff is None:
            k_diff = component.get_parameter("k_diff", part_id = part_id, mechanism = self)
        else:
            k_diff = k_diff
       
        # Sub (Internal) <--> Product (External)
        diffusion_rxn = Reaction.from_massaction(inputs=[substrate],
                                           outputs=[product],
                                           k_forward=k_diff,
                                           k_reverse=k_diff)            
        return [diffusion_rxn]

class Membrane_Protein_Integration(Mechanism):
    """A simple mehanism to integrate into the membrane protein in the membrane.
    Reaction schema for monomers: monomer -> intergral membrane protein
    Reaction schema for oligomer: monomer*[size] -> oligomer -> intergral membrane protein
    """
    def __init__(self, name= "membrane_protein_integration", 
                 mechanism_type="membrane_insertion", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
    
    def update_species(self, integral_membrane_protein, product, complex=None, complex2=None, **keywords):
        
        if complex is None:
            size=integral_membrane_protein.size
            if size > 1:
                complex1=Complex([integral_membrane_protein]*size)
            else: complex1=complex
        else: complex1=complex
                        
        if complex2 is None:
            complex2= None 
        else: complex2=complex2
            
        return [integral_membrane_protein,  product, complex1, complex2]
    
    def update_reactions(self, integral_membrane_protein, product, component=None, part_id=None, complex=None, complex2 = None, 
                         kb1=None, ku1=None,kb2=None,kex=None, kcat=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name
    
        if component is None and (kb1 is None or ku1 is None or kb2 is None or kcat is None or kex is None):
            raise ValueError("Must pass in a Component or values for kb1, ku1, kb2, kcat, and kex.")
        
        if kb1 is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
        else:
            kb1 = kb1

        size=integral_membrane_protein.size
        if size>1:
            if kb2 is None:
                kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
            else:
                kb2 = kb2

            if ku1 is None:
                ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            else:
                ku1= ku1
   
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        else:
            kcat= kcat

        if kex is None:
            kex = component.get_parameter("kex", part_id = part_id, mechanism = self)
        else:
            kex = kex
            
        if complex is None:
            size=integral_membrane_protein.size
            if size > 1:
                complex1=Complex([integral_membrane_protein]*size)
            else: complex1=complex
        else: complex1=complex
        
        #Integration steps based on if protein is monomer or oligomer   
        if size > 1:
            # homo: monomer --> oligomer
            binding_rxn1 = Reaction.from_massaction(inputs=[integral_membrane_protein]*size,
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)

            #oligomer-->integrated
            prophill_negative = ProportionalHillNegative(k=kex, d=complex1, K=kcat, n=4, s1=product)
            integration_rxn1 = Reaction([complex1], [product], propensity_type=prophill_negative)
        else:
            #monomer-->integrated
            prophill_negative = ProportionalHillNegative(k=kex, d=integral_membrane_protein, K=kcat, n=4, s1=product )
            integration_rxn1 = Reaction([integral_membrane_protein], [product], propensity_type=prophill_negative)
        
        if size > 1:
            return [binding_rxn1, integration_rxn1]
        else:
            return [integration_rxn1]

class Simple_Transport(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane channel.
    Does not require energy and has unidirectional transport, following diffusion rules.
    Reaction schema: membrane_channel + substrate <-> membrane_channel + product
    """
    def __init__(self, name= "simple_membrane_protein_transport", 
                 mechanism_type="transport", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
    
    def update_species(self, membrane_channel, substrate, product, **keywords):
        if membrane_channel.attributes[0] !='Passive':
            raise ValueError("Protein is not classified as a channel with passive transport of small molecules. Use mechanism Facilitated_Passive_Transport instead.")
        
        return [membrane_channel, substrate, product]
    def update_reactions(self, membrane_channel, substrate, product, component=None, part_id=None,
                         k_trnsp=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_trnsp is None):
            raise ValueError("Must pass in a Component or values for k_trnsp.")
        if k_trnsp is None:
            k_trnsp = component.get_parameter("k_trnsp", part_id = part_id, mechanism = self)
        else:
            k_trnsp = k_trnsp
       
        # Sub (Internal) <--> Product (External)
        diffusion_rxn = Reaction.from_massaction(inputs=[substrate, membrane_channel],
                                           outputs=[product, membrane_channel],
                                           k_forward=k_trnsp,
                                           k_reverse=k_trnsp)            
        return [diffusion_rxn]

class Facilitated_Transport_MM(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane carrier.
    Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.
    Mechanism for the schema: Sub+MC <--> Sub:MC --> Prod:MC --> Prod + MC
    """
    def __init__(self, name= "facilitated_membrane_protein_transport", 
                 mechanism_type="transport", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
     
    def update_species(self, membrane_carrier, substrate, product, complex=None, complex2 = None, **keywords):
        if complex is None:
            complex1 = Complex([substrate, membrane_carrier])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([product, membrane_carrier])
        else:
            complex2 = complex2

        return [membrane_carrier, substrate, product, complex1, complex2]

    def update_reactions(self, membrane_carrier, substrate, product, component = None, part_id = None, complex=None, complex2 = None, 
                         k1=None, ku1=None, k_trnsp=None, ku2=None,**keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k1 is None or ku1 is None or ku2 is None or k_trnsp is None):
            raise ValueError("Must pass in a Component or values for k1, ku1, ku2, and k_trnsp.")
        
        if k1 is None:
            k1 = component.get_parameter("k1", part_id = part_id, mechanism = self)
        else:
            k1= ParameterEntry("k1", k1)
            
        if ku1 is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
        else:
            ku1 = ku1

        if k_trnsp is None:
            k_trnsp = component.get_parameter("k_trnsp", part_id = part_id, mechanism = self)
        else:
            k_trnsp = k_trnsp

        if ku2 is None:
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku2 = ku2
            
        if complex is None:
            complex1 = Complex([substrate, membrane_carrier])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([product, membrane_carrier])

        # Sub + MC --> Sub:MC
        general = GeneralPropensity(f'k1*{substrate}*{membrane_carrier}*Heaviside({substrate}-{product})-k1*{product}*{membrane_carrier}*Heaviside({substrate}-{product})', propensity_species=[product,substrate,membrane_carrier], propensity_parameters=[k1])
        binding_rxn1 = Reaction([substrate, membrane_carrier], [complex1], propensity_type = general)
                
        # Sub:MC --> Sub + MC
        unbinding_rxn1 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[membrane_carrier, substrate],
                                                k_forward=ku1)
    
        # Sub:MC --> Prod:MC
        transport_rxn = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[complex2],
                                                k_forward=k_trnsp)
        
        # MC:Prod --> MC + Prod
        unbinding_rxn2 = Reaction.from_massaction(inputs=[complex2],
                                                outputs=[product, membrane_carrier],
                                                k_forward=ku2)
                
        return [binding_rxn1,unbinding_rxn1, transport_rxn, unbinding_rxn2]
        
class Primary_Active_Transport_MM(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane carrier.
    Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.
    Mechanism for the schema: Sub+MT <--> Sub:MT + E --> Sub:MT:E --> MT:Prod:E --> Prod + MT:W --> Prod + MT+ W
    """
    def __init__(self, name= "active_membrane_protein_transport", 
                 mechanism_type="transport", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
     
    def update_species(self, membrane_pump, substrate, product, energy , waste, complex=None, complex2 = None, complex3 = None, complex4 = None, **keywords):
              
        nATP=membrane_pump.ATP

        if complex is None:
            complex1 = Complex([substrate, membrane_pump])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([nATP*[energy], complex1])
        else:
            complex2 = complex2        
        if complex3 is None:
            complex3 = Complex([nATP*[energy], product, membrane_pump])
        else:
            complex3 = complex3
        if complex4 is None:
            complex4 = Complex([nATP*[waste], membrane_pump])
        else:
            complex4 = complex4
            
        return [membrane_pump, substrate, product, energy, waste, complex1, complex2, complex3, complex4]

    def update_reactions(self, membrane_pump, substrate, product, energy, waste, component = None, part_id = None, 
                         complex=None, complex2 = None, complex3=None, complex4 = None, complex5 = None, 
                         k1=None, ku1=None, k2=None, ku2=None, k_trnsp=None, ku3=None, ku4=None, **keywords):

        nATP=membrane_pump.ATP
        
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k1 is None or ku1 is None or k2 is None or ku2 is None or k_trnsp is None or ku3 is None or ku4 is None):
            raise ValueError("Must pass in a Component or values for k1, ku1, k2, ku2, k_trnsp, ku3 and ku4.")

        if k1 is None:
            k1 = component.get_parameter("k1", part_id = part_id, mechanism = self)
        else:
            k1 = ParameterEntry("k1", k1)
            
        if ku1 is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
        else:
            ku1 = ku1

        if k2 is None:
            k2 = component.get_parameter("k2", part_id = part_id, mechanism = self)
        else:
            k2 = ParameterEntry("k2", k2)
            
        if ku2 is None:
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku2 = ku2

        if k_trnsp is None:
            k_trnsp = component.get_parameter("k_trnsp", part_id = part_id, mechanism = self)
        else:
            k_trnsp = k_trnsp

        if ku3 is None:
            ku3 = component.get_parameter("ku3", part_id = part_id, mechanism = self)
        else:
            ku3 = ku3

        if ku4 is None:
            ku4 = component.get_parameter("ku4", part_id = part_id, mechanism = self)
        else:
            ku4 = ku4
        
        if complex is None:
            complex1 = Complex([substrate, membrane_pump])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([nATP*[energy], complex1])   
             
        if complex3 is None:
            complex3 = Complex([nATP*[energy], product, membrane_pump])
            
        if complex4 is None:
            complex4 = Complex([nATP*[waste], membrane_pump])

        # Sub + MT<--> Sub:MT
        general1 = GeneralPropensity(f'k1*{substrate}*{membrane_pump}*Heaviside({membrane_pump})', propensity_species=[substrate,membrane_pump], propensity_parameters=[k1])
        binding_rxn1 = Reaction([substrate, membrane_pump], [complex1], propensity_type = general1)
         
        unbinding_rxn1 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[substrate, membrane_pump],
                                                k_forward=ku1)
        
        # Sub:MT + E <--> Sub:MT:E
        general2 = GeneralPropensity(f'k2*{complex1}*{energy}*Heaviside({complex1})', propensity_species=[complex1,energy], propensity_parameters=[k2])
        binding_rxn2 = Reaction([complex1, nATP*[energy]], [complex2], propensity_type = general2)
         
        unbinding_rxn2 = Reaction.from_massaction(inputs=[complex2],
                                                outputs=[complex1, nATP*[energy]],
                                                k_forward=ku2)
        
         # Sub:MT:E --> Prod:MT:E
        transport_rxn = Reaction.from_massaction(inputs=[complex2],
                                                outputs=[complex3],
                                                k_forward=k_trnsp)
        # Prod:MT:E--> Prod+MT:W
        unbinding_rxn3 = Reaction.from_massaction(inputs=[complex3],
                                                outputs=[complex4, product],
                                                k_forward=ku3)
        # MT:W --> MT+W
        unbinding_rxn4 = Reaction.from_massaction(inputs=[complex4],
                                                outputs=[nATP*[waste], membrane_pump],
                                                k_forward=ku4)

        return [binding_rxn1, unbinding_rxn1, binding_rxn2, unbinding_rxn2, transport_rxn, unbinding_rxn3, unbinding_rxn4]