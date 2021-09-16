
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import itertools as it

from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex, Species, WeightedSpecies
from .propensities import ProportionalHillNegative, ProportionalHillPositive, GeneralPropensity, ParameterEntry
from .parameter import ModelParameter, Parameter, ParameterEntry

class Membrane_Protein_Integration(Mechanism):
    """A Mechanism to integrate into the membrane."""
    def __init__(self, name= "membrane_protein_integration", mechanism_type="catalysis", **keywords):
        
        """Initializes the integration of protein into membrane based on the negative hill coefficient.
        
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
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
    
    
    def update_reactions(self, integral_membrane_protein, product, component=None, part_id=None, complex=None, complex2 = None, kd=None, kb=None, ku=None,
                         kcat=None, kex=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name
        ##############################
        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        
        if kd is None:
            kd1 = component.get_parameter("kd1", part_id = part_id, mechanism = self)
        else:
            kd1 = kb

        size=integral_membrane_protein.size
        if size>1:
            if kd is None:
                kd2 = component.get_parameter("kd2", part_id = part_id, mechanism = self)
            else:
                kd2 = kb

            if kb is None:
                kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
            else:
                kb1 = kb
            if ku is None:
                ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            else:
                ku1= ku
        #############################   
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        else:
            kcat= kcat
        #############################   
        if kex is None:
            kex = component.get_parameter("kex", part_id = part_id, mechanism = self)
            
        else:
            kex = kex
            
        ##############################   
        
        if complex is None:
            size=integral_membrane_protein.size
            if size > 1:
                complex1=Complex([integral_membrane_protein]*size)
            else: complex1=complex
        else: complex1=complex
        
        #CREATE OLIGOMER!!!
        
        # homo-->0
        binding_rxn0 = Reaction.from_massaction(inputs=[integral_membrane_protein],
                                                outputs=[],
                                                k_forward=kd1)        
        if size > 1:
            # homo: monomer --> polymer
            binding_rxn1 = Reaction.from_massaction(inputs=[integral_membrane_protein]*size,
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)
        
            # polymer-->0
            binding_rxn2 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[],
                                                k_forward=kd2)

            #poly-->integrated
            prophill_negative = ProportionalHillNegative(k=kex, d=complex1, K=kcat, n=4, s1=product )
            integration_rxn1 = Reaction([complex1], [product], propensity_type=prophill_negative)
        else:
            #homo-->integrated
            prophill_negative = ProportionalHillNegative(k=kex, d=integral_membrane_protein, K=kcat, n=4, s1=product )
            integration_rxn1 = Reaction([integral_membrane_protein], [product], propensity_type=prophill_negative)
        
        if size > 1:
            return [binding_rxn0, binding_rxn1, binding_rxn2, integration_rxn1]
        else:
            return [binding_rxn0, integration_rxn1]

class Passive_Transport(Mechanism):
    """A Mechanism to model the transport of a substrate through a membrane channel.
    Does not require energy and has undirectional transport, following diffusion rules."""
    
    def __init__(self, name= "passive_membrane_protein_transport", mechanism_type="catalysis", **keywords):
        """Initializes a Passive_Membrane_Protein_Transport instance
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
        Mechanism.__init__(self, name, mechanism_type)
        
    def update_species(self, membrane_channel, substrate, product, complex=None, complex2=None, **keywords):
        if complex is None:
            complex1=complex
        else: complex1=complex
            
        if complex2 is None:
            complex2=complex2
        else: complex2=complex2
            
        return [membrane_channel, substrate, product, complex1, complex2]
    
    
    def update_reactions(self, membrane_channel, substrate, product, component=None, part_id=None, complex=None, complex2 = None, 
                         k_channel=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name
        ##############################
        if component is None and (k_channel is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if k_channel is None:
            k_transport = component.get_parameter("kb1", part_id = part_id, mechanism = self)
        else:
            k_transport = k_channel

        if complex is None:
            complex1 = complex
        else:
            complex1 = complex
            
        if complex2 is None:
            complex2=complex2
        else: complex2=complex2
        
     
        # Sub (Internal) <--> Product (External)
        diffusion_rxn = Reaction.from_massaction(inputs=[substrate, membrane_channel],
                                           outputs=[product, membrane_channel],
                                           k_forward=k_transport,
                                           k_reverse=k_transport)            
        
        return [diffusion_rxn]

class Facilitated_Passive_Transport(Mechanism):
    """A Mechanism to model the transport of a substrate through a membrane carrier.
        Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.

       In the Copy RXN version, the membrane carrier is not Consumed
       Sub+MC <--> Sub:MC --> Prod:MC --> Prod + MC
    """

    def __init__(self, name= "passive_membrane_protein_transport", mechanism_type="catalysis", **keywords):
        """Initializes a Passive_Membrane_Protein_Transport instance
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
        Mechanism.__init__(self, name, mechanism_type)
     
    def update_species(self, membrane_carrier, Sub, Prod, complex=None, complex2 = None, **keywords):
        if complex is None:
            complex1 = Complex([Sub, membrane_carrier])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([Prod, membrane_carrier])
        else:
            complex2 = complex2

        return [membrane_carrier, Sub, Prod, complex1, complex2]

    def update_reactions(self, membrane_carrier, Sub, Prod, component = None, part_id = None, complex=None, complex2 = None, kb=None, ku=None,
                          **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or k1 is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
        else:
            kb1= kb
            
        if ku is None:
            ku1 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku1 = ku
            ku2 = ku
            

        if complex is None:
            complex1 = Complex([Sub, membrane_carrier])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([Prod, membrane_carrier])
            

        # Sub + MC --> Sub:MC
        k1 = ParameterEntry("k1", 1e-1)
        general = GeneralPropensity(f'k1*{Sub}*{membrane_carrier}-k1*{Prod}*{membrane_carrier}', propensity_species=[Prod,Sub,membrane_carrier], propensity_parameters=[k1])
        cat_rxn = Reaction([Sub, membrane_carrier], [complex1], propensity_type = general)

        # Sub:MC --> Sub + MC
        binding_rxn1 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[membrane_carrier, Sub],
                                                k_forward=ku1)

        # Sub:MC --> Prod:MC
        binding_rxn2 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[complex2],
                                                k_forward=kb1)
        
        # MC:Prod --> MC + Prod
        binding_rxn3 = Reaction.from_massaction(inputs=[complex2],
                                                outputs=[Prod, membrane_carrier],
                                                k_forward=ku2)


        
        return [cat_rxn, binding_rxn1, binding_rxn2, binding_rxn3]

class Primary_Active_Transport(Mechanism):
    """A Mechanism to model the transport of a substrate through a membrane carrier.
        Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.

       In the Copy RXN version, the membrane carrier is not Consumed
       Sub+MT <--> Sub:MT + E --> Sub:MT:E --> MT:Prod:W --> Prod + MT + W
    """

    def __init__(self, name= "passive_membrane_protein_transport", mechanism_type="catalysis", **keywords):
        """Initializes a Passive_Membrane_Protein_Transport instance
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
        Mechanism.__init__(self, name, mechanism_type)
     
    def update_species(self, membrane_pump, Sub, Prod, Energy , Waste, complex=None, complex2 = None, complex3 = None, complex4 = None, **keywords):
              
        if complex is None:
            complex1 = Complex([Sub, membrane_pump])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([Prod, membrane_pump])
        else:
            complex2 = complex2
        
        nATP=membrane_pump.ATP
        
        if complex3 is None:
            complex3 = Complex([nATP*[Energy], complex1])
        else:
            complex3 = complex3
        if complex4 is None:
            complex4 = Complex([nATP*[Waste], complex2])
        else:
            complex4 = complex4
            
        return [membrane_pump, Sub, Prod, complex1, complex2, complex3, complex4]

    def update_reactions(self, membrane_pump, Sub, Prod, Energy, Waste, component = None, part_id = None, complex=None, complex2 = None,
                         complex3=None, complex4 = None, complex5 = None, kb=None, ku=None,
                         kcat=None, kex=None, **keywords):
        #Get Parameters
        nATP=membrane_pump.ATP
        
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
            kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
            kb3 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
            kb4 = component.get_parameter("kb4", part_id = part_id, mechanism = self)
        else:
            kb1, kb2, kb3, kb4= kb*np.ones(4)
            
        if ku is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
            ku3 = component.get_parameter("ku3", part_id = part_id, mechanism = self)
            ku4 = component.get_parameter("ku4", part_id = part_id, mechanism = self)
            ku5 = component.get_parameter("ku5", part_id = part_id, mechanism = self)
            ku6 = component.get_parameter("ku6", part_id = part_id, mechanism = self)
        else:
            ku1, ku2,ku3, ku4, ku5, ku6= ku*np.ones(6)
            
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        else:
            kcat = kcat
            
        if kex is None:
            kex = component.get_parameter("kex", part_id = part_id, mechanism = self) 
        else:
            kex = kex
        

        if complex is None:
            complex1 = Complex([Sub, membrane_pump])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([Prod, membrane_pump])
             
        if complex3 is None:
            complex3 = Complex([nATP*[Energy], complex1])
            
        if complex4 is None:
            complex4 = Complex([nATP*[Waste], complex2])

        if complex5 is None:
            complex5 = Complex([membrane_pump, nATP*[Energy],])
            
        # print(complex1, complex2, complex3, complex4)

        # Sub + MT<--> Sub:MT
        binding_rxn1 = Reaction.from_massaction(inputs=[Sub, membrane_pump],
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)
        
        # Sub:MT <--> Sub:MT:E
        binding_rxn2 = Reaction.from_massaction(inputs=[complex1, nATP*[Energy]],
                                                outputs=[complex3],
                                                k_forward=kb2,
                                                k_reverse=ku2)
        
         # Sub:MT:E --> Prod:MT:W
        binding_rxn3 = Reaction.from_massaction(inputs=[complex3],
                                                outputs=[complex4],
                                                k_forward=kb3)
        # Prod:MT:W --> Prod:MT+W
        binding_rxn4 = Reaction.from_massaction(inputs=[complex4],
                                                outputs=[complex2, nATP*[Waste]],
                                                k_forward=ku3)
        # Prod:MT --> Prod + MT
        binding_rxn5 = Reaction.from_massaction(inputs=[complex2],
                                                outputs=[Prod, membrane_pump],
                                                k_forward=ku5)
        
        #Energy use without transport
        # MT <--> MT:E
        binding_rxn6 = Reaction.from_massaction(inputs=[membrane_pump, nATP*[Energy]],
                                                outputs=[complex5],
                                                k_forward=kb4,
                                                k_reverse=ku4)
        
         # MT:E --> MT+W
        binding_rxn7 = Reaction.from_massaction(inputs=[complex5],
                                                outputs=[membrane_pump, nATP*[Waste]],
                                                k_forward=ku6)
        
        return [binding_rxn1, binding_rxn2, binding_rxn3, binding_rxn4, binding_rxn5,binding_rxn6, binding_rxn7]