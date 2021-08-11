
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import itertools as it

from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex, Species, WeightedSpecies
from .propensities import ProportionalHillNegative, ProportionalHillPositive

class Membrane_Protein_Integration(Mechanism):
    """A Mechanism to integrate into the mebrane."""
    def __init__(self, name= "membrane_protein_integration", mechanism_type="catalysis", **keywords):
        
        """Initializes the integration of protein into membrane based on the negative hill coefficient.
        
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
        Mechanism.__init__(self, name, mechanism_type)
    
    def update_species(self, membrane_channel, product, complex=None, complex2=None, **keywords):
        if complex is None:
            size=membrane_channel.size
            complex1=Complex([membrane_channel]*size)
        else: complex1=complex
                        
        if complex2 is None:
            complex2= None 
        else: complex2=complex2
            
        return [membrane_channel,  product, complex1, complex2]
    
    
    def update_reactions(self, membrane_channel, product, component=None, part_id=None, complex=None, complex2 = None, kd=None, kb=None, ku=None,
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
            kb1 = kb
        if kd is None:
            kd2 = component.get_parameter("kd2", part_id = part_id, mechanism = self)
        else:
            kd2 = kb
        if kb is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
        else:
            kb1 = kb
        ##############################
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
            size=membrane_channel.size
            complex1=Complex([membrane_channel]*size)
        else:
            complex1 = complex
        
        #CREATE OLIGOMER!!!
        
        # homo-->0
        binding_rxn0 = Reaction.from_massaction(inputs=[membrane_channel],
                                                outputs=[],
                                                k_forward=kd1)
        
        # homo: monomer --> polymer
        binding_rxn1 = Reaction.from_massaction(inputs=[membrane_channel]*size,
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)
        
        # polymer-->0
        binding_rxn2 = Reaction.from_massaction(inputs=[complex1],
                                                outputs=[],
                                                k_forward=kd2)
        
        #Hill Negative
        prophill_negative = ProportionalHillNegative(k=kex, d=complex1, K=kcat, n=4, s1=product )
        integration_rxn1 = Reaction([complex1], [product], propensity_type=prophill_negative)
        
        return [binding_rxn0, binding_rxn1, binding_rxn2, integration_rxn1]


class Passive_Membrane_Protein_Transport(Mechanism):
    """A Mechanism to model the transport of a substrate through a membrane protein"""
    
    def __init__(self, name= "passive_membrane_protein_transport", mechanism_type="catalysis", **keywords):
        """Initializes a Passive_Membrane_Protein_Transport instance
        :param name: name of the Mechanism, default: passive_membrane_protein_transport
        :param mechanism_type: type of the Mechanism, default: passive_transport
        """
        
        Mechanism.__init__(self, name, mechanism_type)
        
    def update_species(self, membrane_channel, substrate, product, complex=None, complex2=None, **keywords):
        if complex is None:
            complex1=Complex([membrane_channel, substrate])
        else: complex1=complex
        if complex2 is None:
            complex2=Complex([membrane_channel, product])
        else: complex2=complex2
        return [membrane_channel, substrate, product, complex1, complex2]
    
    
    def update_reactions(self, membrane_channel, substrate, product, component=None, part_id=None, complex=None, complex2 = None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name
        ##############################
        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
            kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
        else:
            kb1, kb2 = kb
        ##############################
        if ku is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku1, ku2 = ku
        #############################   
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
            kcat_rev = component.get_parameter("kcat_rev", part_id = part_id, mechanism = self)
        else:
            kcat, kcat_rev = kcat
        ##############################    
        if complex is None:
            complex1 = Complex([membrane_channel, substrate])
        else:
            complex1 = complex
        ##############################
        if complex2 == None:
            complex2 = Complex([membrane_channel, product])
        ##############################
        
        # Sub + Protein <--> Sub:Protein (1)
        binding_rxn1 = Reaction.from_massaction(inputs=[membrane_channel, substrate],
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)
        
        if membrane_channel.material_type == 'Passive':
            # Sub:Protein <--> Protein:Prod
            cat_rxn = Reaction.from_massaction(inputs=[complex1],
                                               outputs=[complex2],
                                               k_forward=kcat,
                                               k_reverse=kcat_rev)
            # Protein:Prod <--> Prod + Protein (2)
            binding_rxn2 = Reaction.from_massaction(inputs=[complex2],
                                                    outputs=[membrane_channel, product],
                                                    k_forward=ku2,
                                                    k_reverse=kb2)

        elif membrane_channel.material_type == 'Importer' or 'Exporter':
            # Sub:Protein --> Protein:Prod
            cat_rxn = Reaction.from_massaction(inputs=[complex1],
                                               outputs=[complex2],
                                               k_forward=kcat)
            # Prob:Protein --> Protein + Prod (2)
            binding_rxn2 = Reaction.from_massaction(inputs=[complex2],
                                                    outputs=[membrane_channel, product],
                                                    k_forward=ku2)
        else:
            print('Membrane channel direction not identified.')
            
        
        return [binding_rxn1, binding_rxn2, cat_rxn]