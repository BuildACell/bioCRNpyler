
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import itertools as it

from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex, Species, WeightedSpecies
from .propensities import ProportionalHillNegative, ProportionalHillPositive, GeneralPropensity, ParameterEntry
from .parameter import ModelParameter, Parameter, ParameterEntry


class Membrane_Signaling_Pathway_MM(Mechanism):
    """A mechanism to model a two-component system (TCS) membrane sensor.
    Includes the sensing of signal substrate (SigSub) and the phosphorylation of the response protein (RP) but not include the reporter circuit.
    Mechanism follows Michaelis-Menten Type Reactions.
    Mechanism for the activation of the sensor protein (SP): SP + SigSub <--> SP:SigSub --> SP*
    Mechanism for the auto-phosphorylation: SP* + nATP <--> SP*:nATP --> SP**:nADP --> SP** + nADP
    Mechanism for the phosphorylation of response protein: SP** + RP <--> SP**:RP --> SP*:RP* --> SP*+RP*
    Mechanism for the dephosphorylation of phosphoryled response protein: RP* --> RP + Pi
    """
    def __init__(self, name= "two_component_membrane_signaling", 
                 mechanism_type="membrane_sensor", **keywords):
        Mechanism.__init__(self, name, mechanism_type)
    
    def update_species(self, membrane_sensor_protein, response_protein, 
                       assigned_substrate='P', signal_substrate, energy, waste,
                       complex=None, complex2=None, complex3=None, 
                       complex4=None, complex5=None, complex6=None,
                       complex7=None,**keywords):
        
        if complex is None:
            complex1 = Complex([signal_substrate, membrane_sensor_protein])
        else:
            complex1 = complex

        nATP=membrane_sensor_protein.ATP
        
        if complex2 is None:
            complex2 = Complex([nATP*[energy], complex1])
        else:
            complex2 = complex2
        
        if complex3 is None:
            complex3 = Complex([complex1, nATP*[waste], assigned_substrate])
        else:
            complex3 = complex3

        if complex4 is None:
            complex4 = Complex([complex1, assigned_substrate])
        else:
            complex4 = complex4

        if complex5 is None:
            complex5 = Complex([complex4, response_protein])
        else:
            complex5 = complex5

        if complex6 is None:
            complex6 = Complex([complex1, response_protein, assigned_substrate])
        else:
            complex6 = complex6

        if complex7 is None:
            complex7 = Complex([response_protein, assigned_substrate])
        else:
            complex7 = complex7
            
        return [membrane_sensor_protein, response_protein, assigned_substrate, signal_substrate, energy, waste,
                complex1, complex2, complex3, complex4, complex5, complex6, complex7]
    
    def update_reactions(self, membrane_sensor_protein, response_protein, assigned_substrate, signal_substrate, 
                         energy, waste, component=None, part_id=None, 
                         complex=None, complex2 = None, complex3 = None, complex4 = None,
                         complex5 = None, complex6 = None, complex7 = None,
                         kb1=None, ku1=None, kb2=None, ku2=None, k_hydro=None, ku3=None,  kb4=None, ku4=None,
                         k_phosph=None, ku5=None, ku6=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name
    
        if component is None and (kb1 is None or ku1 is None or kb2 is None or kex is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb1, ku1, kb2, kex, and kcat.")
        
        if kb1 is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
        else:
            kb1 = kb1

        if ku1 is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
        else:
            ku1= ku1

        if kb2 is None:
            kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
        else:
            kb2 = kb2

        if ku2 is None:
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku2= ku2

        if k_hydro is None:
            k_hydro = component.get_parameter("k_hydro", part_id = part_id, mechanism = self)
        else:
            k_hydro= k_hydro

        if ku3 is None:
            ku3 = component.get_parameter("ku3", part_id = part_id, mechanism = self)
        else:
            ku3= ku3

        if kb4 is None:
            kb4 = component.get_parameter("kb4", part_id = part_id, mechanism = self)
        else:
            kb4 = kb4

        if ku4 is None:
            ku4 = component.get_parameter("ku4", part_id = part_id, mechanism = self)
        else:
            ku4= ku4

        if k_phosph is None:
            k_phosph = component.get_parameter("k_phosph", part_id = part_id, mechanism = self)
        else:
            k_phosph= k_phosph

        if ku5 is None:
            ku5 = component.get_parameter("ku5", part_id = part_id, mechanism = self)
        else:
            ku5= ku5
        if ku6 is None:
            ku6 = component.get_parameter("ku6", part_id = part_id, mechanism = self)
        else:
            ku6 = ku6
            
        #Complexes
        if complex is None:
            complex1 = Complex([signal_substrate, membrane_sensor_protein])
        else:
            complex1 = complex

        nATP=membrane_sensor_protein.ATP
        
        if complex2 is None:
            complex2 = Complex([nATP*[energy], complex1])
        else:
            complex2 = complex2
        
        if complex3 is None:
            complex3 = Complex([complex1, nATP*[waste], assigned_substrate])
        else:
            complex3 = complex3

        if complex4 is None:
            complex4 = Complex([complex1, assigned_substrate])
        else:
            complex4 = complex4

        if complex5 is None:
            complex5 = Complex([complex4, response_protein])
        else:
            complex5 = complex5

        if complex6 is None:
            complex6 = Complex([complex1, response_protein, assigned_substrate])
        else:
            complex6 = complex6

        if complex7 is None:
            complex7 = Complex([response_protein, assigned_substrate])
        else:
            complex7 = complex7
        
    #Two-component signal transduction  
        # Activation of membrane sensor: S + P--> P*
        binding_rxn1 = Reaction.from_massaction(inputs=[signal_substrate, membrane_sensor_protein],
                                            outputs=[complex1],
                                            k_forward=kb1,
                                            k_reverse=ku1)
        
        # Auto-phosphorylation membrane sensor:
        ## P* + ATP--> P*:ATP
        binding_rxn2 = Reaction.from_massaction(inputs=[complex1, 2*[energy]],
                                            outputs=[complex2],
                                            k_forward=kb2,
                                            k_reverse=ku2)
        ## P*:ATP--> P*:Pi:ADP
        hydrolysis_rxn1 = Reaction.from_massaction(inputs=[complex2],
                                            outputs=[complex3],
                                            k_forward=k_hydro)
        ## P*:Pi:ADP--> P*:Pi +ADP
        unbinding_rxn3 = Reaction.from_massaction(inputs=[complex3],
                                            outputs=[complex4, 2*[waste]],
                                            k_forward=ku3)
        
        #Phosphorylation of response protein:
        ## P*:Pi + RP --> P*:Pi:RP
        binding_rxn4 = Reaction.from_massaction(inputs=[complex4, response_protein],
                                            outputs=[complex5],
                                            k_forward=kb4,
                                            k_reverse=ku4)
        ## P*:Pi:RP --> P*:RP:Pi
        Phosph_rxn1 = Reaction.from_massaction(inputs=[complex5],
                                            outputs=[complex6],
                                            k_forward=k_phosph)
        ## P*:RP:Pi--> P* + RP:Pi
        unbinding_rxn5 = Reaction.from_massaction(inputs=[complex6],
                                            outputs=[complex7, complex1],
                                            k_forward=ku5)
        ##Dephosphorylation: RP:Pi--> RP + Pi
        unbinding_rxn6 = Reaction.from_massaction(inputs=[complex7],
                                            outputs=[response_protein, assigned_substrate],
                                            k_forward=ku6)
        
        return [binding_rxn1, binding_rxn2, hydrolysis_rxn1, unbinding_rxn3,
                binding_rxn4, Phosph_rxn1, unbinding_rxn5, unbinding_rxn6]