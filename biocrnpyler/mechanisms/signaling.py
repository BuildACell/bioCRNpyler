
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex
from ..core.parameter import Parameter


class Membrane_Signaling_Pathway_MM(Mechanism):
    """A mechanism to model a two-component system (TCS) membrane sensor.
    Includes the sensing of signal substrate (SigSub) and the phosphorylation of the response protein (RP) but not include the reporter circuit.
    Mechanism follows Michaelis-Menten Type Reactions.
    Mechanism for the activation of the sensor protein (SP): SP + SigSub <--> SP:SigSub --> SP*
    Mechanism for the auto-phosphorylation: SP* + nATP <--> SP*:nATP --> SP**:nADP --> SP** + nADP
    Mechanism for the phosphorylation of response protein: SP** + RP <--> SP**:RP --> SP*:RP* --> SP*+RP*
    Mechanism for the dephosphorylation of phosphoryled response protein: RP* --> RP + Pi
    """

    def __init__(self, name="two_component_membrane_signaling",
                 mechanism_type="membrane_sensor", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, membrane_sensor_protein, response_protein, assigned_substrate,
                       signal_substrate, product, energy, waste, complex_dict=None,
                       **keywords):

        nATP = membrane_sensor_protein.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MP'] = Complex(
                [signal_substrate, membrane_sensor_protein])
            # Complex2
            complex_dict['ATP:Activated_MP'] = Complex(
                [nATP*[energy], complex_dict['Activated_MP']])
            # Complex3
            complex_dict['ADP:Activated_MP:Sub'] = Complex(
                [complex_dict['Activated_MP'], nATP*[waste], assigned_substrate])
            # Complex4
            complex_dict['Activated_MP:Sub'] = Complex(
                [complex_dict['Activated_MP'], assigned_substrate])
            # Complex5
            complex_dict['Activated_MP:Sub:RP'] = Complex(
                [complex_dict['Activated_MP:Sub'], response_protein])
            # Complex6
            complex_dict['Activated_MP:RP:Sub'] = Complex(
                [complex_dict['Activated_MP'], response_protein, assigned_substrate])

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [membrane_sensor_protein, response_protein, assigned_substrate, signal_substrate, energy, waste, complex_array]

    def update_reactions(self, membrane_sensor_protein, response_protein, assigned_substrate,
                         signal_substrate, product, energy, waste, complex_dict=None,
                         component=None, part_id=None, **keywords):
        """This always requires the inputs component and part_id to find the relevant parameters"""

        # Get Parameters
        kb_sigMS = component.get_parameter(
            "kb_sigMS", part_id=part_id, mechanism=self)
        ku_sigMS = component.get_parameter(
            "ku_sigMS", part_id=part_id, mechanism=self)
        kb_autoPhos = component.get_parameter(
            "kb_autoPhos", part_id=part_id, mechanism=self)
        ku_autoPhos = component.get_parameter(
            "ku_autoPhos", part_id=part_id, mechanism=self)
        k_hydro = component.get_parameter(
            "k_hydro", part_id=part_id, mechanism=self)
        ku_waste = component.get_parameter(
            "ku_waste", part_id=part_id, mechanism=self)
        kb_phosRP = component.get_parameter(
            "kb_phosRP", part_id=part_id, mechanism=self)
        ku_phosRP = component.get_parameter(
            "ku_phosRP", part_id=part_id, mechanism=self)
        k_phosph = component.get_parameter(
            "k_phosph", part_id=part_id, mechanism=self)
        ku_activeRP = component.get_parameter(
            "ku_activeRP", part_id=part_id, mechanism=self)
        ku_dephos = component.get_parameter(
            "ku_dephos", part_id=part_id, mechanism=self)

        # Complexes
        nATP = membrane_sensor_protein.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MP'] = Complex(
                [signal_substrate, membrane_sensor_protein])
            # Complex2
            complex_dict['ATP:Activated_MP'] = Complex(
                [nATP*[energy], complex_dict['Activated_MP']])
            # Complex3
            complex_dict['ADP:Activated_MP:Sub'] = Complex(
                [complex_dict['Activated_MP'], nATP*[waste], assigned_substrate])
            # Complex4
            complex_dict['Activated_MP:Sub'] = Complex(
                [complex_dict['Activated_MP'], assigned_substrate])
            # Complex5
            complex_dict['Activated_MP:Sub:RP'] = Complex(
                [complex_dict['Activated_MP:Sub'], response_protein])
            # Complex6
            complex_dict['Activated_MP:RP:Sub'] = Complex(
                [complex_dict['Activated_MP'], response_protein, assigned_substrate])

    # Two-component signal transduction
        # Activation of membrane sensor: S + P<--> P*
        binding_rxn1 = Reaction.from_massaction(inputs=[signal_substrate, membrane_sensor_protein],
                                                outputs=[
                                                    complex_dict['Activated_MP']],
                                                k_forward=kb_sigMS,
                                                k_reverse=ku_sigMS)

        # Auto-phosphorylation membrane sensor:
        # P* + ATP<--> P*:ATP
        binding_rxn2 = Reaction.from_massaction(inputs=[complex_dict['Activated_MP'], nATP*[energy]],
                                                outputs=[
                                                    complex_dict['ATP:Activated_MP']],
                                                k_forward=kb_autoPhos,
                                                k_reverse=ku_autoPhos)
        # P*:ATP--> P*:Pi:ADP
        hydrolysis_rxn1 = Reaction.from_massaction(inputs=[complex_dict['ATP:Activated_MP']],
                                                   outputs=[
                                                       complex_dict['ADP:Activated_MP:Sub']],
                                                   k_forward=k_hydro)
        # P*:Pi:ADP--> P*:Pi +ADP
        unbinding_rxn3 = Reaction.from_massaction(inputs=[complex_dict['ADP:Activated_MP:Sub']],
                                                  outputs=[
                                                      complex_dict['Activated_MP:Sub'], nATP*[waste]],
                                                  k_forward=ku_waste)

        # Phosphorylation of response protein:
        # P*:Pi + RP <--> P*:Pi:RP
        binding_rxn4 = Reaction.from_massaction(inputs=[complex_dict['Activated_MP:Sub'], response_protein],
                                                outputs=[
                                                    complex_dict['Activated_MP:Sub:RP']],
                                                k_forward=kb_phosRP,
                                                k_reverse=ku_phosRP)
        # P*:Pi:RP --> P*:RP:Pi
        Phosph_rxn1 = Reaction.from_massaction(inputs=[complex_dict['Activated_MP:Sub:RP']],
                                               outputs=[
                                                   complex_dict['Activated_MP:RP:Sub']],
                                               k_forward=k_phosph)
        # P*:RP:Pi--> P* + RP:Pi
        unbinding_rxn5 = Reaction.from_massaction(inputs=[complex_dict['Activated_MP:RP:Sub']],
                                                  outputs=[
                                                      product, complex_dict['Activated_MP'],],
                                                  k_forward=ku_activeRP)
        # Dephosphorylation: RP:Pi--> RP + Pi
        unbinding_rxn6 = Reaction.from_massaction(inputs=[product],
                                                  outputs=[
                                                      response_protein, assigned_substrate],
                                                  k_forward=ku_dephos)

        return [binding_rxn1, binding_rxn2, hydrolysis_rxn1, unbinding_rxn3,
                binding_rxn4, Phosph_rxn1, unbinding_rxn5, unbinding_rxn6]
