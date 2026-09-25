# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex


class Sensor_TwoComponentSystem(Mechanism):
    r"""Two-component system membrane sensor with Michaelis-Menten kinetics.

    A 'membrane_sensor' mechanism that models a two-component system (TCS)
    for signal transduction across cellular membranes. This mechanism
    includes signal substrate sensing, membrane sensor protein activation,
    auto-phosphorylation via ATP, and phosphorylation of response proteins,
    but does not include downstream reporter circuits.

    The mechanism follows a multi-step Michaelis-Menten kinetic scheme with
    the following reaction pathway:

    1. Activation of membrane sensor protein (MS):
    $$
        'SP' + 'SigSub' <--> 'SP':'SigSub' == 'SP'^*
    $$

    2. Auto-phosphorylation via ATP:
    $$
        'SP'^* + n 'ATP' <--> 'SP'^*:n'ATP' --> 'SP'^{**}:n'ADP'
            --> 'SP'^{**} + n'ADP'
    $$

    3. Phosphorylation of response protein (RP):
    $$
        'SP'^{**} + 'RP' <--> 'SP'^{**}:'RP' --> 'SP'^*:'RP'^*
            --> 'SP'^* + 'RP'^*
    $$

    4. Product formation:
    $$
        2 'RP'^* <--> 'Product'
    $$

    5. Dephosphorylation of phosphorylated response protein:
    $$
        'RP'^* --> 'RP' + 'Pi'
    $$

    Parameters
    ----------
    name : str, default='sensor_two_component_system'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='membrane_sensor'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('membrane_sensor').

    See Also
    --------
    Mechanism : Base class for all mechanisms.
    MichaelisMenten : Enzyme-substrate mechanism with MM kinetics.

    Notes
    -----
    This mechanism models bacterial two-component systems, which
    are common environmental sensing pathways. The sensor protein spans the
    membrane and undergoes conformational changes upon binding external
    signals, leading to autophosphorylation and subsequent phosphotransfer
    to response proteins that regulate gene expression.

    The mechanism requires the membrane sensor protein to have an ATP
    attribute (membrane_sensor.ATP) that specifies the number of
    ATP molecules required for autophosphorylation.

    Required parameters for this mechanism:

    - 'kb_sigMS' : Forward binding rate for signal substrate to membrane
      sensor protein
    - 'ku_sigMS' : Reverse unbinding rate for signal substrate from
      membrane sensor protein
    - 'kb_autoPhos' : Forward binding rate for ATP to activated membrane
      sensor protein
    - 'ku_autoPhos' : Reverse unbinding rate for ATP from activated
      membrane sensor protein
    - 'k_hydro' : ATP hydrolysis rate constant
    - 'ku_waste' : Unbinding rate for ADP waste products
    - 'kb_phosRP' : Forward binding rate for response protein to
      phosphorylated membrane sensor
    - 'ku_phosRP' : Reverse unbinding rate for response protein from
      phosphorylated membrane sensor
    - 'k_phosph' : Phosphotransfer rate constant to response protein
    - 'ku_activeRP' : Unbinding rate for activated response protein
    - 'kb_dimerRP' : Binding rate for product formation from 
        activated response protein
    - 'ku_dimerRP' : Unbinding rate of the product to the activated
        response protein
    - 'ku_dephos' : Dephosphorylation rate constant for phosphorylated
      response protein

    Examples
    --------
    Create a two-component signaling system with default parameters:

    >>> response = bcp.Protein(name='OmpR')
    >>> sensor = bcp.MembraneSensor(
    ...     membrane_sensor='EnvZ',
    ...     response_protein=response.species,
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Osmolarity',
    ...     ATP=2
    ... )
    >>> mechanism = bcp.Sensor_TwoComponentSystem()
    >>> mixture = bcp.Mixture(
    ...     components=[sensor, response],
    ...     mechanisms={'membrane_sensor': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='sensor_two_component_system',
        mechanism_type='membrane_sensor',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        membrane_sensor,
        response_protein,
        assigned_substrate,
        signal_substrate,
        product,
        energy,
        waste,
        complex_dict=None,
        **kwargs,
    ):
        """Generate species for two-component membrane signaling pathway.

        Creates all species involved in the signaling cascade, including the
        membrane sensor protein, response protein, signal substrate, ATP/ADP
        energy species, and all intermediate complexes formed during signal
        transduction and phosphotransfer.

        Parameters
        ----------
        membrane_sensor : Species
            The membrane sensor protein that detects the signal. Must have
            an ATP attribute specifying the number of ATP molecules required
            for autophosphorylation.
        response_protein : Species
            The response protein that receives the phosphate group and
            becomes activated.
        assigned_substrate : Species
            The phosphate substrate (Pi) that is transferred during the
            signaling cascade.
        signal_substrate : Species
            The external signal molecule that activates the membrane sensor
            protein.
        product : Species
            The phosphorylated response protein product (RP*).
        energy : Species
            ATP species used for autophosphorylation.
        waste : Species
            ADP species produced after ATP hydrolysis.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species with keys
            'Activated_MS', 'ATP:Activated_MS', 'ADP:Activated_MS:sub',
            'Activated_MS:sub', 'Activated_MS:sub:RP', 'Activated_RP', and
            'Activated_MS:Activated_RP'. If None, complexes are automatically
            created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing individual species and complex array:
            [membrane_sensor, response_protein, assigned_substrate,
            signal_substrate, product, energy, waste, complex_array]
            where complex_array is a list of all Complex species generated.

        Notes
        -----
        The method creates seven different complex species representing the
        intermediate states of the signaling cascade:

        1. Activated_MS : signal_substrate:membrane_sensor
        2. ATP:Activated_MS : nATP:Activated_MS
        3. ADP:Activated_MS:sub : Activated_MS:nADP:assigned_substrate
        4. Activated_MS:sub : Activated_MS:assigned_substrate
            (phosphorylated sensor)
        5. Activated_MS:sub:RP :
            Activated_MS:assigned_substrate:response_protein
        6. Activated_RP : response_protein:assigned_substrate
            (phosphorylated response protein)
        7. Activated_MS:Activated_RP :
            Activated_MS:(response_protein:assigned_substrate)

        The number of ATP/ADP molecules (nATP) is determined by the
        membrane_sensor.ATP attribute.

        """
        nATP = membrane_sensor.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MS'] = Complex(
                [signal_substrate, membrane_sensor],
                compartment=membrane_sensor.compartment,
            )
            # Complex2
            complex_dict['ATP:Activated_MS'] = Complex(
                [nATP * [energy], complex_dict['Activated_MS']],
                compartment=membrane_sensor.compartment,
            )
            # Complex3
            complex_dict['ADP:Activated_MS:sub'] = Complex(
                [complex_dict['Activated_MS'],
                    nATP * [waste],
                    assigned_substrate,
                ],
                compartment=membrane_sensor.compartment,
            )
            # Complex4
            complex_dict['Activated_MS:sub'] = Complex(
                [complex_dict['Activated_MS'], assigned_substrate],
                compartment=membrane_sensor.compartment,
            )
            # Complex5
            complex_dict['Activated_MS:sub:RP'] = Complex(
                [complex_dict['Activated_MS:sub'], response_protein],
                compartment=membrane_sensor.compartment,
            )
            #Complex 6
            complex_dict['Activated_RP'] = Complex(
                [response_protein, assigned_substrate],
                compartment=membrane_sensor.compartment,
            )
            # Complex7
            complex_dict['Activated_MS:Activated_RP'] = Complex(
                [complex_dict['Activated_MS'], complex_dict['Activated_RP']],
                compartment=membrane_sensor.compartment,
            )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [
            membrane_sensor,
            response_protein,
            assigned_substrate,
            signal_substrate,
            product,
            energy,
            waste,
            complex_array,
        ]

    def update_reactions(
        self,
        membrane_sensor,
        response_protein,
        assigned_substrate,
        signal_substrate,
        product,
        energy,
        waste,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for two-component membrane signaling pathway.

        Creates all nine reactions comprising the complete signaling
        cascade from signal detection through response protein activation
        and dephosphorylation. Reactions follow Michaelis-Menten kinetics
        with reversible binding steps and irreversible catalytic steps.

        Parameters
        ----------
        membrane_sensor : Species
            The membrane sensor protein that detects the signal. Must have
            an ATP attribute specifying the number of ATP molecules required
            for autophosphorylation.
        response_protein : Species
            The response protein that receives the phosphate group and
            becomes activated.
        assigned_substrate : Species
            The phosphate substrate (Pi) that is transferred during the
            signaling cascade.
        signal_substrate : Species
            The external signal molecule that activates the membrane sensor
            protein.
        product : Species
            The phosphorylated response protein product (RP*).
        energy : Species
            ATP species used for autophosphorylation.
        waste : Species
            ADP species produced after ATP hydrolysis.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes
            are automatically created using the same logic as in
            update_species.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of nine reactions representing the complete signaling
            cascade:

            1. Signal binding (reversible)
            2. ATP binding (reversible)
            3. ATP hydrolysis (irreversible)
            4. ADP release (irreversible)
            5. Response protein binding (reversible)
            6. Phosphotransfer (irreversible)
            7. Activated response protein release (irreversible)
            8. Active response protein product formation (reversible)
            9. Response protein dephosphorylation (irreversible)

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme follows this pathway:

        1. SP + SigSub <--> SP:SigSub (rates: 'kb_sigMS', 'ku_sigMS')
        2. SP:SigSub + nATP <--> SP:SigSub:nATP
           (rates: 'kb_autoPhos', 'ku_autoPhos')
        3. SP:SigSub:nATP --> SP:SigSub:Pi:nADP (rate: 'k_hydro')
        4. SP:SigSub:Pi:nADP --> SP:SigSub:Pi + nADP (rate: 'ku_waste')
        5. SP:SigSub:Pi + RP <--> SP:SigSub:Pi:RP
           (rates: 'kb_phosRP', 'ku_phosRP')
        6. SP:SigSub:Pi:RP --> SP:SigSub:RP:Pi (rate: 'k_phosph')
        7. SP:SigSub:RP:Pi --> SP:SigSub + RP:Pi (rate: 'ku_activeRP')
        8. 2 RP:Pi <--> Product (rate: 'kb_dimerRP', 'ku_dimerRP')
        9. RP:Pi --> RP + Pi (rate: 'ku_dephos')

        This method requires both component and part_id parameters to
        retrieve rate constants from the component's parameter database.

        """
        # Get Parameters
        kb_sigMS = component.get_parameter(
            'kb_sigMS', part_id=part_id, mechanism=self
        )
        ku_sigMS = component.get_parameter(
            'ku_sigMS', part_id=part_id, mechanism=self
        )
        kb_autoPhos = component.get_parameter(
            'kb_autoPhos', part_id=part_id, mechanism=self
        )
        ku_autoPhos = component.get_parameter(
            'ku_autoPhos', part_id=part_id, mechanism=self
        )
        k_hydro = component.get_parameter(
            'k_hydro', part_id=part_id, mechanism=self
        )
        ku_waste = component.get_parameter(
            'ku_waste', part_id=part_id, mechanism=self
        )
        kb_phosRP = component.get_parameter(
            'kb_phosRP', part_id=part_id, mechanism=self
        )
        ku_phosRP = component.get_parameter(
            'ku_phosRP', part_id=part_id, mechanism=self
        )
        k_phosph = component.get_parameter(
            'k_phosph', part_id=part_id, mechanism=self
        )
        ku_activeRP = component.get_parameter(
            'ku_activeRP', part_id=part_id, mechanism=self
        )
        kb_dimerRP = component.get_parameter(
            'kb_dimerRP', part_id=part_id, mechanism=self
        )
        ku_dimerRP = component.get_parameter(
            'ku_dimerRP', part_id=part_id, mechanism=self
        )
        ku_dephos = component.get_parameter(
            'ku_dephos', part_id=part_id, mechanism=self
        )

        # Complexes
        nATP = membrane_sensor.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MS'] = Complex(
                [signal_substrate, membrane_sensor],
                compartment=membrane_sensor.compartment,
            )
            # Complex2
            complex_dict['ATP:Activated_MS'] = Complex(
                [nATP * [energy], complex_dict['Activated_MS']],
                compartment=membrane_sensor.compartment,
            )
            # Complex3
            complex_dict['ADP:Activated_MS:sub'] = Complex(
                [complex_dict['Activated_MS'],
                    nATP * [waste],
                    assigned_substrate,
                ],
                compartment=membrane_sensor.compartment,
            )
            # Complex4
            complex_dict['Activated_MS:sub'] = Complex(
                [complex_dict['Activated_MS'], assigned_substrate],
                compartment=membrane_sensor.compartment,
            )
            # Complex5
            complex_dict['Activated_MS:sub:RP'] = Complex(
                [complex_dict['Activated_MS:sub'], response_protein],
                compartment=membrane_sensor.compartment,
            )
            #Complex 6
            complex_dict['Activated_RP'] = Complex(
                [response_protein, assigned_substrate],
                compartment=membrane_sensor.compartment,
            )
            # Complex7
            complex_dict['Activated_MS:Activated_RP'] = Complex(
                [complex_dict['Activated_MS'], complex_dict['Activated_RP']],
                compartment=membrane_sensor.compartment,
            )

        # Two-component signal transduction
        # Activation of membrane sensor: S + P<--> P*
        binding_rxn1 = Reaction.from_massaction(
            inputs=[signal_substrate, membrane_sensor],
            outputs=[complex_dict['Activated_MS']],
            k_forward=kb_sigMS,
            k_reverse=ku_sigMS,
        )

        # Auto-phosphorylation membrane sensor:
        # P* + ATP<--> P*:ATP
        binding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MS'], nATP * [energy]],
            outputs=[complex_dict['ATP:Activated_MS']],
            k_forward=kb_autoPhos,
            k_reverse=ku_autoPhos,
        )
        # P*:ATP--> P*:Pi:ADP
        hydrolysis_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['ATP:Activated_MS']],
            outputs=[complex_dict['ADP:Activated_MS:sub']],
            k_forward=k_hydro,
        )
        # P*:Pi:ADP--> P*:Pi +ADP
        unbinding_rxn3 = Reaction.from_massaction(
            inputs=[complex_dict['ADP:Activated_MS:sub']],
            outputs=[complex_dict['Activated_MS:sub'], nATP * [waste]],
            k_forward=ku_waste,
        )

        # Phosphorylation of response protein:
        # P*:Pi + RP <--> P*:Pi:RP
        binding_rxn4 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MS:sub'], response_protein],
            outputs=[complex_dict['Activated_MS:sub:RP']],
            k_forward=kb_phosRP,
            k_reverse=ku_phosRP,
        )
        # P*:Pi:RP --> P*:RP:Pi
        Phosph_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MS:sub:RP']],
            outputs=[complex_dict['Activated_MS:Activated_RP']],
            k_forward=k_phosph,
        )
        # P*:RP:Pi--> P* + RP:Pi
        unbinding_rxn5 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MS:Activated_RP']],
            outputs=[
                complex_dict['Activated_RP'],
                complex_dict['Activated_MS'],
            ],
            k_forward=ku_activeRP,
        )
        # Product formation: 2 RP:Pi --> Product
        binding_rxn6 = Reaction.from_massaction(
            inputs=[2*[complex_dict['Activated_RP']]],
            outputs=[product],
            k_forward=kb_dimerRP,
            k_reverse=ku_dimerRP
        )
        # Dephosphorylation: RP:Pi--> RP + Pi
        unbinding_rxn6 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_RP']],
            outputs=[response_protein, assigned_substrate],
            k_forward=ku_dephos,
        )

        return [
            binding_rxn1,
            binding_rxn2,
            hydrolysis_rxn1,
            unbinding_rxn3,
            binding_rxn4,
            Phosph_rxn1,
            unbinding_rxn5,
            binding_rxn6,
            unbinding_rxn6,
        ]
