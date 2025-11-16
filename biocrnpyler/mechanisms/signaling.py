# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex


class Membrane_Signaling_Pathway_MM(Mechanism):
    r"""Two-component system membrane sensor with Michaelis-Menten kinetics.

    A 'membrane_sensor' mechanism that models a two-component system (TCS)
    for signal transduction across cellular membranes. This mechanism
    includes signal substrate sensing, membrane sensor protein activation,
    auto-phosphorylation via ATP, and phosphorylation of response proteins,
    but does not include downstream reporter circuits.

    The mechanism follows a multi-step Michaelis-Menten kinetic scheme with
    the following reaction pathway:

    1. Activation of membrane sensor protein (MSP):

    $$SP + SigSub \leftrightarrow SP:SigSub \rightarrow SP^*$$
    
    2. Auto-phosphorylation via ATP:

    $$SP^* + nATP \leftrightarrow SP^*:nATP \rightarrow SP^{**}:nADP \rightarrow SP^{**} + nADP$$
    
    3. Phosphorylation of response protein (RP):

    $$SP^{**} + RP \leftrightarrow SP^{**}:RP \rightarrow SP^*:RP^* \rightarrow SP^* + RP^*$$

    4. Dephosphorylation of phosphorylated response protein:

    $$RP^* \rightarrow RP + Pi$$

    Parameters
    ----------
    name : str, default='two_component_membrane_signaling'
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
    This mechanism models bacterial two-component signaling systems, which
    are common environmental sensing pathways. The sensor protein spans the
    membrane and undergoes conformational changes upon binding external
    signals, leading to autophosphorylation and subsequent phosphotransfer
    to response proteins that regulate gene expression.

    The mechanism requires the membrane sensor protein to have an ATP
    attribute (membrane_sensor_protein.ATP) that specifies the number of
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
    - 'ku_dephos' : Dephosphorylation rate constant for phosphorylated
      response protein

    Examples
    --------
    Create a two-component signaling system with default parameters:

    >>> response = bcp.Protein(name='OmpR')
    >>> sensor = bcp.MembraneSensor(
    ...     membrane_sensor_protein='EnvZ',
    ...     response_protein=response.species,
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Osmolarity',
    ...     ATP=2
    ... )
    >>> mechanism = bcp.Membrane_Signaling_Pathway_MM()
    >>> mixture = bcp.Mixture(
    ...     components=[sensor, response],
    ...     mechanisms={'membrane_sensor': mechanism},
    ...     parameters={
    ...         'kb_sigMS': 1.0, 'ku_sigMS': 0.1,
    ...         'kb_autoPhos': 1.0, 'ku_autoPhos': 0.1,
    ...         'k_hydro': 0.5, 'ku_waste': 1.0,
    ...         'kb_phosRP': 1.0, 'ku_phosRP': 0.1,
    ...         'k_phosph': 0.5, 'ku_activeRP': 1.0,
    ...         'ku_dephos': 0.01
    ...     }
    ... )
    ... mixture.compile_crn()

    """

    def __init__(
        self,
        name='two_component_membrane_signaling',
        mechanism_type='membrane_sensor',
        **kwargs,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self,
        membrane_sensor_protein,
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
        membrane_sensor_protein : Species
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
            'Activated_MSP', 'ATP:Activated_MSP', 'ADP:Activated_MSP:Sub',
            'Activated_MSP:Sub', 'Activated_MSP:Sub:RP', and
            'Activated_MSP:RP:Sub'. If None, complexes are automatically
            created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing individual species and complex array:
            [membrane_sensor_protein, response_protein,
            assigned_substrate, signal_substrate, energy, waste,
            complex_array] where complex_array is a list of all Complex
            species generated.

        Notes
        -----
        The method creates six different complex species representing the
        intermediate states of the signaling cascade:

        1. Activated_MSP : signal_substrate:membrane_sensor_protein
        2. ATP:Activated_MSP : nATP:Activated_MSP
        3. ADP:Activated_MSP:Sub : Activated_MSP:nADP:Pi
        4. Activated_MSP:Sub : Activated_MSP:Pi (phosphorylated sensor)
        5. Activated_MSP:Sub:RP : (Activated_MSP:Pi):response_protein
        6. Activated_MSP:RP:Sub : Activated_MSP:(response_protein:Pi)

        The number of ATP/ADP molecules (nATP) is determined by the
        membrane_sensor_protein.ATP attribute.

        """
        nATP = membrane_sensor_protein.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MSP'] = Complex(
                [signal_substrate, membrane_sensor_protein],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex2
            complex_dict['ATP:Activated_MSP'] = Complex(
                [nATP * [energy], complex_dict['Activated_MSP']],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex3
            complex_dict['ADP:Activated_MSP:Sub'] = Complex(
                [
                    complex_dict['Activated_MSP'],
                    nATP * [waste],
                    assigned_substrate,
                ],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex4
            complex_dict['Activated_MSP:Sub'] = Complex(
                [complex_dict['Activated_MSP'], assigned_substrate],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex5
            complex_dict['Activated_MSP:Sub:RP'] = Complex(
                [complex_dict['Activated_MSP:Sub'], response_protein],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex6
            complex_dict['Activated_MSP:RP:Sub'] = Complex(
                [
                    complex_dict['Activated_MSP'],
                    response_protein,
                    assigned_substrate,
                ],
                compartment=membrane_sensor_protein.compartment,
            )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [
            membrane_sensor_protein,
            response_protein,
            assigned_substrate,
            signal_substrate,
            energy,
            waste,
            complex_array,
        ]

    def update_reactions(
        self,
        membrane_sensor_protein,
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

        Creates all eight reactions comprising the complete signaling
        cascade from signal detection through response protein activation
        and dephosphorylation. Reactions follow Michaelis-Menten kinetics
        with reversible binding steps and irreversible catalytic steps.

        Parameters
        ----------
        membrane_sensor_protein : Species
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
            List of eight reactions representing the complete signaling
            cascade:

            1. Signal binding (reversible)
            2. ATP binding (reversible)
            3. ATP hydrolysis (irreversible)
            4. ADP release (irreversible)
            5. Response protein binding (reversible)
            6. Phosphotransfer (irreversible)
            7. Activated response protein release (irreversible)
            8. Response protein dephosphorylation (irreversible)

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
        8. RP:Pi --> RP + Pi (rate: 'ku_dephos')

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
        ku_dephos = component.get_parameter(
            'ku_dephos', part_id=part_id, mechanism=self
        )

        # Complexes
        nATP = membrane_sensor_protein.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Activated_MSP'] = Complex(
                [signal_substrate, membrane_sensor_protein],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex2
            complex_dict['ATP:Activated_MSP'] = Complex(
                [nATP * [energy], complex_dict['Activated_MSP']],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex3
            complex_dict['ADP:Activated_MSP:Sub'] = Complex(
                [
                    complex_dict['Activated_MSP'],
                    nATP * [waste],
                    assigned_substrate,
                ],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex4
            complex_dict['Activated_MSP:Sub'] = Complex(
                [complex_dict['Activated_MSP'], assigned_substrate],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex5
            complex_dict['Activated_MSP:Sub:RP'] = Complex(
                [complex_dict['Activated_MSP:Sub'], response_protein],
                compartment=membrane_sensor_protein.compartment,
            )
            # Complex6
            complex_dict['Activated_MSP:RP:Sub'] = Complex(
                [
                    complex_dict['Activated_MSP'],
                    response_protein,
                    assigned_substrate,
                ],
                compartment=membrane_sensor_protein.compartment,
            )

        # Two-component signal transduction
        # Activation of membrane sensor: S + P<--> P*
        binding_rxn1 = Reaction.from_massaction(
            inputs=[signal_substrate, membrane_sensor_protein],
            outputs=[complex_dict['Activated_MSP']],
            k_forward=kb_sigMS,
            k_reverse=ku_sigMS,
        )

        # Auto-phosphorylation membrane sensor:
        # P* + ATP<--> P*:ATP
        binding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MSP'], nATP * [energy]],
            outputs=[complex_dict['ATP:Activated_MSP']],
            k_forward=kb_autoPhos,
            k_reverse=ku_autoPhos,
        )
        # P*:ATP--> P*:Pi:ADP
        hydrolysis_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['ATP:Activated_MSP']],
            outputs=[complex_dict['ADP:Activated_MSP:Sub']],
            k_forward=k_hydro,
        )
        # P*:Pi:ADP--> P*:Pi +ADP
        unbinding_rxn3 = Reaction.from_massaction(
            inputs=[complex_dict['ADP:Activated_MSP:Sub']],
            outputs=[complex_dict['Activated_MSP:Sub'], nATP * [waste]],
            k_forward=ku_waste,
        )

        # Phosphorylation of response protein:
        # P*:Pi + RP <--> P*:Pi:RP
        binding_rxn4 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MSP:Sub'], response_protein],
            outputs=[complex_dict['Activated_MSP:Sub:RP']],
            k_forward=kb_phosRP,
            k_reverse=ku_phosRP,
        )
        # P*:Pi:RP --> P*:RP:Pi
        Phosph_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MSP:Sub:RP']],
            outputs=[complex_dict['Activated_MSP:RP:Sub']],
            k_forward=k_phosph,
        )
        # P*:RP:Pi--> P* + RP:Pi
        unbinding_rxn5 = Reaction.from_massaction(
            inputs=[complex_dict['Activated_MSP:RP:Sub']],
            outputs=[
                product,
                complex_dict['Activated_MSP'],
            ],
            k_forward=ku_activeRP,
        )
        # Dephosphorylation: RP:Pi--> RP + Pi
        unbinding_rxn6 = Reaction.from_massaction(
            inputs=[product],
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
            unbinding_rxn6,
        ]
