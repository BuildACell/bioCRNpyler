# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Species

# from ..core.component import Component  # avoid circular import


class BasicIntegration(Mechanism):
    """Basic DNA integration mechanism without enzyme involvement.

    An 'integration' mechanism that models the direct recombination or
    integration of two DNA molecules to form two new DNA products. This
    represents a simplified integration process without explicit enzyme
    involvement, useful for modeling the overall effect of site-specific
    recombination.

    The integration reaction is given by

    DNA1 + DNA2 --> DNA3 + DNA4

    where DNA1 and DNA2 are input DNA molecules, and DNA3 and DNA4 are
    the recombined product DNA molecules.

    Parameters
    ----------
    name : str, default='basic_integration'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='integration'
        Type classification of this mechanism.
    **kwargs
        Additional keyword arguments (unused, maintained for API
        compatibility).

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('integration').

    See Also
    --------
    EnzymeIntegration : Integration with explicit integrase enzyme.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates a single irreversible mass-action reaction
    representing the overall integration process. It does not model the
    detailed molecular mechanism involving integrase binding or
    intermediate complexes.

    Common applications include:

    - Simplified models of site-specific recombination
    - DNA rearrangement in synthetic biology circuits
    - Cassette exchange systems
    - Genome editing with minimal mechanistic detail

    Required parameters for this mechanism:

    - 'kint' : Integration rate constant

    The mechanism does not generate new species (returns empty list from
    update_species) as it models the direct conversion without
    intermediate complexes.

    Examples
    --------
    See 'Specialized_Tutorials/Integrase_Examples'.

    """

    def __init__(
        self,
        name: str = 'basic_integration',
        mechanism_type: str = 'integration',
        **kwargs,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs=None, **kwargs):
        """Generate species for basic integration.

        This method returns an empty list as basic integration does not
        create intermediate species or complexes.

        Parameters
        ----------
        DNA_inputs : list of Species
            List of input DNA species to be integrated.
        DNA_outputs : list of Species, optional
            List of output DNA species after integration. Not used in
            species generation.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            Empty list (no new species generated).

        Notes
        -----
        This mechanism does not generate intermediate complexes. If
        modeling with explicit integrase-DNA complexes is needed, consider
        using a binding mechanism in combination with integration, or use
        the `EnzymeIntegration` mechanism.

        """
        # this doesn't make any species because I use a Binding
        # mechanism for that maybe if we do the tetramerization
        # mechanism then this would do something
        return []

    def update_reactions(
        self,
        DNA_inputs,
        DNA_outputs,
        component=None,
        part_id=None,
        kint=None,
        **kwargs,
    ):
        """Generate integration reaction for basic integration.

        Creates a single irreversible mass-action reaction for DNA
        integration.

        Parameters
        ----------
        DNA_inputs : list of Species
            List of input DNA species to be integrated (typically 2).
        DNA_outputs : list of Species
            List of output DNA species after integration (typically 2).
        component : Component, optional
            Component containing parameter values. Required if kint is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        kint : Parameter or float, optional
            Integration rate constant. If None, retrieved from component
            parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single irreversible mass-action reaction:
            DNA_inputs --> DNA_outputs.

        Raises
        ------
        ValueError
            If component is None and kint is not provided.

        Notes
        -----
        The reaction follows simple mass-action kinetics with rate
        constant 'kint'. All input DNA species are consumed and all output
        DNA species are produced simultaneously.

        """
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        elif kint is None:
            kint = component.get_parameter(
                'kint', part_id=part_id, mechanism=self
            )

        return [
            Reaction.from_massaction(
                inputs=DNA_inputs, outputs=DNA_outputs, k_forward=kint
            )
        ]


class EnzymeIntegration(Mechanism):
    """Enzyme-catalyzed DNA integration mechanism with integrase.

    An 'integration' mechanism that models site-specific recombination
    catalyzed by integrase enzymes. The mechanism explicitly includes the
    integrase as a catalyst that facilitates DNA recombination but is not
    consumed in the reaction. This mechanism uses tetrameric integrase
    complexes (4 integrase molecules) to catalyze the integration.

    The integration reaction is given by

    4*Int + DNA1 + DNA2 --> 4*Int + DNA3 + DNA4

    where Int is the integrase enzyme, DNA1 and DNA2 are input DNA
    molecules, and DNA3 and DNA4 are the recombined product DNA molecules.

    Parameters
    ----------
    name : str, default='enzyme_integration'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='integration'
        Type classification of this mechanism.
    integrase : str, default='Int1'
        Name of the integrase enzyme. A Species with this name and
        material_type 'protein' is created automatically.
    **kwargs
        Additional keyword arguments (unused, maintained for API
        compatibility).

    Attributes
    ----------
    integrase : Species
        The integrase enzyme species with material_type 'protein'.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('integration').

    See Also
    --------
    BasicIntegration : Integration without explicit enzyme.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models site-specific recombination with explicit
    integrase involvement. The integrase acts as a true catalyst,
    appearing on both sides of the reaction and not being consumed.

    The mechanism uses 4 integrase molecules to model the tetrameric
    integrase complex typically formed during site-specific recombination
    (e.g., lambda phage integration, Cre-loxP recombination). This
    reflects the biological reality where integrase monomers form
    tetramers to catalyze strand exchange.

    Common applications include:

    - Lambda phage integration/excision systems
    - Cre-loxP recombination
    - FLP-FRT site-specific recombination
    - Detailed models of integrase-mediated genome editing

    Required parameters for this mechanism:

    - 'kint' : Integration rate constant

    The stoichiometry of 4 integrase molecules reflects biological
    integrase mechanisms where a tetramer is the active form. The
    mechanism does not generate intermediate complexes (returns empty list
    from update_species).

    Examples
    --------
    See 'Specialized_Tutorials/Integrase_Examples'.

    """

    def __init__(
        self,
        name: str = 'enzyme_integration',
        mechanism_type: str = 'integration',
        integrase='Int1',
        **kwargs,
    ):
        # TODO ZAT: remove unused keywords argument
        self.integrase = Species(name=integrase, material_type='protein')
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs=None, **kwargs):
        """Generate species for enzyme integration.

        This method returns an empty list as enzyme integration does not
        create intermediate species or complexes in this implementation.

        Parameters
        ----------
        DNA_inputs : list of Species
            List of input DNA species to be integrated.
        DNA_outputs : list of Species, optional
            List of output DNA species after integration. Not used in
            species generation.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            Empty list (no new species generated).

        Notes
        -----
        This mechanism does not generate intermediate integrase-DNA
        complexes. The integrase-DNA binding and tetramerization steps are
        implicitly captured in the overall rate constant 'kint'. For more
        detailed mechanistic models, consider using explicit binding
        mechanisms before integration.

        """
        # this doesn't make any species because I use a Binding
        # mechanism for that maybe if we do the tetramerization
        # mechanism then this would do something
        return []

    def update_reactions(
        self,
        DNA_inputs,
        DNA_outputs,
        component=None,
        part_id=None,
        kint=None,
        **kwargs,
    ):
        """Generate integration reaction with integrase enzyme.

        Creates a single irreversible mass-action reaction for enzymatic
        DNA integration with tetrameric integrase complex.

        Parameters
        ----------
        DNA_inputs : list of Species
            List of input DNA species to be integrated (typically 2).
        DNA_outputs : list of Species
            List of output DNA species after integration (typically 2).
        component : Component, optional
            Component containing parameter values. Required if kint is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        kint : Parameter or float, optional
            Integration rate constant. If None, retrieved from component
            parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single irreversible mass-action reaction:
            4*integrase + DNA_inputs --> 4*integrase + DNA_outputs.

        Raises
        ------
        ValueError
            If component is None and kint is not provided.

        Notes
        -----
        The reaction includes 4 integrase molecules on both sides,
        modeling the tetrameric integrase complex required for
        site-specific recombination. The integrase is not consumed (acts
        as a true catalyst).

        The stoichiometry reflects biological mechanisms like:

        - Lambda integrase (Int) tetramer formation
        - Cre recombinase synaptic complex
        - Other integrase family enzymes requiring multimers

        """
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")

        elif kint is None:
            kint = component.get_parameter(
                'kint', part_id=part_id, mechanism=self
            )

        return [
            Reaction.from_massaction(
                inputs=[self.integrase] * 4 + DNA_inputs,
                outputs=[self.integrase] * 4 + DNA_outputs,
                k_forward=kint,
            )
        ]
