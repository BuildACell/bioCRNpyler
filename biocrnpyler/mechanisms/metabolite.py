from ..core.mechanism import Mechanism
from ..core.reaction import Reaction


class OneStepPathway(Mechanism):
    """Simple one-step metabolic pathway mechanism.

    A 'metabolic_pathway' mechanism that models the conversion of precursor
    metabolites to product metabolites in a single irreversible step. This
    mechanism can represent metabolic reactions, spontaneous conversions, or
    simplified enzymatic pathways where enzyme dynamics are not explicitly
    modeled.

    The pathway reaction can be any of the following forms:

    - Precursors --> Products (standard conversion)
    - --> Products (creation from nothing)
    - Precursors --> (degradation to nothing)

    where precursors and products can be single species or lists of species.

    Parameters
    ----------
    name : str, default='one_step_pathway'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='metabolic_pathway'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('metabolic_pathway').

    See Also
    --------
    BasicCatalysis : Enzymatic catalysis with explicit enzyme.
    BasicProduction : Catalytic production mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates a single irreversible mass-action reaction
    with rate constant 'k'. It provides flexibility to model various types
    of metabolic processes:

    - Standard conversion: Multiple precursors convert to multiple
      products
    - Creation: Products appear spontaneously (precursor=None)
    - Degradation: Precursors disappear (product=None)

    The mechanism does not explicitly model enzymes or intermediate
    complexes, making it suitable for:

    - Lumped metabolic pathways
    - Spontaneous chemical reactions
    - Simplified models where enzyme dynamics are negligible
    - Material flow in metabolic networks
    - Constitutive processes

    Common applications include:

    - Metabolic flux balance models
    - Simplified biosynthetic pathways
    - Nutrient uptake and consumption
    - Waste product formation
    - Simple chemical transformations

    Required parameters for this mechanism:

    - 'k' : Forward rate constant for the conversion

    The mechanism supports arbitrary stoichiometries through the use of
    lists for precursors and products. Stoichiometry is determined by the
    number of times a species appears in the list.

    Examples
    --------
    Model a simple metabolic conversion:

    >>> g6p = bcp.Metabolite('g6p', precursors=['glucose'], products=[])
    >>> mixture = bcp.Mixture(
    ...     components=[g6p],
    ...     mechanisms={'metabolic_pathway': bcp.OneStepPathway()},
    ...     parameters={'k': 0.1}
    ... )
    >>> mixture.compile_crn()

    Model constitutive metabolite production:

    >>> metabolite = bcp.Metabolite(
    ...     'metabolite', precursors=[None], products=[])
    >>> mixture = bcp.Mixture(
    ...     components=[metabolite],
    ...     mechanisms={'metabolic_pathway': bcp.OneStepPathway()},
    ...     parameters={'k': 1.0}
    ... )
    >>> mixture.compile_crn()

    Model metabolite degradation:

    >>> waste = bcp.Metabolite(
    ...     'waste_product', precursors=[], products=[None])
    >>> mixture = bcp.Mixture(
    ...     components=[waste],
    ...     mechanisms={'metabolic_pathway': bcp.OneStepPathway()},
    ...     parameters={'k': 0.2}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name='one_step_pathway', mechanism_type='metabolic_pathway'
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, precursor, product, **kwargs):
        """Generate species list for metabolic pathway.

        Collects all species involved in the pathway (precursors and
        products) into a single list.

        Parameters
        ----------
        precursor : Species, list of Species, or None
            Precursor species or list of precursor species. If None, the
            pathway represents creation from nothing.
        product : Species, list of Species, or None
            Product species or list of product species. If None, the
            pathway represents degradation to nothing.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            Combined list of all precursor and product species. Returns
            empty list if both are None.

        Notes
        -----
        This method simply aggregates the input species without creating
        new species or complexes. The species should already exist in the
        system before this mechanism is applied.

        """
        species = []
        if precursor is not None:
            species += precursor
        if product is not None:
            species += product
        return species

    def update_reactions(
        self,
        precursor,
        product,
        component=None,
        part_id=None,
        k=None,
        **kwargs,
    ):
        """Generate metabolic pathway reactions.

        Creates a single irreversible mass-action reaction for the
        metabolic conversion of precursors to products.

        Parameters
        ----------
        precursor : Species, list of Species, or None
            Precursor species or list of precursor species. If None, the
            reaction represents creation (no inputs).
        product : Species, list of Species, or None
            Product species or list of product species. If None, the
            reaction represents degradation (no outputs).
        component : Component, optional
            Component containing parameter values. Required if k is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to using
            component's parameter lookup without specific part_id.
        k : Parameter or float, optional
            Forward rate constant for the pathway. If None, retrieved from
            component parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single irreversible mass-action reaction:
            precursors --> products.

        Raises
        ------
        ValueError
            If component is None and k is not provided.

        Notes
        -----
        The reaction follows mass-action kinetics with rate constant 'k'.
        The mechanism supports three modes:

        - Standard: precursors --> products
        - Creation: [] --> products (when precursor is None)
        - Degradation: precursors --> [] (when product is None)

        Multiple species of the same type in the precursor or product list
        determine the stoichiometry of that species in the reaction.

        """
        if precursor is None:
            inputs = []
        else:
            inputs = precursor

        if product is None:
            outputs = []
        else:
            outputs = product

        if component is None and k is None:
            raise ValueError("Must pass in a component or a rate k.")
        elif k is None:
            k = component.get_parameter('k', part_id=part_id, mechanism=self)

        r = Reaction.from_massaction(inputs, outputs, k_forward=k)
        return [r]
