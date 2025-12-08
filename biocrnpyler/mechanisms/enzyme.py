# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex


class BasicCatalysis(Mechanism):
    r"""Basic catalytic mechanism for irreversible substrate conversion.

    A 'catalysis' mechanism where a catalyst (enzyme) converts a substrate
    into a product in a single irreversible step. The catalyst is not
    consumed in the reaction and can continue to catalyze additional
    conversions.

    The catalytic reaction is given by

        $$ 'S' + 'C' --> 'P' + 'C' $$

    where S is the substrate, C is the catalyst (enzyme), and P is the
    product.

    Parameters
    ----------
    name : str, default='basic_catalysis'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='catalysis'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('catalysis').

    See Also
    --------
    BasicProduction : Catalytic production without substrate consumption.
    MichaelisMenten : Two-step enzyme kinetics with complex formation.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates a single irreversible mass-action reaction
    with rate constant 'kcat'. Unlike Michaelis-Menten kinetics, there is
    no explicit enzyme-substrate complex formation; the reaction proceeds
    in a single catalytic step.

    Common applications include:

    - Simplified enzyme kinetics models
    - Catalytic degradation reactions
    - Rate-limiting steps in metabolic pathways

    Required parameters for this mechanism:

    - 'kcat' : Catalytic rate constant for substrate conversion

    Examples
    --------
    Model enzymatic degradation of a substrate:

    >>> enzyme = bcp.Enzyme('E', substrates=['S'], products=['P'])
    >>> mixture = bcp.Mixture(
    ...     components=[enzyme],
    ...     mechanisms={'catalysis': bcp.BasicCatalysis()},
    ...     parameters={'kcat': 1.0}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name: str = 'basic_catalysis', mechanism_type: str = 'catalysis'
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, enzyme, substrate, product=None):
        r"""Generate species for basic catalysis.

        Creates the list of species involved in the catalytic reaction:
        enzyme, substrate, and optionally the product.

        Parameters
        ----------
        enzyme : Species
            The catalyst species that facilitates the reaction.
        substrate : Species
            The substrate species to be converted.
        product : Species, optional
            The product species. If None, only enzyme and substrate are
            returned (useful for degradation reactions where no explicit
            product is tracked).

        Returns
        -------
        list of Species
            List containing [enzyme, substrate] if product is None, or
            [enzyme, substrate, product] otherwise.

        """
        if product is None:
            return [enzyme, substrate]
        else:
            return [enzyme, substrate, product]

    def update_reactions(
        self,
        enzyme,
        substrate,
        product,
        component=None,
        part_id=None,
        kcat=None,
    ):
        """Generate reactions for basic catalysis.

        Creates a single irreversible mass-action reaction for catalytic
        conversion of substrate to product.

        Parameters
        ----------
        enzyme : Species
            The catalyst species that facilitates the reaction.
        substrate : Species
            The substrate species to be converted.
        product : Species
            The product species. Can be None for degradation reactions.
        component : Component, optional
            Component containing parameter values. Required if kcat is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        kcat : Parameter or float, optional
            Catalytic rate constant. If None, retrieved from component
            parameters.

        Returns
        -------
        list of Reaction
            List containing a single irreversible mass-action reaction:
            enzyme + substrate --> enzyme + product.

        Raises
        ------
        ValueError
            If component is None and kcat is not provided.

        Notes
        -----
        The reaction follows mass-action kinetics with rate constant 'kcat'.
        The enzyme appears on both sides of the reaction as it acts as a
        catalyst and is not consumed.

        """
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self
            )

        if product is None:
            product = []

        return [
            Reaction.from_massaction(
                inputs=[enzyme, substrate],
                outputs=[enzyme, product],
                k_forward=kcat,
            )
        ]


class BasicProduction(Mechanism):
    r"""Basic catalytic production mechanism with optional substrate.

    A 'catalysis' mechanism where a catalyst (enzyme) produces a product.
    Optionally, a substrate can be consumed during production, allowing for
    both pure production (C --> P + C) and production with substrate
    consumption (S + C --> P + C).

    The production reaction can be either:
    $$
        'C' --> 'P' + 'C' \quad\text{(pure production, no substrate)}
    $$
    or
    $$
        'S' + 'C' --> 'P' + 'C'
             \quad\text{(production with substrate consumption)}
    $$
    where S is the substrate, C is the catalyst (enzyme), and P is the
    product.

    Parameters
    ----------
    name : str, default='basic_production'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='catalysis'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('catalysis').

    See Also
    --------
    BasicCatalysis : Catalytic conversion requiring a substrate.
    MichaelisMentenCopy : Two-step kinetics preserving the substrate.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates a single irreversible mass-action reaction
    with rate constant 'kcat'. The catalyst is not consumed and appears on
    both sides of the reaction.

    Common applications include:

    - Constitutive gene expression (transcription/translation)
    - Enzymatic synthesis reactions
    - Autocatalytic production processes

    Required parameters for this mechanism:

    - 'kcat' : Catalytic rate constant for product formation

    The flexibility to include or exclude substrates makes this mechanism
    useful for modeling both simple production (e.g., constitutive protein
    expression) and production coupled with substrate consumption (e.g.,
    enzymatic synthesis from precursors).

    Examples
    --------
    Model constitutive protein production from a gene:

    >>> gene = bcp.DNA('gfp')
    >>> protein = bcp.Protein('GFP')
    >>> expression = bcp.Enzyme(gene, substrates=[], products=[protein])
    >>> mixture = bcp.Mixture(
    ...     components=[expression],
    ...     mechanisms={'catalysis': bcp.BasicProduction()},
    ...     parameters={'kcat': 0.01}
    ... )
    >>> mixture.compile_crn()
    Species = dna_gfp, protein_GFP
    Reactions = [
        dna[gene] --> dna[gene]+protein[protein]
    ]

    """

    def __init__(self, name='basic_production', mechanism_type='catalysis'):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, enzyme, substrate=None, product=None):
        """Generate species for basic production.

        Creates the list of species involved in the production reaction:
        enzyme, and optionally substrate and product.

        Parameters
        ----------
        enzyme : Species
            The catalyst species that facilitates production.
        substrate : Species, optional
            The substrate species to be consumed. If None, production
            occurs without substrate consumption.
        product : Species, optional
            The product species. If None, only enzyme (and substrate if
            provided) are returned.

        Returns
        -------
        list of Species
            List containing enzyme and any non-None substrate and product
            species. Order is [enzyme, product, substrate] if all are
            provided.

        """
        species = [enzyme]
        if product is not None:
            species += [product]
        if substrate is not None:
            species += [substrate]

        return species

    def update_reactions(
        self,
        enzyme,
        substrate,
        product,
        component=None,
        part_id=None,
        kcat=None,
    ):
        """Generate reactions for basic production.

        Creates a single irreversible mass-action reaction for catalytic
        production, with or without substrate consumption.

        Parameters
        ----------
        enzyme : Species
            The catalyst species that facilitates production.
        substrate : Species
            The substrate species. Can be None for pure production without
            substrate consumption.
        product : Species
            The product species. Can be None if no explicit product is
            tracked.
        component : Component, optional
            Component containing parameter values. Required if kcat is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        kcat : Parameter or float, optional
            Catalytic rate constant. If None, retrieved from component
            parameters.

        Returns
        -------
        list of Reaction
            List containing a single irreversible mass-action reaction.
            If substrate is None: enzyme --> enzyme + product.
            If substrate is provided: enzyme + substrate --> enzyme +
            product.

        Raises
        ------
        ValueError
            If component is None and kcat is not provided.

        Notes
        -----
        The enzyme appears on both sides of the reaction as it acts as a
        catalyst and is not consumed. The substrate, if provided, is
        consumed in the reaction.

        """
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self
            )

        inputs = [enzyme]
        outputs = [enzyme]
        if product is not None:
            outputs += [product]
        if substrate is not None:
            inputs += [substrate]

        return [
            Reaction.from_massaction(
                inputs=inputs, outputs=outputs, k_forward=kcat
            )
        ]


class MichaelisMenten(Mechanism):
    r"""Standard Michaelis-Menten enzyme kinetics mechanism.

    A 'catalysis' mechanism implementing classical Michaelis-Menten enzyme
    kinetics with explicit enzyme-substrate complex formation. The substrate
    binds reversibly to the enzyme to form a complex, which then
    irreversibly converts to product and releases the enzyme.

    The reaction scheme is

        $$ 'S' + 'E' <--> 'S':'E' --> 'E' + 'P' $$

    where S is the substrate, E is the enzyme, S:E is the enzyme-substrate
    complex, and P is the product.

    Parameters
    ----------
    name : str, default='michaelis_menten'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='catalysis'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('catalysis').

    See Also
    --------
    BasicCatalysis : Single-step catalysis without complex formation.
    MichaelisMentenCopy : Michaelis-Menten preserving substrate.
    MichaelisMentenReversible : Michaelis-Menten with product binding.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates two mass-action reactions:

    1. Reversible binding: S + E <--> S:E (rates 'kb' and 'ku')
    2. Irreversible catalysis: S:E --> E + P (rate 'kcat')

    Common applications include:

    - Enzyme-catalyzed reactions in metabolic pathways
    - Protein degradation by proteases
    - Drug metabolism by cytochrome P450 enzymes
    - Any enzymatic process following Michaelis-Menten kinetics

    Required parameters for this mechanism:

    - 'kb' : Binding rate constant for enzyme-substrate association
    - 'ku' : Unbinding rate constant for enzyme-substrate dissociation
    - 'kcat' : Catalytic rate constant for product formation

    The mechanism can also model degradation reactions by setting product
    to None, resulting in: S + E <--> S:E --> E.

    Examples
    --------
    Model enzyme-catalyzed substrate conversion:

    >>> substrate = bcp.Species('S')
    >>> product = bcp.Species('P')
    >>> enzyme = bcp.Enzyme('E', substrates=[substrate], products=[product])
    >>> mixture = bcp.Mixture(
    ...     components=[enzyme],
    ...     mechanisms={'catalysis': bcp.MichaelisMenten()},
    ...     parameters={'kb': 1.0, 'ku': 0.1, 'kcat': 0.5}
    ... )
    >>> mixture.compile_crn()
    Species = protein_E, S, P, complex_S_protein_E_
    Reactions = [
       S+protein[E] <--> complex[S:protein[E]]
       complex[S:protein[E]] --> P+protein[E]
    ]

    Model enzymatic degradation:

    >>> degradase = bcp.Protein('degradase')
    >>> target = bcp.Protein('target')
    >>> degrader = bcp.Enzyme(degradase, substrates=[target], products=[])
    >>> mixture = bcp.Mixture(
    ...     components=[degrader],
    ...     mechanisms={'catalysis': bcp.MichaelisMenten()},
    ...     parameters={'kb': 1.0, 'ku': 0.1, 'kcat': 0.2}
    ... )

    """

    def __init__(self, name='michaelis_menten', mechanism_type='catalysis'):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, enzyme, substrate, product=None, complex=None):
        """Generate species for Michaelis-Menten kinetics.

        Creates the species involved in Michaelis-Menten enzyme kinetics:
        enzyme, substrate, enzyme-substrate complex, and optionally the
        product.

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the reaction.
        substrate : Species
            The substrate species to be converted.
        product : Species, optional
            The product species. If None, only enzyme, substrate, and
            complex are returned (useful for degradation reactions).
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).

        Returns
        -------
        list of Species
            List containing [enzyme, substrate, complex] if product is
            None, or [enzyme, substrate, product, complex] otherwise.

        Notes
        -----
        The complex is automatically generated as a Complex object
        containing the substrate and enzyme if not explicitly provided.

        """
        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex
        if product is None:
            return [enzyme, substrate, complexS]
        else:
            return [enzyme, substrate, product, complexS]

    def update_reactions(
        self,
        enzyme,
        substrate,
        product,
        component=None,
        part_id=None,
        complex=None,
        kb=None,
        ku=None,
        kcat=None,
    ):
        r"""Generate reactions for Michaelis-Menten kinetics.

        Creates two mass-action reactions implementing Michaelis-Menten
        enzyme kinetics: reversible enzyme-substrate binding and
        irreversible catalytic conversion.

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the reaction.
        substrate : Species
            The substrate species to be converted.
        product : Species
            The product species. Can be None for degradation reactions.
        component : Component, optional
            Component containing parameter values. Required if kb, ku, or
            kcat are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).
        kb : Parameter or float, optional
            Forward binding rate constant. If None, retrieved from
            component parameters.
        ku : Parameter or float, optional
            Reverse unbinding rate constant. If None, retrieved from
            component parameters.
        kcat : Parameter or float, optional
            Catalytic rate constant. If None, retrieved from component
            parameters.

        Returns
        -------
        list of Reaction
            List containing two reactions:
            [binding_reaction, catalysis_reaction].

        Raises
        ------
        ValueError
            If component is None and any of kb, ku, or kcat is not
            provided.

        Notes
        -----
        The mechanism generates the following reactions:

        1. S + E <--> S:E (binding, rates 'kb' and 'ku')
        2. S:E --> E + P (catalysis, rate 'kcat')

        For degradation (product is None):

        2. S:E --> E (degradation, rate 'kcat')

        """
        # Get parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat."
            )
        if kb is None:
            kb = component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            )
        if ku is None:
            ku = component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            )
        if kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self
            )

        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex

        # substrate + Enz <--> substrate:Enz
        binding_rxn = Reaction.from_massaction(
            inputs=[substrate, enzyme],
            outputs=[complexS],
            k_forward=kb,
            k_reverse=ku,
        )
        if product is not None:
            # substrate:Enz --> Enz + product
            cat_rxn = Reaction.from_massaction(
                inputs=[complexS], outputs=[product, enzyme], k_forward=kcat
            )
        else:  # degradation Reaction
            # substrate:Enz --> Enz
            cat_rxn = Reaction.from_massaction(
                inputs=[complexS], outputs=[enzyme], k_forward=kcat
            )
        return [binding_rxn, cat_rxn]


class MichaelisMentenReversible(Mechanism):
    r"""Reversible Michaelis-Menten kinetics with product binding.

    A 'catalysis' mechanism implementing Michaelis-Menten enzyme kinetics
    where the product can also bind reversibly to the enzyme. Both the
    substrate and product form distinct enzyme complexes, and the catalytic
    step itself is reversible.

    The reaction scheme is

        $$ 'S' + 'E' <--> 'S':'E' <--> 'E':'P' <--> 'E' + 'P', $$

    where S is the substrate, E is the enzyme, S:E is the enzyme-substrate
    complex, E:P is the enzyme-product complex, and P is the product.

    Parameters
    ----------
    name : str, default='michaelis_menten_reverse_binding'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='catalysis'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('catalysis').

    See Also
    --------
    MichaelisMenten : Standard Michaelis-Menten with irreversible catalysis.
    MichaelisMentenCopy : Michaelis-Menten preserving substrate.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates three mass-action reactions:

    1. Reversible substrate binding: S + E <--> S:E (rates 'kb1' and 'ku1')

    2. Reversible product binding: P + E <--> E:P (rates 'kb2' and 'ku2')

    3. Reversible catalysis: S:E <--> E:P (rates 'kcat' and 'kcat_rev')

    Common applications include:

    - Reversible enzymatic reactions near equilibrium
    - Bidirectional metabolic pathways
    - Reactions where product inhibition is significant
    - Detailed kinetic models requiring thermodynamic consistency

    Required parameters for this mechanism:

    - 'kb1' : Forward binding rate for substrate-enzyme association
    - 'ku1' : Reverse unbinding rate for substrate-enzyme dissociation
    - 'kb2' : Forward binding rate for product-enzyme association
    - 'ku2' : Reverse unbinding rate for product-enzyme dissociation
    - 'kcat' : Forward catalytic rate constant (S:E --> E:P)
    - 'kcat_rev' : Reverse catalytic rate constant (E:P --> S:E)

    This mechanism is particularly useful when modeling reactions close to
    equilibrium where the reverse reaction and product binding cannot be
    neglected.

    Examples
    --------
    Model a reversible enzymatic conversion:

    >>> enzyme = bcp.Species('E', material_type='protein')
    >>> substrate = bcp.Species('S')
    >>> product = bcp.Species('P')
    >>> comp = bcp.Enzyme(
    ...     enzyme, substrates=[substrate], products=[product],
    ...     mechanisms={'catalysis': bcp.MichaelisMentenReversible()},
    ...     parameters={
    ...         'kb1': 2.0, 'ku1': 0.5,
    ...         'kb2': 1.5, 'ku2': 0.3,
    ...         'kcat': 1.0, 'kcat_rev': 0.4
    ...     }
    ... )
    >>> mixture = bcp.Mixture(components=[comp])
    >>> mixture.compile_crn()
    Species = protein_E, S, P, complex_S_protein_E_, complex_P_protein_E_
    Reactions = [
        S+protein[E] <--> complex[S:protein[E]]
        P+protein[E] <--> complex[P:protein[E]]
        complex[S:protein[E]] <--> complex[P:protein[E]]
    ]

    """

    def __init__(
        self,
        name='michaelis_menten_reverse_binding',
        mechanism_type='catalysis',
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self, enzyme, substrate, product, complex=None, complex2=None
    ):
        """Generate species for reversible Michaelis-Menten kinetics.

        Creates the species involved in reversible Michaelis-Menten enzyme
        kinetics: enzyme, substrate, product, enzyme-substrate complex, and
        enzyme-product complex.

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the reaction.
        substrate : Species
            The substrate species.
        product : Species
            The product species.
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).
        complex2 : Species, optional
            Pre-specified enzyme-product complex. If None, automatically
            creates a Complex([product, enzyme]).

        Returns
        -------
        list of Species
            List containing [enzyme, substrate, product, complex1,
            complex2] where complex1 is S:E and complex2 is E:P.

        Notes
        -----
        Both complexes are automatically generated if not explicitly
        provided. The enzyme-substrate complex contains [substrate, enzyme]
        and the enzyme-product complex contains [product, enzyme].

        """
        if complex is None:
            complex1 = Complex([substrate, enzyme])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([product, enzyme])
        else:
            complex2 = complex2
        return [enzyme, substrate, product, complex1, complex2]

    def update_reactions(
        self,
        enzyme,
        substrate,
        product,
        component=None,
        part_id=None,
        complex=None,
        complex2=None,
        kb=None,
        ku=None,
        kcat=None,
    ):
        r"""Generate reactions for reversible Michaelis-Menten kinetics.

        Creates three mass-action reactions implementing reversible
        Michaelis-Menten enzyme kinetics with product binding: substrate
        binding, product binding, and reversible catalysis.

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the reaction.
        substrate : Species
            The substrate species.
        product : Species
            The product species.
        component : Component, optional
            Component containing parameter values. Required if kb, ku, or
            kcat are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).
        complex2 : Species, optional
            Pre-specified enzyme-product complex. If None, automatically
            creates a Complex([product, enzyme]).
        kb : tuple of (float or Parameter), optional
            Tuple of (kb1, kb2) binding rate constants. If None, kb1 and
            kb2 retrieved separately from component parameters.
        ku : tuple of (float or Parameter), optional
            Tuple of (ku1, ku2) unbinding rate constants. If None, ku1 and
            ku2 retrieved separately from component parameters.
        kcat : tuple of (float or Parameter), optional
            Tuple of (kcat, kcat_rev) catalytic rate constants. If None,
            kcat and kcat_rev retrieved separately from component
            parameters.

        Returns
        -------
        list of Reaction
            List containing three reactions: [substrate_binding_reaction,
            product_binding_reaction, catalysis_reaction].

        Raises
        ------
        ValueError
            If component is None and any of kb, ku, or kcat is not
            provided.

        Notes
        -----
        The mechanism generates the following reactions:

        1. S + E <--> S:E (binding, rates 'kb1' and 'ku1')
        2. P + E <--> E:P (binding, rates 'kb2' and 'ku2')
        3. S:E <--> E:P (catalysis, rates 'kcat' and 'kcat_rev')

        When providing parameters directly (not via component), kb, ku, and
        kcat should be tuples of two values each.

        """
        # Get parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat."
            )
        if kb is None:
            kb1 = component.get_parameter(
                'kb1', part_id=part_id, mechanism=self
            )
            kb2 = component.get_parameter(
                'kb2', part_id=part_id, mechanism=self
            )
        else:
            kb1, kb2 = kb
        if ku is None:
            ku1 = component.get_parameter(
                'ku1', part_id=part_id, mechanism=self
            )
            ku2 = component.get_parameter(
                'ku2', part_id=part_id, mechanism=self
            )
        else:
            ku1, ku2 = ku
        if kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self
            )
            kcat_rev = component.get_parameter(
                'kcat_rev', part_id=part_id, mechanism=self
            )
        else:
            kcat, kcat_rev = kcat

        if complex is None:
            complex1 = Complex([substrate, enzyme])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([product, enzyme])

        # substrate + Enz <--> substrate:Enz
        binding_rxn1 = Reaction.from_massaction(
            inputs=[substrate, enzyme],
            outputs=[complex1],
            k_forward=kb1,
            k_reverse=ku1,
        )

        binding_rxn2 = Reaction.from_massaction(
            inputs=[product, enzyme],
            outputs=[complex2],
            k_forward=kb2,
            k_reverse=ku2,
        )

        # substrate:Enz --> Enz:product
        cat_rxn = Reaction.from_massaction(
            inputs=[complex1],
            outputs=[complex2],
            k_forward=kcat,
            k_reverse=kcat_rev,
        )

        return [binding_rxn1, binding_rxn2, cat_rxn]


class MichaelisMentenCopy(Mechanism):
    r"""Michaelis-Menten kinetics with substrate preservation.

    A 'copy' mechanism implementing Michaelis-Menten enzyme kinetics where
    the substrate is not consumed during the reaction. Instead, the
    substrate acts as a template that is copied or read, producing a
    product while preserving the original substrate.

    The reaction scheme is

        $$ 'S' + 'E' <--> 'S':'E' --> 'S' + 'E' + 'P' $$

    where S is the substrate (template), E is the enzyme, S:E is the
    enzyme-substrate complex, and P is the product.

    Parameters
    ----------
    name : str, default='michaelis_menten_copy'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='copy'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('copy').

    See Also
    --------
    MichaelisMenten : Standard Michaelis-Menten consuming substrate.
    BasicProduction : Simpler production without complex formation.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism generates two mass-action reactions:

    1. Reversible binding: S + E <--> S:E (rates 'kb' and 'ku')
    2. Catalytic copying: S:E --> S + E + P (rate 'kcat')

    Common applications include:

    - Gene transcription (DNA template produces RNA)
    - Translation (mRNA template produces protein)
    - DNA replication
    - Any process where a template is read without being consumed

    Required parameters for this mechanism:

    - 'kb' : Binding rate constant for enzyme-substrate association
    - 'ku' : Unbinding rate constant for enzyme-substrate dissociation
    - 'kcat' : Catalytic rate constant for product formation

    The key difference from standard Michaelis-Menten is that the substrate
    appears on both sides of the catalytic step, making it a true copying
    or templating mechanism rather than a conversion.

    Examples
    --------
    Model translation with component:

    >>> mrna = bcp.Species('mRNA')
    >>> ribosome = bcp.Species('Ribo')
    >>> protein = bcp.species('GFP')
    >>> comp = bcp.Enzyme(
    ...     ribosome, substrates=[mrna], products=[protein],
    ...     parameters={'kb': 2.0, 'ku': 0.2, 'kcat': 0.1}
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[comp],
    ...     mechanisms={'catalysis': bcp.MichaelisMentenCopy()},
    ... )
    >>> mixture.compile_crn()
    Species = Ribo, mRNA, complex_Ribo_mRNA_, GFP
    Reactions = [
        mRNA+Ribo <--> complex[Ribo:mRNA]
        complex[Ribo:mRNA] --> mRNA+GFP+Ribo
    ]

    """

    def __init__(self, name='michaelis_menten_copy', mechanism_type='copy'):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, enzyme, substrate, complex=None, product=None):
        """Generate species for copy-type Michaelis-Menten kinetics.

        Creates the species involved in copy-type Michaelis-Menten enzyme
        kinetics: enzyme, substrate (template), enzyme-substrate complex,
        and optionally the product(s).

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the copying reaction.
        substrate : Species
            The substrate (template) species that is copied but not
            consumed.
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).
        product : Species or list of Species, optional
            The product species or list of products. If None, only enzyme,
            substrate, and complex are returned.

        Returns
        -------
        list of Species
            List containing [enzyme, substrate, complex] if product is
            None. If product is provided, returns [enzyme, substrate,
            complex, product] for single product or [enzyme, substrate,
            complex] + product for list of products.

        Notes
        -----
        This method can handle multiple products by accepting product as a
        list. This is useful for modeling processes like transcription
        where multiple transcript copies may be produced.

        """
        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex

        if product is None:
            return [enzyme, substrate, complexS]
        elif isinstance(product, list):
            return [enzyme, substrate, complexS] + product
        else:
            return [enzyme, substrate, product, complexS]

    def update_reactions(
        self,
        enzyme,
        substrate,
        product,
        component=None,
        part_id=None,
        complex=None,
        kb=None,
        ku=None,
        kcat=None,
    ):
        r"""Generate reactions for copy-type Michaelis-Menten kinetics.

        Creates two mass-action reactions implementing copy-type
        Michaelis-Menten enzyme kinetics: reversible enzyme-substrate
        binding and catalytic copying that preserves the substrate.

        Parameters
        ----------
        enzyme : Species
            The enzyme species that catalyzes the copying reaction.
        substrate : Species
            The substrate (template) species that is copied but not
            consumed.
        product : Species
            The product species.
        component : Component, optional
            Component containing parameter values. Required if kb, ku, or
            kcat are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Species, optional
            Pre-specified enzyme-substrate complex. If None, automatically
            creates a Complex([substrate, enzyme]).
        kb : Parameter or float, optional
            Forward binding rate constant. If None, retrieved from
            component parameters.
        ku : Parameter or float, optional
            Reverse unbinding rate constant. If None, retrieved from
            component parameters.
        kcat : Parameter or float, optional
            Catalytic rate constant. If None, retrieved from component
            parameters.

        Returns
        -------
        list of Reaction
            List containing two reactions:
            [binding_reaction, catalysis_reaction].

        Raises
        ------
        ValueError
            If component is None and any of kb, ku, or kcat is not
            provided.

        Notes
        -----
        The mechanism generates the following reactions:

        1. S + E <--> S:E (binding, rates 'kb' and 'ku')
        2. S:E --> S + E + P (copying, rate 'kcat')

        The key feature is that the substrate appears on both sides of the
        catalytic reaction, ensuring it is not consumed. This makes the
        reaction a true template-based copying mechanism.

        """
        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex

        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if kb is None and component is not None:
            kb = component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            )
        if ku is None and component is not None:
            ku = component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            )
        if kcat is None and component is not None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self
            )
        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat."
            )
        # substrate + Enz <--> substrate:Enz
        binding_rxn = Reaction.from_massaction(
            inputs=[substrate, enzyme],
            outputs=[complexS],
            k_forward=kb,
            k_reverse=ku,
        )

        # substrate:Enz --> Enz + product + substrate
        cat_rxn = Reaction.from_massaction(
            inputs=[complexS],
            outputs=[substrate, product, enzyme],
            k_forward=kcat,
        )

        return [binding_rxn, cat_rxn]
