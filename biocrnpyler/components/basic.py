#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from ..core.component import Component
from ..core.reaction import Reaction
from ..core.species import Complex, Species


class DNA(Component):
    """DNA sequence component with specified length.

    A `DNA` component represents a DNA sequence with a given length in base
    pairs. This is a basic component that produces no reactions on its own
    but can be used as a building block for more complex genetic constructs.

    Parameters
    ----------
    name : str
        Name of the DNA sequence.
    length : int, default=0
        Length of the DNA sequence in base pairs.
    attributes : list of str, optional
        List of attribute tags to associate with the DNA species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    species : Species
        The DNA species object with material_type='dna'.

    See Also
    --------
    RNA : RNA sequence component.
    Protein : Protein sequence component.
    Component : Base class for biomolecular components.

    Notes
    -----
    The DNA component produces no reactions by itself. It is typically used
    as a building block in more complex components like DNA assemblies,
    promoters, or genes that can be transcribed.

    Examples
    --------
    Create a simple DNA sequence:

    >>> dna = bcp.DNA(name='my_gene', length=1000)
    >>> dna.get_species()
    dna_my_gene

    Create DNA with attributes:

    >>> promoter = bcp.DNA(
    ...     name='pLac',
    ...     length=100,
    ...     attributes=['inducible', 'strong']
    ... )

    """

    def __init__(self, name, length=0, attributes=None, **kwargs):
        self.species = self.set_species(
            name, material_type='dna', attributes=attributes
        )
        # self.length = length  # RMM, 14 Sep '25: was not being set
        super().__init__(name=name, **kwargs)

    def get_species(self) -> Species:
        """Get the DNA species.

        Returns
        -------
        Species
            The DNA species object with material_type='dna'.

        """
        return self.species

    def update_species(self) -> List[Species]:
        """Generate species associated with the DNA component.

        Returns
        -------
        list of Species
            List containing only the DNA species itself, as DNA produces no
            additional species.

        """
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        """Generate reactions associated with the DNA component.

        Returns
        -------
        list
            Empty list, as DNA produces no reactions.

        """
        return []


class RNA(Component):
    """RNA sequence component with specified length.

    An `RNA` component represents an RNA sequence with a given length in
    base pairs. This is a basic component that produces no reactions on its
    own but can be used to represent mRNA, tRNA, rRNA, or other RNA
    molecules.

    Parameters
    ----------
    name : str
        Name of the RNA sequence.
    length : int, default=0
        Length of the RNA sequence in base pairs.
    attributes : list of str, optional
        List of attribute tags to associate with the RNA species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    species : Species
        The RNA species object with material_type='rna'.

    See Also
    --------
    DNA : DNA sequence component.
    Protein : Protein sequence component.
    Component : Base class for biomolecular components.

    Notes
    -----
    The RNA component produces no reactions by itself. It is typically used
    to represent transcripts or can be part of more complex components that
    implement translation or RNA degradation mechanisms.

    Examples
    --------
    Create a simple RNA sequence:

    >>> rna = bcp.RNA(name='my_transcript', length=500)
    >>> rna.get_species()
    rna_my_transcript

    Create mRNA with attributes:

    >>> mrna = bcp.RNA(
    ...     name='gfp_mrna',
    ...     length=750,
    ...     attributes=['coding', 'stable']
    ... )

    """

    def __init__(self, name: str, length=0, attributes=None, **kwargs):
        self.species = self.set_species(
            name, material_type='rna', attributes=attributes
        )
        # self.length = length  # RMM, 14 Sep '25: was not being set
        super().__init__(name=name, **kwargs)

    def get_species(self) -> Species:
        """Get the RNA species.

        Returns
        -------
        Species
            The RNA species object with material_type='rna'.

        """
        return self.species

    def update_species(self) -> List[Species]:
        """Generate species associated with the RNA component.

        Returns
        -------
        list of Species
            List containing only the RNA species itself, as RNA produces no
            additional species.

        """
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        """Generate reactions associated with the RNA component.

        Returns
        -------
        list
            Empty list, as RNA produces no reactions.

        """
        return []


class Protein(Component):
    """Protein component with specified length.

    A `Protein` component represents a protein or peptide with a given
    length in amino acids. This is a basic component that produces no
    reactions on its own but can be used to represent enzymes, transcription
    factors, structural proteins, or any other protein molecules.

    Parameters
    ----------
    name : str
        Name of the protein.
    length : int, default=0
        Length of the protein in number of amino acids.
    attributes : list of str, optional
        List of attribute tags to associate with the protein species. Common
        attributes include degradation tags (e.g., 'ssrAtagged') or functional
        properties (e.g., 'fluorescent').
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    species : Species
        The protein species object with material_type='protein'.

    See Also
    --------
    DNA : DNA sequence component.
    RNA : RNA sequence component.
    Enzyme : Enzymatic protein component.
    Component : Base class for biomolecular components.

    Notes
    -----
    The Protein component produces no reactions by itself. It is typically
    used to represent translation products or can be part of more complex
    components like enzymes or binding complexes.

    Examples
    --------
    Create a simple protein:

    >>> protein = bcp.Protein(name='GFP', length=238)
    >>> protein.get_species()
    protein_GFP

    Create a protein with degradation tag:

    >>> protein = bcp.Protein(
    ...     name='LacI',
    ...     length=360,
    ...     attributes=['ssrAtagged']
    ... )

    """

    def __init__(self, name: str, length=0, attributes=None, **kwargs):
        self.species = self.set_species(
            name, material_type='protein', attributes=attributes
        )
        # self.length = length  # RMM, 14 Sep '25: was not being set
        super().__init__(name=name, **kwargs)

    def get_species(self) -> Species:
        """Get the protein species.

        Returns
        -------
        Species
            The protein species object with material_type='protein'.

        """
        return self.species

    def update_species(self) -> List[Species]:
        """Generate species associated with the protein component.

        Returns
        -------
        list of Species
            List containing only the protein species itself, as Protein
            produces no additional species.

        """
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        """Generate reactions associated with the protein component.

        Returns
        -------
        list
            Empty list, as Protein produces no reactions.

        """
        return []


class Metabolite(Component):
    """Metabolic compound that can be produced, utilized, or degraded.

    A `Metabolite` component represents a metabolic compound that
    participates in biochemical pathways. It can have precursors (species
    that are converted into this metabolite) and products (species that this
    metabolite is converted into). The component uses a 'metabolic_pathway'
    mechanism to generate production and degradation reactions.

    Parameters
    ----------
    name : str
        Name of the metabolite.
    attributes : list of str, optional
        List of attribute tags to associate with the metabolite species.
    precursors : list of Species, str, Component, or None, optional
        List of chemical species that are directly transformed into this
        metabolite via the production mechanism. None represents
        constitutive production (production from nothing).
    products : list of Species, str, Component, or None, optional
        List of chemical species produced from this metabolite via the
        degradation mechanism. None represents total degradation
        (degradation to nothing).
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    species : Species
        The metabolite species object with material_type='metabolite'.
    precursors : list of Species or None
        List of precursor species. None values represent constitutive
        production.
    products : list of Species or None
        List of product species. None values represent total degradation.

    See Also
    --------
    Enzyme : Enzymatic component for catalysis.
    Component : Base class for biomolecular components.

    Notes
    -----
    The Metabolite component looks for a 'metabolic_pathway' mechanism but
    will not throw an error if it is not found. If the mechanism is present:

    - Production reactions are generated from precursors to the metabolite
    - Degradation reactions are generated from the metabolite to products

    None is a valid precursor/product representing constitutive
    production/degradation.

    Examples
    --------
    Create a metabolite with constitutive production and degradation:

    >>> atp = bcp.Metabolite(
    ...     name='ATP',
    ...     precursors=[None],
    ...     products=[None]
    ... )

    Create a metabolite with specific precursor and product:

    >>> adp = bcp.Metabolite(
    ...     name='ADP',
    ...     precursors=['ATP'],
    ...     products=['AMP']
    ... )

    Use with a mixture and metabolic pathway mechanism:

    >>> from biocrnpyler.mechanisms import OneStepPathway
    >>> mixture = bcp.Mixture(
    ...     components=[atp],
    ...     mechanisms={'metabolic_pathway': OneStepPathway()},
    ...     parameters={'k': 0.1}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        name: str,
        attributes=None,
        precursors=None,
        products=None,
        **kwargs,
    ):
        self.species = self.set_species(
            name, material_type='metabolite', attributes=attributes
        )

        # Set percursor species list
        self.precursors = []
        if precursors is not None:
            for p in precursors:
                if p is None:
                    # Valid precursor representing constuitive production
                    self.precursors.append(None)
                else:
                    self.precursors.append(self.set_species(p))

        # Set product species list
        self.products = []
        if products is not None:
            for p in products:
                if p is None:
                    # None is a valid product representing total degradation
                    self.products.append(None)
                else:
                    self.products.append(self.set_species(p))

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self) -> Species:
        """Get the metabolite species.

        Returns
        -------
        Species
            The metabolite species object with material_type='metabolite'.

        """
        return self.species

    def update_species(self) -> List[Species]:
        """Generate species for metabolite production and degradation.

        Uses the 'metabolic_pathway' mechanism (if present) to generate
        species for production reactions (from precursors to metabolite) and
        degradation reactions (from metabolite to products).

        Returns
        -------
        list of Species
            List of species including the metabolite itself and any
            additional species generated by the metabolic_pathway mechanism.
            If no mechanism is present, returns only the metabolite species.

        """
        species = [self.get_species()]
        mech_pathway = self.get_mechanism(
            'metabolic_pathway', optional_mechanism=True
        )

        if mech_pathway is not None:
            if len(self.precursors) > 0:
                species += mech_pathway.update_species(
                    precursor=self.precursors,
                    product=[self.get_species()],
                    component=self,
                    part_id=self.name + '_production',
                )
            if len(self.products) > 0:
                species += mech_pathway.update_species(
                    precursor=[self.get_species()],
                    product=self.products,
                    component=self,
                    part_id=self.name + '_degradation',
                )
        return species

    def update_reactions(self) -> List:
        """Generate reactions for metabolite production and degradation.

        Uses the 'metabolic_pathway' mechanism (if present) to generate
        production reactions (from precursors to metabolite) and degradation
        reactions (from metabolite to products).

        Returns
        -------
        list of Reaction
            List of reactions including production and degradation pathways.
            If no mechanism is present, returns an empty list.

        """
        reactions = []
        mech_pathway = self.get_mechanism(
            'metabolic_pathway', optional_mechanism=True
        )

        if mech_pathway is not None:
            if len(self.precursors) > 0:
                reactions += mech_pathway.update_reactions(
                    precursor=self.precursors,
                    product=[self.get_species()],
                    component=self,
                    part_id=self.name + '_production',
                )
            if len(self.products) > 0:
                reactions += mech_pathway.update_reactions(
                    precursor=[self.get_species()],
                    product=self.products,
                    component=self,
                    part_id=self.name + '_degradation',
                )
        return reactions


class ChemicalComplex(Component):
    """Complex formed by binding of two or more molecular species.

    A `ChemicalComplex` component represents a molecular complex formed when
    two or more species bind together. The complex automatically inherits
    attributes from its constituent species. The component uses a 'binding'
    mechanism to generate binding and unbinding reactions.

    Parameters
    ----------
    species : list of Species, str, or Component
        List of species that form the complex. Must contain at least two
        elements. Each element can be a `Species` object, string name, or
        `Component` with an associated species.
    name : str, optional
        Name of the complex. If None, a name is automatically generated
        from the constituent species names.
    material_type : str, default='complex'
        Material type identifier for the complex species. Can be customized
        for specific complex types.
    attributes : list of str, optional
        List of attribute tags to associate with the complex species. The
        complex also inherits attributes from its constituent species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    species : Complex
        The complex species object created from the constituent species.
    internal_species : list of Species
        List of individual species that make up the complex.

    See Also
    --------
    Component : Base class for biomolecular components.
    Species : Chemical species representation.
    Complex : Species subclass for molecular complexes.

    Notes
    -----
    The ChemicalComplex component uses a 'binding' mechanism which must be
    provided by the containing mixture. The binding mechanism generates:

    - Forward binding reactions (species --> complex)
    - Reverse unbinding reactions (complex --> species)

    The first species in the list is treated as the 'bindee' and remaining
    species are treated as 'binders' in the binding mechanism.

    Examples
    --------
    Create a simple protein-DNA complex:

    >>> complex = bcp.ChemicalComplex(
    ...     species=['TF_protein', 'DNA_promoter'],
    ...     name='TF_bound'
    ... )

    Create an enzyme-substrate complex:

    >>> complex = bcp.ChemicalComplex(
    ...     species=['protein_E', 'S'],
    ...     name='ES_complex'
    ... )

    Use with a mixture and binding mechanism:

    >>> from biocrnpyler.mechanisms import One_Step_Binding
    >>> mixture = bcp.Mixture(
    ...     components=[complex],
    ...     mechanisms={'binding': One_Step_Binding()},
    ...     parameters={'kb': 1.0, 'ku': 0.1}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        species: List[Species],
        name: str = None,
        material_type='complex',
        attributes=None,
        **kwargs,
    ):
        if not isinstance(species, list) or len(species) < 2:
            raise ValueError(
                f"Invalid Species {species}. Species must be a list of "
                "Species, strings, or Component objects."
            )

        self.internal_species = []  # a list of species inside the complex

        for s in species:
            self.internal_species.append(self.set_species(s))
        if attributes is None:
            attributes = []
        self.species = Complex(
            species=self.internal_species,
            name=name,
            material_type=material_type,
            attributes=attributes,
        )

        if name is None:
            name = self.species.name

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self) -> List[Species]:
        """Get the complex species.

        Returns
        -------
        Complex
            The complex species object containing all constituent species.

        """
        return self.species

    def update_species(self) -> List[Species]:
        """Generate species for complex binding reactions.

        Uses the 'binding' mechanism to generate all species needed for
        binding and unbinding reactions, including the individual species
        and the complex.

        Returns
        -------
        list of Species
            List of all species generated by the binding mechanism,
            typically including the constituent species and the complex
            species.

        """
        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        species = mech_b.update_species(
            binder,
            bindee,
            complex_species=self.get_species(),
            component=self,
            part_id=self.name,
        )
        return species

    def update_reactions(self) -> List[Reaction]:
        """Generate binding and unbinding reactions for the complex.

        Uses the 'binding' mechanism to generate reactions for complex
        formation (binding) and dissociation (unbinding).

        Returns
        -------
        list of Reaction
            List of reactions generated by the binding mechanism, typically
            including forward binding and reverse unbinding reactions.

        """
        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        reactions = mech_b.update_reactions(
            binder,
            bindee,
            complex_species=self.get_species(),
            component=self,
            part_id=self.name,
        )
        return reactions


class Enzyme(Component):
    """Enzyme that catalyzes conversion of substrates to products.

    An `Enzyme` component represents an enzyme that catalyzes the conversion
    of one or more substrates into one or more products. The enzyme itself
    is not consumed in the reaction. This component uses a 'catalysis'
    mechanism to generate the appropriate chemical reactions.

    Parameters
    ----------
    enzyme : Species, str, or Component
        The enzyme species that catalyzes the reaction. Can be a `Species`
        object, a string name (creates new protein Species), or a
        `Component` with an associated species.
    substrates : list of Species, str, or Component
        List of substrate species that are consumed by the enzymatic
        reaction. Each element can be a `Species` object, string name, or
        `Component`.
    products : list of Species, str, or Component
        List of product species that are produced by the enzymatic
        reaction. Each element can be a `Species` object, string name, or
        `Component`.
    attributes : list of str, optional
        List of attribute tags to associate with the enzyme species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    enzyme : Species
        The enzyme species object.
    substrates : list of Species
        List of substrate species objects.
    products : list of Species
        List of product species objects.

    See Also
    --------
    Component : Base class for biomolecular components.
    Metabolite : Component for metabolic compounds.
    ChemicalComplex : Component for molecular complexes.

    Notes
    -----
    The `Enzyme` component assumes all substrates are converted to all
    products in a single enzymatic step:

        S1 + S2 + ... + SN + E --> P1 + P2 + ... + PM + E

    For enzymes that catalyze multiple distinct reactions, create separate
    `Enzyme` components with the same internal enzyme species.

    The component uses a mechanism called 'catalysis' which must be
    provided by the containing mixture. Common catalysis mechanisms include
    Michaelis-Menten kinetics and other enzymatic rate laws.

    Examples
    --------
    Create a simple enzyme that converts substrate S to product P:

    >>> enzyme = bcp.Enzyme(
    ...     enzyme='E',
    ...     substrates=['S'],
    ...     products=['P']
    ... )
    >>> enzyme.get_species()
    protein_E

    Create an enzyme with multiple substrates and products:

    >>> enzyme = bcp.Enzyme(
    ...     enzyme='Kinase',
    ...     substrates=['ATP', 'Protein'],
    ...     products=['ADP', 'Protein_P']
    ... )

    Use with a mixture and Michaelis-Menten mechanism:

    >>> from biocrnpyler.mechanisms import MichaelisMenten
    >>> mixture = bcp.Mixture(
    ...     components=[enzyme],
    ...     mechanisms={'catalysis': MichaelisMenten()},
    ...     parameters={'kb': 0.1, 'ku': 0.01, 'kcat': 1.0}
    ... )
    >>> crn = mixture.compile_crn()

    """

    # TODO: implement multiple substrates and multiple products

    def __init__(
        self,
        enzyme: Union[Species, str, Component],
        substrates: List[Union[Species, str, Component]],
        products: List[Union[Species, str, Component]],
        attributes=None,
        **kwargs,
    ):
        self.enzyme = self.set_species(
            enzyme, material_type='protein', attributes=attributes
        )
        self.substrates = substrates
        self.products = products

        Component.__init__(self=self, name=self.enzyme.name, **kwargs)

    @property
    def substrates(self) -> List:
        """List of substrate species for the enzymatic reaction.

        Returns
        -------
        list of Species

        """
        return self._substrates

    @substrates.setter
    def substrates(
        self, new_substrates: List[Union[Species, str, Component]]
    ):
        """Set the substrate species list.

        Parameters
        ----------
        new_substrates : Species, str, Component, or list
            Substrate(s) to set. Can be a single species or a list. Each
            element is converted to a `Species` object.

        Notes
        -----
        Automatically converts single substrate to a list and converts all
        elements to `Species` objects using `set_species`.

        """
        if not isinstance(new_substrates, list):
            new_substrates = [new_substrates]
        # convert the new substrates to Species
        self._substrates = [self.set_species(s) for s in new_substrates]

    @property
    def products(self) -> List:
        """List of product species for the enzymatic reaction.

        Returns
        -------
        list of Species

        """
        return self._products

    @products.setter
    def products(self, new_products: List[Union[Species, str, Component]]):
        """Set the product species list.

        Parameters
        ----------
        new_products : Species, str, Component, or list
            Product(s) to set. Can be a single species or a list. Each
            element is converted to a `Species` object.

        Notes
        -----
        Automatically converts single product to a list and converts all
        elements to `Species` objects using `set_species`.

        """
        if not isinstance(new_products, list):
            new_products = [new_products]
        # convert the new products to Products
        self._products = [self.set_species(p) for p in new_products]

    def get_species(self) -> Species:
        """Get the enzyme species.

        Returns
        -------
        Species
            The enzyme species object that catalyzes the reaction.

        """
        return self.enzyme

    def update_species(self) -> List[Species]:
        """Generate species required for enzymatic catalysis.

        Uses the 'catalysis' mechanism to generate all species needed for
        the enzymatic reaction, including enzyme, substrates, products, and
        any intermediate complexes.

        Returns
        -------
        list of Species
            List of all species generated by the catalysis mechanism,
            typically including enzyme, substrates, products, and
            enzyme-substrate complexes.

        """
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_species(
            enzyme=self.enzyme,
            substrate=self.substrates,
            product=self.products,
        )

    def update_reactions(self) -> List[Reaction]:
        """Generate reactions for enzymatic catalysis.

        Uses the 'catalysis' mechanism to generate all reactions needed for
        the enzymatic conversion of substrates to products.

        Returns
        -------
        list of Reaction
            List of all reactions generated by the catalysis mechanism,
            typically including substrate binding, catalysis, and product
            release steps.

        """
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_reactions(
            enzyme=self.enzyme,
            substrate=self.substrates,
            product=self.products,
            component=self,
            part_id=self.name,
        )
