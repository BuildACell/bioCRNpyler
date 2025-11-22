#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from ...core.component import Component
from ...core.mechanism import Mechanism
from ...core.mixture import Mixture
from ...core.parameter import ParameterDatabase
from ...core.reaction import Reaction
from ...core.species import Species
from ..basic import DNA, RNA, Protein
from .promoter import Promoter
from .rbs import RBS


class DNAassembly(DNA):
    """High-level representation of a gene expression construct.

    A DNAassembly represents a complete gene expression unit combining a
    promoter region, ribosome binding site (RBS), coding sequence, and the RNA
    and protein products. This class provides a convenient interface for
    modeling the central dogma pathway: DNA --> RNA --> Protein, where
    the promoter controls transcription and the RBS controls translation.

    Parameters
    ----------
    name : str
        Name of the DNA assembly.
    dna : DNA, str, or None, optional
        The DNA species or name for the assembly. If None, a DNA species
        with `name` is created automatically.
    promoter : Promoter, str, or None, optional
        The promoter component or name controlling transcription. If None,
        no transcription occurs. If a string, a default Promoter is created.
    transcript : RNA, str, bool, or None, optional
        The RNA transcript produced by transcription. If None, an RNA
        species with `name` is created. If False, no transcript is created
        (used in expression mixtures for direct translation).
    rbs : RBS, str, or None, optional
        The ribosome binding site component or name controlling translation.
        If None, no translation occurs. If a string, a default RBS is
        created.
    protein : Protein, str, or None, optional
        The protein product of translation. If None, a Protein species with
        `name` is created automatically.
    length : int, optional
        Length of the DNA sequence in base pairs.
    attributes : list of str, optional
        List of attribute tags for the assembly and its species.
    mechanisms : dict or list, optional
        Custom mechanisms for this assembly, overriding mixture defaults.
    compartment : Compartment, optional
        The compartment containing this assembly and its products.
    parameters : dict, optional
        Parameter values specific to this assembly.
    initial_concentration : float, optional
        Initial concentration of the DNA species.
    **kwargs
        Additional keyword arguments passed to the parent `DNA` class.

    Attributes
    ----------
    dna : Species
        The DNA species representing the genetic construct.
    promoter : Promoter or None
        The promoter component controlling transcription.
    rbs : RBS or None
        The ribosome binding site controlling translation.
    transcript : Species or None
        The RNA transcript produced by transcription.
    protein : Species or None
        The protein product of translation.

    See Also
    --------
    DNA : Base class for DNA components.
    Promoter : Component representing transcriptional control elements.
    RBS : Component representing ribosome binding sites.
    RNA : Base class for RNA components.
    Protein : Base class for protein components.

    Notes
    -----
    The DNAassembly automatically coordinates its sub-components (promoter,
    RBS) by propagating updates to mechanisms, parameters, and mixtures. When
    mechanisms or parameters are added to the assembly, they are also added
    to the promoter and RBS (but never overwrite existing values in those
    components).

    The 'transcription' mechanism is used by the promoter to generate the
    species and reactions for transcript and the 'translation' mechanism is
    used by the RBS to generate the species and reactions for ribosome binding
    and protein production.

    For expression mixtures where transcription is bypassed, set
    `transcript=False` to enable direct translation from DNA to protein.  In
    this case, the 'transcription' mechanism will be used to generate the
    protein.

    Examples
    --------
    Create a simple constitutive gene expression construct:

    >>> # Basic assembly with automatic species creation
    >>> gene = bcp.DNAassembly(
    ...     name='gene_gfp',
    ...     promoter='pconst',
    ...     rbs='rbs1'
    ... )
    >>> gene.dna
    dna_gene_gfp
    >>> gene.transcript
    rna_gene_gfp
    >>> gene.protein
    protein_gene_gfp

    Create an assembly with custom species names:

    >>> gene = bcp.DNAassembly(
    ...     name='gene_reporter',
    ...     promoter='p_lac',
    ...     rbs='rbs_strong',
    ...     transcript='mRNA_gfp',
    ...     protein='protein_gfp'
    ... )

    Create an expression construct (no transcript):

    >>> gene = bcp.DNAassembly(
    ...     name='gene_direct',
    ...     promoter='p_const',
    ...     transcript=False,
    ...     protein='protein_x'
    ... )

    """

    def __init__(
        self,
        name: str,
        dna=None,
        promoter=None,
        transcript=None,
        rbs=None,
        protein=None,
        length=None,
        attributes=None,
        mechanisms=None,
        compartment=None,
        parameters=None,
        initial_concentration=None,
        **kwargs,
    ):
        """Initialize a DNAassembly.

        See class docstring for parameter descriptions.

        """
        self.promoter = None
        self.rbs = None
        self.transcript = None

        # This has to be called at the end so mechanisms are set for
        # the promoter, RBS, etc.
        DNA.__init__(
            self,
            name,
            length=length,
            mechanisms=mechanisms,
            parameters=parameters,
            initial_concentration=initial_concentration,
            attributes=attributes,
            compartment=compartment,
            **kwargs,
        )

        self.update_dna(dna, attributes=attributes)
        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_promoter(
            promoter, transcript=self.transcript, protein=self.protein
        )
        self.update_rbs(rbs, transcript=self.transcript, protein=self.protein)

    def get_species(self):
        """Get the primary DNA species of this assembly.

        Returns
        -------
        Species
            The DNA species representing this genetic construct.

        """
        return self.dna

    def set_mixture(self, mixture: Mixture) -> None:
        """Set the mixture containing this component and its sub-components.

        Also propagates the mixture reference to the promoter and RBS
        components if they exist.

        Parameters
        ----------
        mixture : Mixture
            The mixture object that contains this assembly.

        """
        self.mixture = mixture
        if self.promoter is not None:
            self.promoter.set_mixture(mixture)
        if self.rbs is not None:
            self.rbs.set_mixture(mixture)

    def update_dna(self, dna: Union[None, DNA, str], attributes=None) -> None:
        """Set or update the DNA species for this assembly.

        Creates a DNA species from the provided input and updates the DNA
        references in the promoter and RBS components if they exist.

        Parameters
        ----------
        dna : DNA, str, or None
            The DNA component, species name, or None. If None, creates a DNA
            species using the assembly's name. If a string, creates a new DNA
            species with that name. If a DNA object, uses it directly.
        attributes : list of str, optional
            Attribute tags to add to the DNA species.

        Notes
        -----
        This method automatically updates the `dna` attribute of the promoter
        and RBS components to maintain consistency across the assembly.

        """
        if dna is None:
            self.dna = self.set_species(
                self.name,
                material_type='dna',
                attributes=attributes,
                compartment=self.compartment,
            )
        else:
            self.dna = self.set_species(
                dna,
                material_type='dna',
                attributes=attributes,
                compartment=self.compartment,
            )

        if self.promoter is not None:
            self.promoter.dna = self.dna
        if self.rbs is not None:
            self.rbs.dna = self.dna

    def update_transcript(
        self, transcript: Union[None, RNA, str, bool], attributes=None
    ) -> None:
        """Set or update the RNA transcript for this assembly.

        Creates an RNA species from the provided input and updates the
        transcript references in the promoter and RBS components if they
        exist.

        Parameters
        ----------
        transcript : RNA, str, bool, or None
            The RNA component, species name, False, or None. If None, creates
            an RNA species using the assembly's name. If a string, creates a
            new RNA species with that name. If an RNA object, uses it
            directly. If False, sets transcript to None (used for expression
            mixtures without transcription).
        attributes : list of str, optional
            Attribute tags to add to the RNA species.

        Notes
        -----
        Setting `transcript=False` is used in expression mixtures where
        translation occurs directly from DNA without an explicit RNA
        intermediate.

        This method automatically updates the `transcript` attribute of the
        promoter and RBS components to maintain consistency.

        """
        if transcript is None:
            self.transcript = self.set_species(
                self.name,
                material_type='rna',
                attributes=attributes,
                compartment=self.compartment,
            )

        # this is used for expression mixtures where there is no transcript!
        elif transcript is False:
            self.transcript = None
        else:
            self.transcript = self.set_species(
                transcript,
                material_type='rna',
                attributes=attributes,
                compartment=self.compartment,
            )

        if self.promoter is not None:
            self.promoter.transcript = self.transcript
        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(
        self, protein: Union[None, Protein, str], attributes=None
    ) -> None:
        """Set or update the protein product for this assembly.

        Creates a Protein species from the provided input and updates the
        protein references in the promoter and RBS components if they exist.

        Parameters
        ----------
        protein : Protein, str, or None
            The Protein component, species name, or None. If None, creates a
            Protein species using the assembly's name. If a string, creates a
            new Protein species with that name. If a Protein object, uses it
            directly.
        attributes : list of str, optional
            Attribute tags to add to the Protein species.

        Notes
        -----
        This method automatically updates the `protein` attribute of the
        promoter and RBS components to maintain consistency across the
        assembly.

        """
        if protein is None:
            self.protein = self.set_species(
                self.name,
                material_type='protein',
                attributes=attributes,
                compartment=self.compartment,
            )
        else:
            self.protein = self.set_species(
                protein,
                material_type='protein',
                attributes=attributes,
                compartment=self.compartment,
            )

        if self.rbs is not None:
            self.rbs.protein = self.protein
        if self.promoter is not None:
            self.promoter.protein = self.protein

    def update_promoter(
        self,
        promoter: Union[Protein, str],
        transcript: RNA = None,
        protein: Protein = None,
    ) -> None:
        """Set or update the promoter component for this assembly.

        Creates a Promoter component from the provided input and propagates
        the assembly's parameters, mixture, and mechanisms to the promoter.

        Parameters
        ----------
        promoter : Promoter, str, or None
            The Promoter component, promoter name, or None. If None, no
            promoter is created. If a string, creates a default Promoter with
            that name using `Promoter.from_promoter`. If a Promoter object,
            uses it directly.
        transcript : RNA, optional
            The RNA transcript to associate with the promoter. If provided,
            updates the assembly's transcript before creating the promoter.
        protein : Protein, optional
            The protein product to associate with the promoter (used for some
            regulatory mechanisms).

        Notes
        -----
        This method automatically:

        - Propagates the assembly's parameter database to the promoter
        - Sets the promoter's mixture reference
        - Adds the assembly's mechanisms to the promoter (without overwriting
          existing promoter mechanisms)

        """
        if transcript is not None:
            self.update_transcript(transcript)

        if promoter is not None:
            self.promoter = Promoter.from_promoter(
                name=promoter,
                assembly=self,
                transcript=self.transcript,
                protein=protein,
            )
        else:
            self.promoter = None

        if self.promoter is not None:
            self.promoter.update_parameters(
                parameter_database=self.parameter_database,
                overwrite_parameters=False,
            )
            self.promoter.set_mixture(self.mixture)
            self.promoter.add_mechanisms(
                self.mechanisms, optional_mechanism=True
            )

    def update_rbs(
        self,
        rbs: Union[RBS, str],
        transcript: RNA = None,
        protein: Protein = None,
    ) -> None:
        """Set or update the ribosome binding site component.

        Creates an RBS component from the provided input and propagates the
        assembly's parameters, mixture, and mechanisms to the RBS.

        Parameters
        ----------
        rbs : RBS, str, or None
            The RBS component, RBS name, or None. If None, no RBS is created.
            If a string, creates a default RBS with that name using
            `RBS.from_rbs`. If an RBS object, uses it directly.
        transcript : RNA, optional
            The RNA transcript containing the RBS. If provided, updates the
            assembly's transcript before creating the RBS.
        protein : Protein, optional
            The protein product of translation. If provided, updates the
            assembly's protein before creating the RBS.

        Notes
        -----
        This method automatically:

        - Propagates the assembly's parameter database to the RBS
        - Sets the RBS's mixture reference
        - Adds the assembly's mechanisms to the RBS (without overwriting
          existing RBS mechanisms)

        """
        if protein is not None:
            self.update_protein(protein)

        if transcript is not None:
            self.update_transcript(transcript)

        if rbs is not None:
            self.rbs = RBS.from_rbs(
                name=rbs,
                assembly=self,
                protein=self.protein,
                transcript=self.transcript,
            )
        else:
            self.rbs = None

        if self.rbs is not None:
            self.rbs.update_parameters(
                parameter_database=self.parameter_database,
                overwrite_parameters=False,
            )
            self.rbs.set_mixture(self.mixture)
            self.rbs.add_mechanisms(self.mechanisms, optional_mechanism=True)

    def update_species(self) -> List[Species]:
        """Generate all species associated with this assembly.

        Collects species from the DNA, promoter, and RBS components during
        CRN compilation.

        Returns
        -------
        list of Species
            List containing the DNA species and all species generated by the
            promoter and RBS components.

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn` to collect all chemical species generated by
        this assembly.

        """
        # :return: list of Species that a DNAassemlby instance holds
        species = []
        species.append(self.dna)
        if self.promoter is not None:
            species += self.promoter.update_species()

        if self.rbs is not None:
            species += self.rbs.update_species()

        # deg_mech = self.get_mechanism(
        #    "rna_degradation", optional_mechanism=True)
        # if deg_mech is not None and self.promoter is not None \
        #    and self.transcript is not None:
        #    species += deg_mech.update_species(
        #        rna=self.transcript, component=self.promoter,
        #        part_id=self.transcript.name)

        return species

    def update_reactions(self) -> List[Reaction]:
        """Generate all reactions associated with this assembly.

        Collects reactions from the promoter and RBS components during CRN
        compilation.

        Returns
        -------
        list of Reaction
            List of all reactions generated by the promoter and RBS
            components, including transcription, translation, and regulatory
            reactions.

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn` to collect all chemical reactions generated by
        this assembly.

        """
        # :return: list of Reactions that a DNAassemlby instance holds.
        reactions = []
        if self.promoter is not None:
            reactions += self.promoter.update_reactions()

        if self.rbs is not None:
            reactions += self.rbs.update_reactions()

        # deg_mech = self.get_mechanism(
        #     "rna_degradation", optional_mechanism=True)
        # if deg_mech is not None and self.promoter is not None \
        #    and self.transcript is not None:
        #    reactions += deg_mech.update_reactions(
        #        rna=self.transcript, component=self.promoter,
        #        part_id=self.transcript.name)

        return reactions

    def update_parameters(
        self,
        parameter_file: str = None,
        parameters: ParameterDatabase = None,
        overwrite_parameters: bool = True,
    ) -> None:
        """Update parameters for the assembly and its sub-components.

        Propagates parameter updates to the DNA assembly itself and to the
        promoter and RBS components if they exist.

        Parameters
        ----------
        parameter_file : str, optional
            Path to a CSV or TSV parameter file to load.
        parameters : ParameterDatabase, optional
            ParameterDatabase object to merge with the assembly's parameters.
        overwrite_parameters : bool, default=True
            If True, new parameter values overwrite existing ones. If False,
            existing parameters are preserved.

        Notes
        -----
        This method calls `update_parameters` on:

        1. The parent DNA class (updating the DNA's parameters)
        2. The promoter component (if it exists)
        3. The RBS component (if it exists)

        """
        DNA.update_parameters(
            self=self,
            parameter_file=parameter_file,
            parameters=parameters,
            overwrite_parameters=overwrite_parameters,
        )

        if self.promoter is not None:
            self.promoter.update_parameters(
                parameter_file=parameter_file,
                parameters=parameters,
                overwrite_parameters=overwrite_parameters,
            )

        if self.rbs is not None:
            self.rbs.update_parameters(
                parameter_file=parameter_file,
                parameters=parameters,
                overwrite_parameters=overwrite_parameters,
            )

    def add_mechanism(
        self,
        mechanism: Mechanism,
        mech_type: str = None,
        overwrite: bool = False,
        optional_mechanism: bool = False,
    ) -> None:
        """Add a mechanism to the assembly and its sub-components.

        Adds the mechanism to the assembly's mechanism dictionary and
        propagates it to the promoter and RBS components without overwriting
        their existing mechanisms.

        Parameters
        ----------
        mechanism : Mechanism
            The mechanism object to add.
        mech_type : str, optional
            The mechanism type key. If None, uses the mechanism's
            `mechanism_type` attribute.
        overwrite : bool, default=False
            If True, overwrites existing mechanisms with the same type in the
            assembly. If False, raises ValueError for duplicate types.
        optional_mechanism : bool, default=False
            If True, suppresses ValueError when a mechanism key conflict
            occurs in the assembly and `overwrite` is False.

        Notes
        -----
        The mechanism is always added to the promoter and RBS with
        `optional_mechanism=True`, meaning it will never overwrite existing
        mechanisms in those components even if `overwrite=True`.

        """
        Component.add_mechanism(
            self,
            mechanism,
            mech_type=mech_type,
            overwrite=overwrite,
            optional_mechanism=optional_mechanism,
        )

        if self.promoter is not None:
            self.promoter.add_mechanism(
                mechanism, mech_type=mech_type, optional_mechanism=True
            )

        if self.rbs is not None:
            self.rbs.add_mechanism(
                mechanism, mech_type=mech_type, optional_mechanism=True
            )

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self)
        if self.promoter is not None:
            txt += '\n\t' + repr(self.promoter)
            txt += '\n\ttranscript = ' + repr(self.transcript)
        if self.rbs is not None:
            txt += '\n\t' + repr(self.rbs)
            txt += '\n\tprotein = ' + repr(self.protein)
        return txt
