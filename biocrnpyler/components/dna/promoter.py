#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import itertools as it
from warnings import warn

from ...core.species import Species
from ...mechanisms.binding import (
    CombinatorialCooperativeBinding,
    OneStepCooperativeBinding,
)
from ...mechanisms.txtl import (
    NegativeHillTranscription,
    PositiveHillTranscription,
)
from ..basic import DNA, RNA
from .construct import DNApart


class Promoter(DNApart):
    """Basic promoter component for constitutive transcription.

    A promoter represents a DNA regulatory element that controls transcription
    of an RNA transcript. This base class implements constitutive
    (unregulated) transcription. The component uses the 'transcription'
    mechanism to generate species and reactions during CRN compilation. The
    promoter must be included in a `DNAassembly` or `DNAconstruct` to function
    properly.

    Parameters
    ----------
    name : str
        Name of the promoter.
    assembly : DNAassembly, optional
        The DNA assembly containing this promoter. If provided, the assembly's
        name is used to create the default transcript.
    transcript : RNA, str, or None, optional
        The RNA transcript produced by this promoter. If None and `assembly`
        is provided, creates an RNA species using the assembly's name. Can be
        a list of transcripts for multi-cistronic operons.
    length : int, default=0
        Length of the promoter sequence in base pairs.
    mechanisms : dict or list, optional
        Custom mechanisms for this promoter, overriding mixture defaults.
    parameters : dict, optional
        Parameter values specific to this promoter.
    protein : Protein, str, list, or None, optional
        Protein product(s) for expression mixtures where transcription is
        bypassed. Can be a single protein, list of proteins, or None.
    dna_to_bind : DNA or Species, optional
        The DNA species that serves as the transcription template. If None,
        uses the assembly's DNA when available.
    **kwargs
        Additional keyword arguments passed to the parent `DNApart` class.

    Attributes
    ----------
    transcript : Species, list of Species, or None
        The RNA transcript(s) produced by transcription.
    protein : list of Species or None
        Protein product(s) for expression systems.
    length : int
        Length of the promoter in base pairs.
    dna_to_bind : Species or None
        The DNA species used as transcription template.

    See Also
    --------
    RegulatedPromoter : Promoter with independent regulator binding.
    ActivatablePromoter_Hill : Promoter with Hill function activation.
    RepressiblePromoter_Hill : Promoter with Hill function repression.
    CombinatorialPromoter : Promoter with combinatorial regulation.
    DNAassembly : Container for promoters and genetic constructs.

    Notes
    -----
    Promoters cannot have initial concentrations set directly. Initial
    conditions must be set on the containing `DNAassembly` or `DNAconstruct`.

    Examples
    --------
    Create a basic constitutive promoter:

    >>> promoter = bcp.Promoter(
    ...     name='pconst',
    ...     transcript='mRNA_gfp'
    ... )

    Create a promoter within an assembly:

    >>> assembly = bcp.DNAassembly(name='gene_x')
    >>> promoter = bcp.Promoter(
    ...     name='p_lac',
    ...     assembly=assembly
    ... )

    """

    def __init__(
        self,
        name,
        assembly=None,
        transcript=None,
        length=0,
        mechanisms=None,
        parameters=None,
        protein=None,
        dna_to_bind=None,
        **kwargs,
    ):
        self._dna_bind = dna_to_bind
        self.length = length

        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type='rna')
        elif isinstance(transcript, list):
            self.transcript = transcript
        else:
            self.transcript = self.set_species(
                transcript, material_type='rna'
            )

        if isinstance(protein, str):
            self.protein = [
                self.set_species(protein, material_type='protein')
            ]
        elif isinstance(protein, Species):
            self.protein = [protein]
        elif isinstance(protein, list):
            self.protein = protein
        else:
            self.protein = None
        # Promoter should not have initial conditions. These need to
        # be in DNAAssembly or DNAConstruct
        if (
            'initial_concentration' in kwargs.values()
            and kwargs['initial_concentration'] is not None
        ):
            raise AttributeError(
                "Cannot set initial_concentration of a Promoter. Must set "
                "initial_concentration for the DNAassembly or DNAConstruct"
            )
        if (
            'initial_condition_dictionary' in kwargs.values()
            and kwargs['initial_condition_dictionary'] is not None
        ):
            raise AttributeError(
                "Cannot set initial_condition_dictionary of a Promoter. Must "
                "set initial_condition_dictionary for the DNAassembly or "
                "DNAconstruct."
            )

        DNApart.__init__(
            self,
            name=name,
            mechanisms=mechanisms,
            parameters=parameters,
            assembly=assembly,
            **kwargs,
        )

    def update_species(self):
        """Generate species associated with this promoter.

        Calls the transcription mechanism to generate species for constitutive
        transcription from the DNA template.

        Returns
        -------
        list of Species
            List of species generated by the transcription mechanism,
            including RNA polymerase-DNA complexes and transcripts.

        """
        mech_tx = self.get_mechanism('transcription')
        species = []

        species += mech_tx.update_species(
            dna=self.dna_to_bind,
            transcript=self.transcript,
            protein=self.get_protein_for_expression(),
            component=self,
            part_id=self.name,
        )

        return species

    @property
    def dna_to_bind(self):
        """Species or None: DNA species used as transcription template.

        If not explicitly set, defaults to the assembly's DNA species.

        """
        if self._dna_bind is None:
            if self.assembly is None:
                return None
            else:
                self._dna_bind = self.assembly.dna
        return self._dna_bind

    @dna_to_bind.setter
    def dna_to_bind(self, value):
        """Set the DNA species used as transcription template.

        Parameters
        ----------
        value : Species or None
            The DNA species to use for transcription.

        """
        self._dna_bind = value

    def get_species(self):
        """Get the primary species associated with this promoter.

        Returns
        -------
        None
            Promoters do not have a primary species; they are part of a DNA
            assembly.

        """
        return None

    def update_reactions(self):
        """Generate reactions associated with this promoter.

        Calls the transcription mechanism to generate reactions for
        constitutive transcription from the DNA template.

        Returns
        -------
        list of Reaction
            List of reactions generated by the transcription mechanism,
            including RNA polymerase binding and transcript production.

        """
        mech_tx = self.get_mechanism('transcription')
        reactions = []

        reactions += mech_tx.update_reactions(
            dna=self.dna_to_bind,
            component=self,
            part_id=self.name,
            complex=None,
            transcript=self.transcript,
            protein=self.get_protein_for_expression(),
        )
        return reactions

    def update_component(self, internal_species=None, **kwargs):
        """Create a copy of the promoter with updated DNA binding target.

        Used for component enumeration when promoters are part of larger
        constructs that need to be duplicated with different species.

        Parameters
        ----------
        internal_species : Species, optional
            The new DNA species to use as the binding target.
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        Promoter or None
            A shallow copy of this promoter with the updated `dna_to_bind`
            attribute. Returns None if parent is RNA.

        Raises
        ------
        TypeError
            If parent is neither DNA nor RNA construct.

        """
        if isinstance(self.parent, RNA):
            return None
        elif isinstance(self.parent, DNA):
            out_component = copy.copy(self)
            out_component.dna_to_bind = internal_species
            return out_component
        else:
            raise TypeError(
                f"Unknown parent class {type(self.parent)}, expect either "
                "DNAconstruct or RNAconstruct"
            )

    # Used for expression mixtures where transcripts are replaced by proteins
    def get_protein_for_expression(self):
        """Get protein product for expression mixtures.

        In expression mixtures, transcription may be bypassed and translation
        may occur directly from DNA. This method returns the protein product
        when the gene is expressed.

        Returns
        -------
        list of Species or None
            The protein product(s) if transcript is None, otherwise None.

        Notes
        -----
        This is used by expression mixtures where the transcript species is
        omitted and translation occurs directly.

        """
        return self.protein

    @classmethod
    def from_promoter(cls, name, assembly, transcript, protein):
        """Create a promoter instance from another promoter or string.

        Factory method for creating promoters from various input types.

        Parameters
        ----------
        name : Promoter or str
            Either a string name for a new promoter, or an existing Promoter
            object to copy.
        assembly : DNAassembly
            The assembly containing this promoter.
        transcript : RNA or str
            The RNA transcript produced by this promoter.
        protein : Protein, str, or list
            The protein product(s) for expression mixtures.

        Returns
        -------
        Promoter
            A new Promoter instance. If `name` is a Promoter, returns a
            deep copy with updated assembly, transcript, and protein
            attributes.

        Raises
        ------
        TypeError
            If `name` is neither a string nor a Promoter.

        """
        if isinstance(name, Promoter):
            promoter_instance = copy.deepcopy(name)
            promoter_instance.assembly = assembly
            promoter_instance.transcript = transcript
            promoter_instance.protein = protein
        elif isinstance(name, str):
            promoter_instance = cls(
                name=name,
                assembly=assembly,
                transcript=transcript,
                protein=protein,
            )
        else:
            raise TypeError(
                "Promoter can be initialized from string or another "
                f"promoter! We got {type(name)}"
            )
        return promoter_instance


class RegulatedPromoter(Promoter):
    """Promoter with simple independent regulatory binding.

    A regulated promoter allows multiple regulatory proteins (activators or
    repressors) to bind independently to the promoter DNA. Each regulator
    binds independently, and transcription can occur from both bound and
    unbound states with different rates. The component uses the 'binding'
    mechanism (`OneStepCooperativeBinding` by default) to generate
    DNA-regulator complexes and the 'transcription' mechanism to generate
    transcription reactions from each regulatory state.

    Parameters
    ----------
    name : str
        Name of the promoter.
    regulators : Species, str, or list
        Regulator species that bind to the promoter. Can be a single
        regulator or a list. Each regulator binds independently.
    leak : bool, default=True
        If True, allows transcription from the unbound promoter state (leak
        expression). If False, only bound states transcribe.
    assembly : DNAassembly, optional
        The assembly containing this promoter.
    transcript : RNA or str, optional
        The RNA transcript produced by this promoter.
    length : int, default=0
        Length of the promoter in base pairs.
    mechanisms : dict or list, optional
        Custom mechanisms for this promoter.
    parameters : dict, optional
        Parameter values specific to this promoter.
    **kwargs
        Additional keyword arguments passed to parent class.

    Attributes
    ----------
    regulators : list of Species
        List of protein species that regulate this promoter.
    leak : bool
        Whether leak transcription (from unbound state) is allowed.
    complexes : list of Species
        List of DNA-regulator complexes generated during compilation.

    See Also
    --------
    Promoter : Base promoter class.
    ActivatablePromoter_Hill : Hill function-based activation.
    RepressiblePromoter_Hill : Hill function-based repression.
    CombinatorialPromoter : Combinatorial regulation logic.

    Notes
    -----
    Each regulator binds independently, creating multiple DNA-protein
    complexes. Transcription can occur from each complex with different
    parameters identified by part_id.

    The leak behavior allows modeling of constitutive expression that occurs
    even without regulator binding.

    Examples
    --------
    Create a promoter with a single regulator:

    >>> promoter = bcp.RegulatedPromoter(
    ...     name='p_reg',
    ...     regulators='protein_TF',
    ...     leak=True
    ... )

    Create a promoter with multiple independent regulators:

    >>> promoter = bcp.RegulatedPromoter(
    ...     name='p_multi',
    ...     regulators=['protein_TF1', 'protein_TF2'],
    ...     leak=False
    ... )

    """

    def __init__(
        self,
        name: str,
        regulators,
        leak=True,
        assembly=None,
        transcript=None,
        length=0,
        mechanisms=None,
        parameters=None,
        **kwargs,
    ):
        Promoter.__init__(
            self,
            name=name,
            assembly=assembly,
            transcript=transcript,
            length=length,
            mechanisms=mechanisms,
            parameters=parameters,
            **kwargs,
        )

        if not isinstance(regulators, list):
            regulators = [regulators]

        self.regulators = []
        for regulator in regulators:
            self.regulators += [
                self.set_species(regulator, material_type='protein')
            ]

        self.leak = leak

        self.add_mechanism(OneStepCooperativeBinding(), 'binding')
        self.complexes = []

    def update_species(self):
        """Generate species for regulated transcription.

        Creates species for regulator-DNA binding and transcription from each
        regulatory state. Generates DNA-regulator complexes and identifies
        which complexes can transcribe.

        Returns
        -------
        list of Species
            List containing all DNA-regulator complexes and species generated
            by the transcription mechanism for each regulatory state.

        Notes
        -----
        The method generates:

        - DNA-regulator binding complexes for each regulator
        - Transcription-related species for unbound DNA (if leak is True)
        - Transcription-related species for each DNA-regulator complex

        Complexes are stored in `self.complexes` for use in
        `update_reactions`.

        """
        mech_tx = self.get_mechanism('transcription')
        mech_b = self.get_mechanism('binding')
        species = []

        self.complexes = []
        if self.leak is not False:
            species += mech_tx.update_species(
                dna=self.dna_to_bind,
                transcript=self.transcript,
                component=self,
                part_id=self.name + '_leak',
                protein=self.get_protein_for_expression(),
            )

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]

            species_b = mech_b.update_species(
                regulator,
                self.dna_to_bind,
                part_id=self.name + '_' + regulator.name,
                component=self,
                protein=self.get_protein_for_expression(),
            )
            species += species_b

            # Find complexes containing DNA and the regulator
            # DNA should be *not* a part of an OrderedPolymer for this to work
            # dna_simple = copy.deepcopy(self.dna_to_bind)
            # dna_simple.remove()
            for s in species_b:
                if s.contains_species_monomer(
                    self.dna_to_bind
                ) and s.contains_species_monomer(regulator):
                    self.complexes += [s]

                    species += mech_tx.update_species(
                        dna=s,
                        transcript=self.transcript,
                        protein=self.get_protein_for_expression(),
                        part_id=self.name + '_' + regulator.name,
                        component=self,
                    )
        return species

    def update_reactions(self):
        """Generate reactions for regulated transcription.

        Creates binding reactions for each regulator and transcription
        reactions for each regulatory state (bound and unbound).

        Returns
        -------
        list of Reaction
            List containing all binding reactions for regulators and
            transcription reactions for each DNA state.

        Notes
        -----
        Reactions are generated for:

        - Regulator binding and unbinding to DNA
        - Transcription from unbound DNA (if leak is True)
        - Transcription from each DNA-regulator complex

        """
        reactions = []
        mech_tx = self.get_mechanism('transcription')
        mech_b = self.get_mechanism('binding')

        if self.leak is not False:
            reactions += mech_tx.update_reactions(
                dna=self.dna_to_bind,
                component=self,
                part_id=self.name + '_leak',
                transcript=self.transcript,
                protein=self.get_protein_for_expression(),
            )

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]

            reactions += mech_b.update_reactions(
                regulator,
                self.dna_to_bind,
                component=self,
                part_id=self.name + '_' + regulator.name,
                protein=self.protein,
            )

            reactions += mech_tx.update_reactions(
                dna=complex_,
                component=self,
                part_id=self.name + '_' + regulator.name,
                transcript=self.transcript,
                protein=self.get_protein_for_expression(),
            )

        return reactions


class ActivatablePromoter_Hill(Promoter):
    r"""Promoter with Hill function-based activation.

    An activatable promoter models transcriptional activation by a single
    regulator species using Hill function kinetics. The component uses a
    'transcription' mechanism (`PositiveHillTranscription`) to generate
    species and reactions where the transcription rate increases with
    activator concentration following cooperative binding dynamics.

    Parameters
    ----------
    name : str
        Name of the promoter.
    activator : Species or str
        The activator protein species that enhances transcription.
    transcript : RNA or str, optional
        The RNA transcript produced by this promoter.
    leak : bool, default=False
        If True, allows basal transcription without activator. If False, no
        transcription occurs without activator binding.
    **kwargs
        Additional keyword arguments passed to parent `Promoter` class.

    Attributes
    ----------
    activator : Species
        The activator protein that enhances transcription.
    leak : bool
        Whether basal (leak) transcription is allowed.

    See Also
    --------
    RepressiblePromoter_Hill : Hill function-based repression.
    RegulatedPromoter : Independent regulator binding.
    PositiveHillTranscription : Mechanism for Hill activation.

    Notes
    -----
    The activation follows a Hill function:

    .. math::

        \text{rate} = k_{\text{max}} \frac{[A]^n}{K_d^n + [A]^n}
            + k_{\text{leak}}

    where [A] is activator concentration, n is the Hill coefficient, and
    :math:`K_d` is the dissociation constant.

    Examples
    --------
    Create an activatable promoter:

    >>> promoter = bcp.ActivatablePromoter_Hill(
    ...     name='p_ara',
    ...     activator='protein_AraC',
    ...     leak=True
    ... )

    """

    def __init__(
        self, name, activator, transcript=None, leak=False, **kwargs
    ):
        # Always call the superclass __init__() with **kwargs passed through
        Promoter.__init__(self, name=name, transcript=transcript, **kwargs)

        # Set the Regulator
        # Component.set_species(
        #     species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species,
        # string, or Component.
        self.activator = self.set_species(activator)

        self.leak = leak  # toggles whether or not there is a leak reaction

        # Non-default Mechanisms are added to the Component with
        # .add_mechanism
        self.add_mechanism(
            PositiveHillTranscription(), 'transcription', overwrite=True
        )

    def update_species(self, **kwargs):
        """Generate species for Hill function activation.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the transcription
            mechanism.

        Returns
        -------
        list of Species
            List of species generated by the transcription mechanism for
            Hill-based activation.

        """
        # Mechanisms are stored in an automatically created
        # dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.get_mechanism('transcription')

        species = []  # A list of species must be returned
        species += mech_tx.update_species(
            dna=self.dna_to_bind,
            transcript=self.transcript,
            regulator=self.activator,
            part_id=self.name + '_' + self.activator.name,
            leak=self.leak,
            component=self,
            protein=self.get_protein_for_expression(),
        )

        return species

    def update_reactions(self, **kwargs):
        """Generate reactions for Hill function activation.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the transcription
            mechanism.

        Returns
        -------
        list of Reaction
            List of reactions for Hill-based transcriptional activation.

        """
        mech_tx = self.get_mechanism('transcription')

        reactions = []  # a list of reactions must be returned

        reactions += mech_tx.update_reactions(
            dna=self.dna_to_bind,
            transcript=self.transcript,
            regulator=self.activator,
            component=self,
            part_id=self.name + '_' + self.activator.name,
            leak=self.leak,
            protein=self.get_protein_for_expression(),
        )

        return reactions


class RepressiblePromoter_Hill(Promoter):
    r"""Promoter with Hill function-based repression.

    A repressible promoter models transcriptional repression by a single
    regulator species using Hill function kinetics. The component uses a
    'transcription' mechanism (`NegativeHillTranscription`) to generate
    species and reactions where the transcription rate decreases with
    repressor concentration following cooperative binding dynamics.

    Parameters
    ----------
    name : str
        Name of the promoter.
    repressor : Species or str
        The repressor protein species that inhibits transcription.
    transcript : RNA or str, optional
        The RNA transcript produced by this promoter.
    leak : bool, default=False
        If True, allows residual transcription even at high repressor
        concentrations. If False, transcription is fully repressed.
    **kwargs
        Additional keyword arguments passed to parent `Promoter` class.

    Attributes
    ----------
    repressor : Species
        The repressor protein that inhibits transcription.
    leak : bool
        Whether leak transcription at high repressor is allowed.

    See Also
    --------
    ActivatablePromoter_Hill : Hill function-based activation.
    RegulatedPromoter : Independent regulator binding.
    NegativeHillTranscription : Mechanism for Hill repression.

    Notes
    -----
    The repression follows a Hill function:

    .. math::

        \text{rate} = k_{\text{max}} \frac{K_d^n}{K_d^n + [R]^n}
            + k_{\text{leak}}

    where [R] is repressor concentration, n is the Hill coefficient, and
    :math:`K_d` is the dissociation constant.

    Examples
    --------
    Create a repressible promoter:

    >>> promoter = bcp.RepressiblePromoter_Hill(
    ...     name='p_lac',
    ...     repressor='protein_LacI',
    ...     leak=False
    ... )

    """

    def __init__(
        self, name, repressor, transcript=None, leak=False, **kwargs
    ):
        # Always call the superclass __init__() with **kwargs passed through
        Promoter.__init__(self, name=name, transcript=transcript, **kwargs)

        # Set the Regulator
        # Component.set_species(
        #    species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species,
        # string, or Component.
        self.repressor = self.set_species(repressor)

        self.leak = leak  # toggles whether or not there is a leak reaction

        # Mechanisms are inherited from the Mixture unless set
        # specifically in self.default_mechanisms.
        self.add_mechanism(
            NegativeHillTranscription(), 'transcription', overwrite=True
        )

    def update_species(self, **kwargs):
        """Generate species for Hill function repression.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the transcription
            mechanism.

        Returns
        -------
        list of Species
            List of species generated by the transcription mechanism for
            Hill-based repression.

        """
        # Mechanisms are stored in an automatically created
        # dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.get_mechanism('transcription')

        species = []  # A list of species must be returned
        species += mech_tx.update_species(
            dna=self.dna_to_bind,
            transcript=self.transcript,
            regulator=self.repressor,
            component=self,
            part_id=self.name + '_' + self.repressor.name,
            leak=self.leak,
            protein=self.get_protein_for_expression(),
            **kwargs,
        )

        return species

    def update_reactions(self, **kwargs):
        """Generate reactions for Hill function repression.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the transcription
            mechanism.

        Returns
        -------
        list of Reaction
            List of reactions for Hill-based transcriptional repression.

        """
        mech_tx = self.get_mechanism('transcription')

        reactions = []  # a list of reactions must be returned

        reactions += mech_tx.update_reactions(
            dna=self.dna_to_bind,
            transcript=self.transcript,
            regulator=self.repressor,
            component=self,
            part_id=self.name + '_' + self.repressor.name,
            leak=self.leak,
            protein=self.get_protein_for_expression(),
            **kwargs,
        )
        return reactions


class CombinatorialPromoter(Promoter):
    """Promoter with combinatorial regulatory logic.

    A combinatorial promoter allows multiple regulators to bind cooperatively,
    where transcription behavior depends on the specific combination of bound
    regulators. The component uses the 'binding' mechanism
    (`CombinatorialCooperativeBinding`) to generate all possible DNA-regulator
    complexes and the 'transcription' mechanism to generate reactions for each
    regulatory state. This enables complex logic gates (AND, OR, NOR, etc.)
    and multi-input regulatory functions.

    Parameters
    ----------
    name : str
        Name of the promoter.
    regulators : Species, str, or list
        List of regulator species that can bind to the promoter. Regulators
        can bind in various combinations.
    leak : bool, default=False
        If True, allows transcription from promoter states not in
        `tx_capable_list` (including unbound state). If False, only states in
        `tx_capable_list` transcribe.
    assembly : DNAassembly, optional
        The assembly containing this promoter.
    transcript : RNA or str, optional
        The RNA transcript produced by this promoter.
    length : int, default=0
        Length of the promoter in base pairs.
    mechanisms : dict or list, optional
        Custom mechanisms for this promoter.
    parameters : dict, optional
        Parameter values specific to this promoter.
    protein : Protein, str, list, or None, optional
        Protein product(s) for expression mixtures.
    tx_capable_list : list of list, optional
        List specifying which combinations of bound regulators enable
        transcription. Each element is a list of regulator names (strings or
        Species). If None, all combinations enable transcription.
    cooperativity : dict, optional
        Dictionary mapping regulator names to their cooperativity values
        (Hill coefficients) for binding, e.g., {'regulator1': 2,
        'regulator2': 1}.
    **kwargs
        Additional keyword arguments passed to parent class.

    Attributes
    ----------
    regulators : list of Species
        Sorted list of protein regulators (sorted for consistency).
    cooperativity : dict or None
        Cooperativity values for each regulator.
    tx_capable_list : list of set
        List of regulator combinations (as sets) that enable transcription.
    leak : bool
        Whether leak transcription is allowed for non-capable states.
    complex_combinations : dict
        Dictionary mapping combinations to their DNA-regulator complexes.
    tx_capable_complexes : list of Species
        List of DNA complexes that can transcribe.
    leak_complexes : list of Species
        List of DNA complexes that transcribe with leak parameters.

    See Also
    --------
    RegulatedPromoter : Independent regulatory binding.
    CombinatorialCooperativeBinding : Binding mechanism used by this
        promoter.

    Notes
    -----
    Only combinations in `tx_capable_list` transcribe with full rate;
    others transcribe with leak rate (if leak is True) or not at all.

    Regulators are automatically sorted alphabetically to ensure consistent
    ordering when checking combinations.

    Examples
    --------
    Create an AND gate promoter (transcribes only when both bound):

    >>> promoter = bcp.CombinatorialPromoter(
    ...     name='p_and',
    ...     regulators=['TF1', 'TF2'],
    ...     tx_capable_list=[['TF1', 'TF2']],
    ...     leak=False
    ... )

    Create an OR gate promoter (transcribes when either is bound):

    >>> promoter = bcp.CombinatorialPromoter(
    ...     name='p_or',
    ...     regulators=['TF1', 'TF2'],
    ...     tx_capable_list=[['TF1'], ['TF2'], ['TF1', 'TF2']],
    ...     leak=False
    ... )

    Create a promoter with cooperative binding:

    >>> promoter = bcp.CombinatorialPromoter(
    ...     name='p_coop',
    ...     regulators=['TF1', 'TF2'],
    ...     cooperativity={'TF1': 2, 'TF2': 1},
    ...     tx_capable_list=[['TF1', 'TF2']]
    ... )

    """

    def __init__(
        self,
        name,
        regulators,
        leak=False,
        assembly=None,
        transcript=None,
        length=0,
        mechanisms=None,
        parameters=None,
        protein=None,
        tx_capable_list=None,
        cooperativity=None,
        **kwargs,
    ):
        Promoter.__init__(
            self,
            name=name,
            assembly=assembly,
            transcript=transcript,
            length=length,
            mechanisms=mechanisms,
            parameters=parameters,
            protein=protein,
            **kwargs,
        )

        if not isinstance(regulators, list):
            # you could give one string as a regulator
            regulators = [regulators]
        self.cooperativity = cooperativity
        self.regulators = []
        for regulator in regulators:
            self.regulators += [
                self.set_species(regulator, material_type='protein')
            ]

        # after we've sanitized the inputs, then sort
        self.regulators = sorted(self.regulators)
        # now let's work out the tx_capable_list
        if tx_capable_list is None:
            # if nothing is passed, that means everything transcribes
            allcomb = []
            for r in range(1, len(self.regulators) + 1):
                # make all combinations of regulators
                allcomb += [
                    set(a)
                    for a in it.combinations(
                        [a.name for a in self.regulators], r
                    )
                ]
            self.tx_capable_list = allcomb
        elif isinstance(tx_capable_list, list):
            # if the user passed a list then the user knows what they want
            newlist = []
            # this part converts any species in the tx_capable_list into a
            # string
            for element in tx_capable_list:
                sublist = []
                for specie in element:
                    if isinstance(specie, Species):
                        sublist += [specie.name]
                    else:
                        sublist += [specie]
                newlist += [sublist]
            self.tx_capable_list = [set(a) for a in newlist]

        self.leak = leak
        self.complex_combinations = {}
        self.tx_capable_complexes = []
        self.leak_complexes = []
        self.add_mechanism(
            CombinatorialCooperativeBinding(), 'binding', overwrite=True
        )

    def update_species(self):
        """Generate species for combinatorial regulatory logic.

        Creates all possible DNA-regulator complexes through cooperative
        binding and identifies which complexes enable transcription based on
        `tx_capable_list`.

        Returns
        -------
        list of Species
            List containing:

            - The unbound DNA
            - All regulator species
            - All DNA-regulator binding complexes
            - Transcription-related species for capable complexes
            - Transcription-related species for leak complexes (if leak is
              True)

        Notes
        -----
        The method classifies DNA-regulator complexes into two categories:

        - `tx_capable_complexes`: Complexes that match combinations in
          `tx_capable_list` and transcribe with full rate parameters
        - `leak_complexes`: Complexes not in `tx_capable_list` that
          transcribe with leak parameters (if leak is True)

        Complexes are stored in these attributes for use by
        `update_reactions`.

        """
        mech_tx = self.get_mechanism('transcription')
        mech_b = self.get_mechanism('binding')
        # set the tx_capable_complexes to nothing because we havent updated
        # species yet!
        self.tx_capable_complexes = []
        self.leak_complexes = []
        species = [self.dna_to_bind] + self.regulators
        self.complexes = []
        bound_species = mech_b.update_species(
            self.regulators,
            self.dna_to_bind,
            component=self,
            part_id=self.name,
            cooperativity=self.cooperativity,
            protein=self.protein,
        )
        # above is all the species with DNA bound to regulators. Now, we need
        # to extract only the ones which are transcribable

        if self.leak is not False:
            # this part takes care of the promoter not bound to anything
            species += mech_tx.update_species(
                dna=self.dna_to_bind,
                transcript=self.transcript,
                protein=self.get_protein_for_expression(),
                component=self,
                part_id=self.name + '_leak',
            )
            # self.leak_complexes += []

        for bound_complex in bound_species:
            species_inside = []
            for regulator in self.regulators:
                if bound_complex.contains_species_monomer(
                    regulator
                ) and bound_complex.contains_species_monomer(
                    self.dna_to_bind
                ):
                    species_inside += [regulator.name]
            if set(species_inside) in [set(a) for a in self.tx_capable_list]:
                # only the transcribable complexes in tx_capable_list get
                # transcription reactions
                tx_capable_species = mech_tx.update_species(
                    dna=bound_complex,
                    transcript=self.transcript,
                    protein=self.get_protein_for_expression(),
                    component=self,
                    part_id=self.name,
                )
                species += tx_capable_species
                self.tx_capable_complexes += [bound_complex]
            else:
                # in this case there's a combination of regulators which does
                # not feature in tx_capable_list.  This means:
                # 1) this complex does nothing so don't add it
                # 2) we said we wanted leak, so then you should add this,
                #    but with the "_leak" parameters (that happens in
                #    update_reactions)
                if self.leak is not False:
                    leak_species = mech_tx.update_species(
                        dna=bound_complex,
                        transcript=self.transcript,
                        protein=self.get_protein_for_expression(),
                        component=self,
                        part_id=self.name + '_leak',
                    )
                    species += leak_species
                    self.leak_complexes += [bound_complex]
        species += bound_species
        return species

    def update_reactions(self):
        """Generate reactions for combinatorial regulatory logic.

        Creates binding reactions for all regulator combinations and
        transcription reactions for capable and leak complexes.

        Returns
        -------
        list of Reaction
            List containing:

            - Cooperative binding reactions for all regulator combinations
            - Transcription reactions from unbound DNA (if leak is True)
            - Transcription reactions from capable complexes
            - Transcription reactions from leak complexes (if leak is True)

        Warns
        -----
        UserWarning
            If no complexes can transcribe after calling `update_species`.

        Notes
        -----
        This method automatically calls `update_species` if the complex
        lists are empty, ensuring that species generation occurs before
        reaction generation.

        Each transcription reaction uses a unique part_id that identifies the
        regulatory state, constructed from the promoter name and bound
        regulators (e.g., 'p_name_TF1_TF2_RNAP').

        """
        reactions = []
        mech_tx = self.get_mechanism('transcription')
        mech_b = self.get_mechanism('binding')

        if self.leak is not False:
            # once again the DNA not bound to anything gets special treatment
            reactions += mech_tx.update_reactions(
                dna=self.dna_to_bind,
                component=self,
                part_id=self.name + '_leak',
                transcript=self.transcript,
                protein=self.protein,
            )
        # the binding reactions happen no matter what
        reactions += mech_b.update_reactions(
            self.regulators,
            self.dna_to_bind,
            component=self,
            part_id=self.name,
            cooperativity=self.cooperativity,
            protein=self.protein,
        )
        if (
            self.tx_capable_complexes is None
        ) or self.tx_capable_complexes == []:
            # this could mean we haven't run update_species() yet
            self.update_species()

            if self.tx_capable_complexes == []:
                if self.leak_complexes is None or self.leak_complexes == []:
                    # if it's still zero after running update_species
                    # then we could be in trouble
                    warn(
                        "nothing can transcribe from combinatorial promoter "
                        f"{self.name}"
                    )

        if len(self.tx_capable_complexes) > 0:
            for specie in self.tx_capable_complexes:
                tx_partid = self.name
                for part in specie.species_set:
                    # construct the name of the promoter with regulators bound
                    if part.material_type == 'dna':
                        # the DNA doesn't matter
                        pass
                    else:
                        # put in the regulators!
                        tx_partid += '_' + part.name
                if tx_partid[0] == '_':
                    # this will only happen if the name of the dna is ""
                    tx_partid = tx_partid[1:]
                # if it's bound to RNAP then it transcribes, right?
                tx_partid = tx_partid + '_RNAP'
                reactions += mech_tx.update_reactions(
                    dna=specie,
                    component=self,
                    part_id=tx_partid,
                    transcript=self.transcript,
                    protein=self.get_protein_for_expression(),
                )
        if len(self.leak_complexes) > 0:
            for specie in self.leak_complexes:
                # in this case every reaction uses the "promoter_leak" partid
                leak_partid = self.name + '_leak'
                reactions += mech_tx.update_reactions(
                    dna=specie,
                    component=self,
                    part_id=leak_partid,
                    transcript=self.transcript,
                    protein=self.get_protein_for_expression(),
                )

        return reactions


# Legacy classes
RepressiblePromoter = RepressiblePromoter_Hill
ActivatablePromoter = ActivatablePromoter_Hill
