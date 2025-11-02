# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..components.basic import Protein
from ..components.dna.assembly import DNAassembly
from ..core.chemical_reaction_network import ChemicalReactionNetwork
from ..core.mechanism import EmptyMechanism
from ..core.mixture import Mixture
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import BasicCatalysis, MichaelisMenten
from ..mechanisms.global_mechanisms import Degradation_mRNA_MM, Dilution
from ..mechanisms.txtl import (
    OneStepGeneExpression,
    SimpleTranscription,
    SimpleTranslation,
    Transcription_MM,
    Translation_MM,
)


class ExpressionDilutionMixture(Mixture):
    """In vivo gene expression with dilution but without cellular machinery.

    A simplified mixture that models gene expression as a single direct
    reaction from DNA to protein, without explicitly representing
    transcription and translation as separate processes or cellular machinery
    (ribosomes, polymerases). This mixture lumps transcription and
    translation into a single 'expression' reaction and includes global
    dilution to model cell growth and division effects on all non-DNA
    species.

    This mixture is appropriate for coarse-grained models of in vivo gene
    expression where mRNA dynamics are negligible and growth dilution is
    important.

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and `Mechanism` objects as values, a
        list of `Mechanism` objects, or a single `Mechanism`.
    components : list of Component or Component, optional
        Components to include in the mixture. Components are deep-copied when
        added to prevent modification of original objects.
    parameters : dict, optional
        Dictionary of parameter values. Keys follow the format
        (mechanism, part_id, param_name).
    compartment : Compartment, optional
        Default compartment for all components and species in this mixture.
    parameter_file : str, optional
        Path to a CSV or TSV file containing parameters to load.
    overwrite_parameters : bool, default=False
        If True, parameters from file/dict overwrite existing parameters.
        If False, existing parameters are preserved.
    global_mechanisms : dict, list, or GlobalMechanism, optional
        Global mechanisms that apply to all species after component
        compilation (e.g., dilution, global degradation). Can be a dict,
        list, or single `GlobalMechanism`.
    species : list of Species or Species, optional
        Additional species to add directly to the CRN without going through
        component compilation.
    initial_condition_dictionary : dict, optional
        Dictionary mapping species to initial concentration values. Deprecated
        in favor of using parameters with mechanism='initial concentration'.
    global_component_enumerators : list, optional
        List of global component enumerators for advanced component generation
        patterns (e.g., creating all pairwise interactions).
    global_recursion_depth : int, default=4
        Maximum recursion depth for global component enumeration during
        compilation.
    local_recursion_depth : int, optional
        Maximum recursion depth for local component enumeration. If None,
        defaults to `global_recursion_depth + 2`.

    Attributes
    ----------
    name : str
        Name of the mixture.
    compartment : Compartment or None
        Default compartment for the mixture.
    components : list of Component
        List of components in the mixture (deep copies of added components).
    mechanisms : dict
        Dictionary of default mechanisms, keyed by mechanism type (str).
    global_mechanisms : dict
        Dictionary of global mechanisms, keyed by mechanism type (str).
    parameter_database : ParameterDatabase
        Database storing all parameters for this mixture.
    added_species : list of Species
        List of species added directly to the mixture.
    global_component_enumerators : list
        List of global component enumerators.
    global_recursion_depth : int
        Recursion depth for global component enumeration.
    local_recursion_depth : int
        Recursion depth for local component enumeration.
    crn : ChemicalReactionNetwork or None
        The compiled CRN, created by calling `compile_crn`.

    See Also
    --------
    ExpressionExtract : Expression without dilution.
    SimpleTxTlDilutionMixture : TX-TL with dilution.
    TxTlDilutionMixture : TX-TL with machinery and dilution.
    Mixture : Base class for all mixtures.

    Notes
    -----
    Default mechanisms included:

    - 'transcription' : `OneStepGeneExpression` - Single-step gene
      expression (DNA --> DNA + Protein) without intermediate mRNA
    - 'translation' : `EmptyMechanism` - Dummy mechanism that generates no
      reactions (translation is disabled)
    - 'catalysis' : `BasicCatalysis` - Simple catalytic reactions without
      explicit enzyme binding
    - 'binding' : `One_Step_Binding` - Simple multi-species binding
    - 'dilution' : `Dilution` - Global dilution mechanism (Species --> ∅)
      applied to all non-DNA species to model growth/division

    Key features of this mixture:

    - No explicit transcription or translation steps
    - No cellular machinery (RNAP, ribosomes, RNases)
    - No intermediate mRNA species
    - Global dilution of all species except DNA
    - Models growth dilution effects in vivo
    - Simplified parameter space

    When compiled, this mixture automatically disables transcript generation
    in DNAassemblies that produce proteins, routing expression directly from
    DNA to protein.

    Common applications include:

    - In vivo gene circuit modeling with growth effects
    - Steady-state gene expression in growing cells
    - Models where mRNA dynamics are negligible
    - High-level circuit design with dilution

    Examples
    --------
    Create an in vivo expression mixture with dilution for GFP:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.ExpressionDilutionMixture(
    ...     name='cell_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/cell_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(self, name='', **kwargs):
        Mixture.__init__(self, name=name, **kwargs)

        # Create default mechanisms for Gene Expression
        dummy_translation = EmptyMechanism(
            name='dummy_translation', mechanism_type='translation'
        )
        mech_expression = OneStepGeneExpression()
        mech_cat = BasicCatalysis()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_expression.mechanism_type: mech_expression,
            dummy_translation.mechanism_type: dummy_translation,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)

        # Create global mechanism for dilution
        dilution_mechanism = Dilution(
            name='dilution', filter_dict={'dna': False}, default_on=True
        )
        global_mechanisms = {'dilution': dilution_mechanism}
        self.add_mechanisms(global_mechanisms, overwrite=None)

    def compile_crn(self, **kwargs) -> ChemicalReactionNetwork:
        """Compile CRN with transcript generation disabled in gene expression.

        Overrides the parent `compile_crn` method to automatically disable
        transcript generation in DNAassemblies that produce proteins. This
        ensures that gene expression proceeds directly from DNA to protein
        without intermediate mRNA species.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the parent Mixture
            `compile_crn` method.

        Returns
        -------
        ChemicalReactionNetwork
            Compiled chemical reaction network with expression and dilution
            reactions.

        Notes
        -----
        This method automatically modifies DNAassemblies before compilation:

        - For assemblies with a protein product, sets transcript to False
        - RNA-only assemblies (no protein) are not affected
        - Mechanisms receive protein instead of transcript when transcript
          is disabled

        This behavior enables the single-step expression mechanism to route
        production directly to protein.

        See `Mixture.compile_crn
        <biocrnpyler.core.mixture.Mixture.compile_crn>` for a more detailed
        description of the parent method behavior.

        """
        for component in self.components:
            if isinstance(component, DNAassembly):
                # Only turn off transcription for an Assembly that makes a
                # Protein.  Some assemblies might only make RNA!
                if component.protein is not None:
                    # This will turn off transcription and set
                    # Promoter.transcript = False Mechanisms that recieve no
                    # transcript but a protein will use the protein instead.
                    component.update_transcript(False)

        # Call the superclass function
        return Mixture.compile_crn(self, **kwargs)


class SimpleTxTlDilutionMixture(Mixture):
    """In vivo TX-TL with simple mechanisms and continuous dilution.

    A mixture that models transcription and translation as separate catalytic
    reactions without explicitly representing cellular machinery (RNAP,
    ribosomes). This mixture uses simple mass-action kinetics where DNA and
    mRNA act as catalysts for transcript and protein production,
    respectively. Includes global dilution to model cell growth and division
    effects, plus separate RNA degradation to model endonuclease activity.

    This mixture is appropriate for in vivo gene expression models where
    machinery is not limiting but explicit TX-TL steps and growth dilution
    are important.

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and `Mechanism` objects as values, a
        list of `Mechanism` objects, or a single `Mechanism`.
    components : list of Component or Component, optional
        Components to include in the mixture. Components are deep-copied when
        added to prevent modification of original objects.
    parameters : dict, optional
        Dictionary of parameter values. Keys follow the format
        (mechanism, part_id, param_name).
    compartment : Compartment, optional
        Default compartment for all components and species in this mixture.
    parameter_file : str, optional
        Path to a CSV or TSV file containing parameters to load.
    overwrite_parameters : bool, default=False
        If True, parameters from file/dict overwrite existing parameters.
        If False, existing parameters are preserved.
    global_mechanisms : dict, list, or GlobalMechanism, optional
        Global mechanisms that apply to all species after component
        compilation (e.g., dilution, global degradation). Can be a dict,
        list, or single `GlobalMechanism`.
    species : list of Species or Species, optional
        Additional species to add directly to the CRN without going through
        component compilation.
    initial_condition_dictionary : dict, optional
        Dictionary mapping species to initial concentration values. Deprecated
        in favor of using parameters with mechanism='initial concentration'.
    global_component_enumerators : list, optional
        List of global component enumerators for advanced component generation
        patterns (e.g., creating all pairwise interactions).
    global_recursion_depth : int, default=4
        Maximum recursion depth for global component enumeration during
        compilation.
    local_recursion_depth : int, optional
        Maximum recursion depth for local component enumeration. If None,
        defaults to `global_recursion_depth + 2`.

    Attributes
    ----------
    name : str
        Name of the mixture.
    compartment : Compartment or None
        Default compartment for the mixture.
    components : list of Component
        List of components in the mixture (deep copies of added components).
    mechanisms : dict
        Dictionary of default mechanisms, keyed by mechanism type (str).
    global_mechanisms : dict
        Dictionary of global mechanisms, keyed by mechanism type (str).
    parameter_database : ParameterDatabase
        Database storing all parameters for this mixture.
    added_species : list of Species
        List of species added directly to the mixture.
    global_component_enumerators : list
        List of global component enumerators.
    global_recursion_depth : int
        Recursion depth for global component enumeration.
    local_recursion_depth : int
        Recursion depth for local component enumeration.
    crn : ChemicalReactionNetwork or None
        The compiled CRN, created by calling `compile_crn`.

    See Also
    --------
    ExpressionDilutionMixture : Single-step expression with dilution.
    TxTlDilutionMixture : TX-TL with machinery and dilution.
    SimpleTxTlExtract : TX-TL without dilution.
    Mixture : Base class for all mixtures.

    Notes
    -----
    Default mechanisms included:

    - 'transcription' : `SimpleTranscription` - Simple catalytic
      transcription (DNA --> DNA + mRNA) without explicit RNAP binding
    - 'translation' : `SimpleTranslation` - Simple catalytic translation
      (mRNA --> mRNA + Protein) without explicit ribosome binding
    - 'catalysis' : `BasicCatalysis` - Simple catalytic reactions without
      explicit enzyme binding
    - 'binding' : `One_Step_Binding` - Simple multi-species binding
    - 'dilution' : `Dilution` - Global dilution mechanism (Species --> ∅)
      applied to all non-DNA species to model growth/division
    - 'rna_degradation' : `Dilution` - Separate RNA degradation mechanism
      (mRNA --> ∅) applied to all RNA species to model endonuclease
      activity

    Key features of this mixture:

    - Explicit transcription and translation steps
    - Intermediate mRNA species
    - Simple mass-action kinetics (no enzyme binding)
    - No cellular machinery (RNAP, ribosomes)
    - Global dilution of all non-DNA species
    - Separate RNA degradation (faster than dilution)
    - Models growth effects in vivo

    Common applications include:

    - In vivo gene circuit modeling with growth
    - Models where machinery is not limiting
    - Constitutive or weakly regulated promoters in growing cells
    - mRNA dynamics with degradation and dilution

    Examples
    --------
    Create a simple in vivo TX-TL mixture with dilution for GFP:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.SimpleTxTlDilutionMixture(
    ...     name='cell_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/cell_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(self, name='', **kwargs):
        # Always call the superclass __init__ with **kwargs
        Mixture.__init__(self, name=name, **kwargs)

        # Create TxTl Mechanisms
        # Transcription will not involve machinery
        simple_transcription = SimpleTranscription()
        simple_translation = SimpleTranslation()
        mech_cat = BasicCatalysis()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            simple_transcription.mechanism_type: simple_transcription,
            simple_translation.mechanism_type: simple_translation,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)

        # Global Dilution Mechanisms
        # By Default Species are diluted S-->0 Unless:
        # They are of type 'dna'
        # They have the attribute 'machinery'
        dilution_mechanism = Dilution(
            filter_dict={'dna': False}, default_on=True
        )
        deg_mrna = Dilution(
            name='rna_degradation',
            filter_dict={'rna': True},
            default_on=False,
        )

        global_mechanisms = {
            'dilution': dilution_mechanism,
            'rna_degradation': deg_mrna,
        }
        self.add_mechanisms(global_mechanisms, overwrite=None)


class TxTlDilutionMixture(Mixture):
    """In vivo TX-TL with explicit machinery, dilution, and background load.

    A mixture that models transcription and translation with explicit
    representation of RNA polymerase (RNAP), ribosomes, and RNases for in
    vivo contexts. This mixture uses Michaelis-Menten kinetics for TX-TL,
    explicitly tracking enzyme-substrate binding and catalysis. Includes
    global dilution to model cell growth effects and a background load
    component representing endogenous cellular processes that compete for
    shared machinery.

    Unlike `TxTlExtract`, this mixture includes dilution for non-DNA and
    non-machinery species. Machinery components (RNAP, ribosomes, RNases) are
    protected from dilution via the 'machinery' attribute. This model does
    not include explicit energy species.

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    rnap : str, default='RNAP'
        Name for the RNA polymerase protein species.
    ribosome : str, default='Ribo'
        Name for the ribosome protein species.
    rnaase : str, default='RNAase'
        Name for the ribonuclease protein species.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and `Mechanism` objects as values, a
        list of `Mechanism` objects, or a single `Mechanism`.
    components : list of Component or Component, optional
        Components to include in the mixture. Components are deep-copied when
        added to prevent modification of original objects.
    parameters : dict, optional
        Dictionary of parameter values. Keys follow the format
        (mechanism, part_id, param_name).
    compartment : Compartment, optional
        Default compartment for all components and species in this mixture.
    parameter_file : str, optional
        Path to a CSV or TSV file containing parameters to load.
    overwrite_parameters : bool, default=False
        If True, parameters from file/dict overwrite existing parameters.
        If False, existing parameters are preserved.
    global_mechanisms : dict, list, or GlobalMechanism, optional
        Global mechanisms that apply to all species after component
        compilation (e.g., dilution, global degradation). Can be a dict,
        list, or single `GlobalMechanism`.
    species : list of Species or Species, optional
        Additional species to add directly to the CRN without going through
        component compilation.
    initial_condition_dictionary : dict, optional
        Dictionary mapping species to initial concentration values. Deprecated
        in favor of using parameters with mechanism='initial concentration'.
    global_component_enumerators : list, optional
        List of global component enumerators for advanced component generation
        patterns (e.g., creating all pairwise interactions).
    global_recursion_depth : int, default=4
        Maximum recursion depth for global component enumeration during
        compilation.
    local_recursion_depth : int, optional
        Maximum recursion depth for local component enumeration. If None,
        defaults to `global_recursion_depth + 2`.

    Attributes
    ----------
    name : str
        Name of the mixture.
    rnap : Protein
        RNA polymerase component with 'machinery' attribute.
    ribosome : Protein
        Ribosome component with 'machinery' attribute.
    rnaase : Protein
        Ribonuclease component with 'machinery' attribute.
    compartment : Compartment or None
        Default compartment for the mixture.
    components : list of Component
        List of components in the mixture (deep copies of added components).
    mechanisms : dict
        Dictionary of default mechanisms, keyed by mechanism type (str).
    global_mechanisms : dict
        Dictionary of global mechanisms, keyed by mechanism type (str).
    parameter_database : ParameterDatabase
        Database storing all parameters for this mixture.
    added_species : list of Species
        List of species added directly to the mixture.
    global_component_enumerators : list
        List of global component enumerators.
    global_recursion_depth : int
        Recursion depth for global component enumeration.
    local_recursion_depth : int
        Recursion depth for local component enumeration.
    crn : ChemicalReactionNetwork or None
        The compiled CRN, created by calling `compile_crn`.

    See Also
    --------
    SimpleTxTlDilutionMixture : TX-TL without machinery, with dilution.
    TxTlExtract : TX-TL with machinery, without dilution.
    ExpressionDilutionMixture : Single-step expression with dilution.
    Mixture : Base class for all mixtures.

    Notes
    -----
    This mixture automatically adds the following components:

    - RNA polymerase (RNAP) with 'machinery' attribute
    - Ribosome with 'machinery' attribute
    - Ribonuclease (RNase) with 'machinery' attribute
    - Background processes DNAassembly representing cellular load

    Default mechanisms included:

    - 'transcription' : `Transcription_MM` - Michaelis-Menten transcription
      with explicit RNAP binding (DNA + RNAP <--> DNA:RNAP --> DNA + RNAP +
      mRNA)
    - 'translation' : `Translation_MM` - Michaelis-Menten translation with
      explicit ribosome binding (mRNA + Rib <--> mRNA:Rib --> mRNA + Rib +
      Protein)
    - 'rna_degradation' : `Degradation_mRNA_MM` - Global RNA degradation by
      RNase using Michaelis-Menten kinetics
    - 'catalysis' : `MichaelisMenten` - General Michaelis-Menten enzyme
      catalysis
    - 'binding' : `One_Step_Binding` - Simple multi-species binding
    - 'dilution' : `Dilution` - Global dilution mechanism (Species --> ∅)
      applied to all species except DNA and machinery

    Key features of this mixture:

    - Explicit modeling of transcription and translation machinery
    - Resource competition (genes and background processes compete for
      machinery)
    - Enzyme sequestration in complexes
    - RNA degradation dynamics
    - Global dilution modeling cell growth
    - Machinery protected from dilution
    - Background load representing cellular processes
    - Suitable for modeling in vivo gene expression with resource limits

    Background processes:

    - Implemented as a DNAassembly component ('cellular_processes')
    - Represents endogenous genes competing for machinery
    - Uses average promoter and RBS parameters
    - Creates realistic loading effects on available machinery
    - Does not model effects of loading on cell growth rate

    Common applications include:

    - In vivo gene circuit modeling with growth
    - Resource allocation in growing cells
    - Gene expression burden studies
    - Competition between heterologous and endogenous genes
    - Synthetic biology in cellular contexts

    Examples
    --------
    Create an in vivo TX-TL mixture with machinery and dilution for GFP:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.TxTlDilutionMixture(
    ...     name='cell_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/cell_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self, name='', rnap='RNAP', ribosome='Ribo', rnaase='RNAase', **kwargs
    ):
        Mixture.__init__(self, name=name, **kwargs)

        # Create Components for TxTl machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        self.rnap.add_attribute('machinery')
        self.ribosome.add_attribute('machinery')
        self.rnaase.add_attribute('machinery')

        # DNAassmbly represents background processes / loading in a cell
        background_parameters = {
            ('transcription', None, 'ku'): 50,
            ('transcription', None, 'kb'): 500,
            ('transcription', None, 'ktx'): 0.1,
            ('translation', None, 'ku'): 5,
            ('translation', None, 'kb'): 500,
            ('translation', None, 'ktl'): 0.1,
            ('rna_degradation', None, 'ku'): 50,
            ('rna_degradation', None, 'kb'): 500,
            ('rna_degradation', None, 'kdeg'): 0.1,
        }
        BackgroundProcesses = DNAassembly(
            name='cellular_processes',
            promoter='average_promoter',
            rbs='average_rbs',
            parameters=background_parameters,
        )

        default_components = [
            self.rnap,
            self.ribosome,
            self.rnaase,
            BackgroundProcesses,
        ]
        self.add_components(default_components)

        # Create TxTl Mechansisms
        mech_tx = Transcription_MM(rnap=self.rnap.get_species())
        mech_tl = Translation_MM(ribosome=self.ribosome.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

        # Create Global Dilution Mechanisms
        dilution_mechanism = Dilution(
            filter_dict={'dna': False, 'machinery': False}, default_on=True
        )
        mech_rna_deg = Degradation_mRNA_MM(nuclease=self.rnaase.get_species())

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            'dilution': dilution_mechanism,
        }

        self.add_mechanisms(default_mechanisms, overwrite=None)
