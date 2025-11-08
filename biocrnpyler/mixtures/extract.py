# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..components.basic import Metabolite, Protein
from ..components.dna.assembly import DNAassembly
from ..core.chemical_reaction_network import ChemicalReactionNetwork
from ..core.mechanism import EmptyMechanism
from ..core.mixture import Mixture
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import BasicCatalysis, MichaelisMenten
from ..mechanisms.global_mechanisms import Degradation_mRNA_MM, Dilution
from ..mechanisms.metabolite import OneStepPathway
from ..mechanisms.txtl import (
    Energy_Transcription_MM,
    Energy_Translation_MM,
    OneStepGeneExpression,
    SimpleTranscription,
    SimpleTranslation,
    Transcription_MM,
    Translation_MM,
)


class ExpressionExtract(Mixture):
    """Gene expression extract without explicit TX-TL machinery.

    A simplified mixture that models gene expression as a single direct
    reaction from DNA to protein, without explicitly representing
    transcription and translation as separate processes. This extract lumps
    transcription and translation into a single 'expression' reaction,
    eliminating intermediate mRNA species and cellular machinery (ribosomes,
    polymerases, etc.).

    This extract is appropriate for coarse-grained models where mRNA dynamics
    are negligible and computational efficiency is prioritized over
    mechanistic detail.

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and mechanism objects as values, a
        list of mechanism objects, or a single `Mechanism`.
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
    SimpleTxTlExtract : TX-TL with separate transcription and translation.
    TxTlExtract : TX-TL with explicit machinery.
    OneStepGeneExpression : Mechanism used for expression.
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

    Key features of this extract:

    - No explicit transcription or translation steps
    - No cellular machinery (RNAP, ribosomes, RNases)
    - No intermediate mRNA species
    - Simplified parameter space (single 'kexpress' rate)
    - Fast compilation and simulation

    When compiled, this extract automatically disables transcript generation
    in DNA assemblies that produce proteins, routing expression directly from
    DNA to protein.

    Common applications include:

    - High-level gene circuit modeling
    - Steady-state or quasi-steady-state analyses
    - Rapid prototyping of genetic designs
    - Models where mRNA dynamics are negligible

    Examples
    --------
    Create an expression mixture for GFP production:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.ExpressionExtract(
    ...     name='expression_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/extract_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(self, name='', **kwargs):
        # always call the superlcass Mixture.__init__(...)
        Mixture.__init__(self, name=name, **kwargs)

        # Create default Expression Mechanisms
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

    def compile_crn(self, **kwargs) -> ChemicalReactionNetwork:
        """Compile CRN with transcript generation disabled in gene expression.

        Overrides the parent `compile_crn` method to automatically disable
        transcript generation in DNA assemblies that produce proteins. This
        ensures that gene expression proceeds directly from DNA to protein
        without intermediate mRNA species.

        Parameters
        ----------
        **kwargs
            Additional keyword arguments passed to the parent Mixture
            `compile_crn <biocrnpyler.core.mixture.Mixture.compile_crn>`
            method.

        Returns
        -------
        ChemicalReactionNetwork
            Compiled chemical reaction network with expression reactions.

        Notes
        -----
        This method automatically modifies DNA assemblies before compilation:

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
                # Only turn off transcription for an Assembly that
                # makes a Protein.  Some assemblies might only make
                # RNA!
                if component.protein is not None:
                    # This will turn off transcription and set
                    # Promoter.transcript = False Mechanisms that
                    # recieve no transcript but a protein will use the
                    # protein instead.
                    component.update_transcript(False)

        # Call the superclass function
        return Mixture.compile_crn(self, **kwargs)


class SimpleTxTlExtract(Mixture):
    """TX-TL extract with simple transcription and translation mechanisms.

    A mixture that models transcription and translation as separate catalytic
    reactions without explicitly representing cellular machinery (RNAP,
    ribosomes, RNases). This extract uses simple mass-action kinetics where
    DNA and mRNA act as catalysts for transcript and protein production,
    respectively. Unlike `ExpressionExtract`, this mixture includes explicit
    mRNA species and separate TX-TL steps. Unlike `TxTlExtract`, it does not
    model enzyme binding or resource competition.

    This extract includes global RNA degradation via dilution.

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and mechanism objects as values, a
        list of mechanism objects, or a single `Mechanism`.
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
    ExpressionExtract : Single-step expression without transcripts.
    TxTlExtract : TX-TL with explicit machinery.
    SimpleTranscription : Mechanism used for transcription.
    SimpleTranslation : Mechanism used for translation.
    Mixture : Base class for all mixtures.

    Notes
    -----
    Default mechanisms included:

    - 'transcription' : `SimpleTranscription` - Simple catalytic
      transcription (DNA --> DNA + mRNA) without explicit RNAP binding
    - 'translation' : `SimpleTranslation` - Simple catalytic
      translation (mRNA --> mRNA + Protein) without explicit ribosome binding
    - 'rna_degradation' : `Dilution` - Global RNA degradation mechanism
      (mRNA --> ∅) applied to all RNA species
    - 'catalysis' : `BasicCatalysis` - Simple catalytic reactions without
      explicit enzyme binding
    - 'binding' : `One_Step_Binding` - Simple multi-species binding

    Key features of this extract:

    - Explicit transcription and translation steps
    - Intermediate mRNA species
    - Simple mass-action kinetics (no enzyme binding)
    - No cellular machinery (RNAP, ribosomes)
    - Global RNA degradation
    - Faster simulation than Michaelis-Menten models

    Common applications include:

    - Gene circuit modeling with explicit TX-TL
    - Models where machinery is not limiting
    - Constitutive or weakly regulated promoters
    - Rapid prototyping with mRNA dynamics

    Examples
    --------
    Create a simple TX-TL mixture for GFP expression:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.SimpleTxTlExtract(
    ...     name='simple_txtl_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/extract_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(self, name='', **kwargs):
        # Always call the superlcass Mixture.__init__(...)
        Mixture.__init__(self, name=name, **kwargs)

        # TxTl Mechanisms
        mech_tx = SimpleTranscription()
        mech_tl = SimpleTranslation()
        mech_cat = BasicCatalysis()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=False)

        # global mechanisms for dilution and rna degradation
        mech_rna_deg_global = Dilution(
            name='rna_degradation',
            filter_dict={'rna': True},
            default_on=False,
        )
        global_mechanisms = {'rna_degradation': mech_rna_deg_global}
        self.add_mechanisms(global_mechanisms, overwrite=None)


class TxTlExtract(Mixture):
    """TX-TL extract with explicit transcription and translation machinery.

    A mixture that models transcription and translation with explicit
    representation of RNA polymerase (RNAP), ribosomes, and RNases. This
    extract uses Michaelis-Menten kinetics for transcription and translation,
    explicitly tracking enzyme-substrate binding and catalysis. Unlike
    `SimpleTxTlExtract`, this mixture models resource competition and enzyme
    sequestration effects.

    This model does not include explicit energy species. For energy-aware
    modeling, use `EnergyTxTlExtract`.

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
        mechanism types (str) as keys and mechanism objects as values, a
        list of mechanism objects, or a single `Mechanism`.
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
    crn : ChemicalReactionNetwork or None
        The compiled CRN, created by calling `compile_crn`.
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
    rnap : Protein
        RNA polymerase component.
    ribosome : Protein
        Ribosome component.
    rnaase : Protein
        Ribonuclease component.

    See Also
    --------
    SimpleTxTlExtract : TX-TL without explicit machinery.
    EnergyTxTlExtract : TX-TL with explicit energy consumption.
    ExpressionExtract : Combined TX-TL without transcripts.
    Mixture : Base class for all mixtures.

    Notes
    -----
    This mixture automatically adds the following components:

    - RNA polymerase (RNAP)
    - Ribosome
    - Ribonuclease (RNase)

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

    Key features of this mixture:

    - Explicit modeling of transcription and translation machinery
    - Resource competition effects (multiple genes compete for RNAP)
    - Enzyme sequestration in complexes
    - RNA degradation dynamics
    - Suitable for modeling TX-TL systems with limited machinery

    Common applications include:

    - Cell-free TX-TL systems
    - Resource allocation in gene circuits
    - Gene expression burden studies
    - Synthetic biology prototyping

    Examples
    --------
    Create a TX-TL mixture for GFP expression:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.TxTlExtract(
    ...     name='txtl_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/extract_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self, name='', rnap='RNAP', ribosome='Ribo', rnaase='RNAase', **kwargs
    ):
        # Always call the superlcass Mixture.__init__(...)
        Mixture.__init__(self, name=name, **kwargs)

        # create default Components to represent cellular machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        default_components = [self.rnap, self.ribosome, self.rnaase]
        self.add_components(default_components)

        # Create default TxTl Mechanisms
        mech_tx = Transcription_MM(rnap=self.rnap.get_species())
        mech_tl = Translation_MM(ribosome=self.ribosome.get_species())
        mech_rna_deg = Degradation_mRNA_MM(nuclease=self.rnaase.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)


class EnergyTxTlExtract(Mixture):
    """TX-TL cell extract with explicit machinery and energy consumption.

    A mixture that models transcription and translation with explicit
    representation of RNA polymerase (RNAP), ribosomes, RNases, and energy
    carrier molecules. This extract uses Michaelis-Menten kinetics with
    length-dependent fuel consumption to model realistic TX-TL energetics.
    Unlike `TxTlExtract`, this mixture explicitly tracks NTPs, amino acids,
    and fuel species (e.g., 3PGA for NTP regeneration).

    Energy usage for transcription and translation is length-dependent,
    reflecting the stoichiometric consumption of NTPs and amino acids during
    biopolymer synthesis.

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
    ntps : str, default='NTPs'
        Name for the nucleotide triphosphate species (lumped NTPs).
    ndps : str, default='NDPs'
        Name for the nucleotide diphosphate species (lumped NDPs).
    amino_acids : str, default='amino_acids'
        Name for the amino acid species (lumped amino acids).
    fuel : str, default='Fuel_3PGA'
        Name for the fuel species used for NTP regeneration (e.g., 3PGA).
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and mechanism objects as values, a
        list of mechanism objects, or a single `Mechanism`.
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
        RNA polymerase component.
    ribosome : Protein
        Ribosome component.
    rnaase : Protein
        Ribonuclease component.
    amino_acids : Metabolite
        Amino acid metabolite component.
    fuel : Metabolite
        Fuel metabolite component for ATP regeneration.
    ndps : Metabolite
        Nucleotide diphosphate metabolite component.
    ntps : Metabolite
        Nucleotide triphosphate metabolite component with fuel-dependent
        regeneration.
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
    TxTlExtract : TX-TL with machinery but no energy.
    SimpleTxTlExtract : TX-TL without machinery or energy.
    Energy_Transcription_MM : Mechanism for energy-consuming transcription.
    Energy_Translation_MM : Mechanism for energy-consuming translation.
    Mixture : Base class for all mixtures.

    Notes
    -----
    This mixture automatically adds the following components:

    - RNA polymerase (RNAP)
    - Ribosome
    - Ribonuclease (RNase)
    - Amino acids (lumped)
    - NTPs (nucleotide triphosphates, lumped)
    - NDPs (nucleotide diphosphates, lumped)
    - Fuel (e.g., 3PGA for ATP regeneration)

    Default mechanisms included:

    - 'transcription' : `Energy_Transcription_MM` - Michaelis-Menten
      transcription with length-dependent NTP consumption (DNA + RNAP <-->
      DNA:RNAP; NTP + DNA:RNAP --> DNA + RNAP + mRNA + NDP)
    - 'translation' : `Energy_Translation_MM` - Michaelis-Menten translation
      with length-dependent amino acid and NTP consumption (mRNA + Rib <-->
      mRNA:Rib; AA + NTP + mRNA:Rib --> mRNA + Rib + Protein + NDP)
    - 'rna_degradation' : `Degradation_mRNA_MM` - Global RNA degradation by
      RNase using Michaelis-Menten kinetics
    - 'catalysis' : `MichaelisMenten` - General Michaelis-Menten enzyme
      catalysis
    - 'binding' : `One_Step_Binding` - Simple multi-species binding
    - 'pathway' : `OneStepPathway` - Metabolite conversion (added to NTPs
      and fuel components)

    Key features of this mixture:

    - Explicit modeling of transcription and translation machinery
    - Length-dependent energy consumption
    - NTP regeneration from fuel species
    - Resource competition and depletion effects
    - Realistic modeling of TX-TL resource limits
    - Energy-dependent expression dynamics

    Energy model details:

    - Transcription consumes L NTPs per mRNA of length L
    - Translation consumes L amino acids and 4L NTPs per protein of length L
    - Fuel species regenerates NTPs from NDPs
    - Different nucleotides and amino acids are lumped together

    Common applications include:

    - Cell-free TX-TL systems with limited resources
    - Models of energy-limited gene expression
    - Resource allocation and burden studies
    - TX-TL system optimization
    - Metabolic coupling with gene expression

    Examples
    --------
    Create an energy-aware TX-TL mixture for GFP expression:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.EnergyTxTlExtract(
    ...     name='energy_txtl_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/extract_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        name='',
        rnap='RNAP',
        ribosome='Ribo',
        rnaase='RNAase',
        ntps='NTPs',
        ndps='NDPs',
        amino_acids='amino_acids',
        fuel='Fuel_3PGA',
        **kwargs,
    ):
        Mixture.__init__(self, name=name, **kwargs)

        # create default Components to represent cellular machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)
        self.amino_acids = Metabolite(amino_acids)
        # fuel is degraded into things other than ATP as well
        self.fuel = Metabolite(fuel)
        self.ndps = Metabolite(ndps)  # NDPs
        self.ntps = Metabolite(
            ntps, precursors=[self.fuel, self.ndps], products=[self.ndps]
        )  # fuel becomes ATP, and ATP is degraded

        # These mechanisms are Component specific and only added to
        # the NTPs metabolite
        mech_pathway = OneStepPathway()
        self.ntps.add_mechanisms(mech_pathway, overwrite=None)
        self.fuel.add_mechanisms(mech_pathway, overwrite=None)

        default_components = [
            self.rnap,
            self.ribosome,
            self.rnaase,
            self.amino_acids,
            self.ntps,
            self.fuel,
        ]
        self.add_components(default_components)

        # Create default TxTl Mechanisms
        mech_tx = Energy_Transcription_MM(
            rnap=self.rnap.get_species(),
            fuels=[self.ntps.get_species()],
            wastes=[],
        )
        mech_tl = Energy_Translation_MM(
            ribosome=self.ribosome.get_species(),
            fuels=4 * [self.ntps.get_species()]
            + [self.amino_acids.get_species()],
            wastes=4 * [self.ndps.get_species()],
        )
        mech_rna_deg = Degradation_mRNA_MM(nuclease=self.rnaase.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)
