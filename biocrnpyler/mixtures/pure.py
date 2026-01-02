# pure.py - mixture model for PURE
# RMM, 20 Sep 2025

from ..components.basic import Metabolite, Protein
from ..core.mixture import Mixture
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import MichaelisMenten
from ..mechanisms.txtl import (
    Energy_Transcription_MM,
    Energy_Translation_MM,
    SimpleTranscription,
    SimpleTranslation,
    Transcription_MM,
    Translation_MM,
)


class PURE(Mixture):
    """PURE cell-free protein synthesis system with customizable mechanisms.

    A mixture that models the PURE (Protein synthesis Using Recombinant
    Elements) reconstituted cell-free transcription-translation system with
    customizable mechanisms.

    Parameters
    ----------
    name : str, default='PURE'
        Name identifier for the mixture.
    include_machinery : bool, default=True
        Include components and mechanisms for RNAP and ribosomes.
    include_resources : bool, default=True
        Include components and mechanisms for NTP and amino acid utilization.
        Requires `include_machinery` to be true.
    include_fuel : bool, default=True
        Include fuel component and energy utilization mechanism. Requires that
        `include_machinery` and `include_resources` be True.
    rnap : str, default='RNAP'
        Name for the RNA polymerase protein species.
    ribosome : str, default='Ribo'
        Name for the ribosome protein species.
    ntps : str, default='NTPs'
        Name for the nucleotide triphosphate species (lumped NTPs excluding
        ATP).
    ndps : str, default='NDPs'
        Name for the nucleotide diphosphate species (lumped NDPs).
    amino_acids : str, default='AAs'
        Name for the amino acid species (lumped amino acids).
    fuel : str, default='ATP'
        Name for the primary energy carrier species (ATP).
    parameter_file : str, default='mixtures/pure_parameters.tsv'
        Path to file containing default parameter values for the PURE
        system.
    **kwargs
        Additional keyword arguments passed to the parent Mixture class.

    Attributes
    ----------
    rnap : Protein
        RNA polymerase component.
    ribosome : Protein
        Ribosome component.
    ntps : Metabolite
        Nucleotide triphosphate metabolite component (excluding ATP).
    amino_acids : Metabolite
        Amino acid metabolite component.
    fuel : Metabolite
        Fuel metabolite component (ATP).
    name : str
        Name of the mixture.

    See Also
    --------
    EnergyTxTlExtract : TX-TL with fuel regeneration.
    TxTlExtract : TX-TL with machinery but no energy.
    Energy_Transcription_MM : Mechanism for energy-consuming transcription.
    Energy_Translation_MM : Mechanism for energy-consuming translation.
    Mixture : Base class for all mixtures.

    Notes
    -----
    This mixture automatically adds the following components:

    - RNA polymerase (RNAP)
    - Ribosome
    - Amino acids (lumped)
    - NTPs (nucleotide triphosphates excluding ATP, lumped)
    - NDPs (nucleotide diphosphates, lumped)
    - Fuel (ATP for energy)

    Default mechanisms included:

    - 'transcription' : `Energy_Transcription_MM` - Michaelis-Menten
      transcription with length-dependent ATP and NTP consumption
    - 'translation' : `Energy_Translation_MM` - Michaelis-Menten translation
      with length-dependent amino acid and ATP consumption
    - 'catalysis' : `MichaelisMenten` - General Michaelis-Menten enzyme
      catalysis for user-defined enzymatic reactions
    - 'binding' : `One_Step_Binding` - Simple multi-species binding for
      forming complexes

    Key features of this mixture:

    - Explicit modeling of PURE system components
    - Length-dependent energy consumption (realistic stoichiometry)
    - No fuel regeneration mechanisms (finite resource pool)
    - Resource competition effects (genes compete for RNAP and ribosomes)
    - Resource depletion dynamics (ATP, NTPs, amino acids deplete)
    - Enzyme sequestration in complexes
    - Separate tracking of ATP vs other NTPs
    - Suitable for modeling batch-mode PURE reactions

    Energy model details:

    - Transcription: Consumes L NTPs and L ATPs per mRNA of length L
    - Translation: Consumes L amino acids and 4L ATPs per protein of length
      L (4 ATPs per amino acid reflect GTP hydrolysis during elongation)
    - No regeneration: ATP, NTPs, and amino acids are consumed but not
      regenerated
    - Energy depletion: Expression stops when resources are exhausted
    - Length parameter L: Represents gene/protein length in appropriate
      units
    - Lumped species: Different nucleotides lumped into NTPs, different
      amino acids lumped into single species
    - Separate ATP: ATP tracked separately from other NTPs for independent
      energy accounting

    Differences from `EnergyTxTlExtract`:

    - No fuel regeneration pathway (no NTP regeneration from 3PGA or other
      fuel sources)
    - ATP modeled as separate fuel species rather than included in NTPs
    - Default parameter file points to PURE-specific parameters
    - Intended for modeling finite-resource batch reactions
    - More realistic for in vitro PURE systems

    Common applications include:

    - PURE cell-free TX-TL systems
    - Resource-limited gene expression modeling
    - TX-TL system optimization with fixed resource budgets
    - Batch mode TX-TL reactions
    - Energy budget and resource allocation studies
    - Multi-gene expression burden analysis
    - In vitro synthetic biology applications

    Examples
    --------
    Create a PURE mixture for GFP expression:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.BasicPURE(
    ...     name='pure_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/pure_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        include_machinery=True,
        include_resources=True,
        include_fuel=True,
        name='PURE',
        rnap='RNAP',
        ribosome='Ribo',
        ntps='NTPs',
        ndps='NDPs',
        amino_acids='AAs',
        fuel='ATP',
        parameter_file='mixtures/pure_parameters.tsv',
        **kwargs,
    ):
        Mixture.__init__(
            self, name=name, parameter_file=parameter_file, **kwargs
        )

        # Start with basic mechansims for transcription and translation
        default_mechanisms = {
            'transcription': SimpleTranscription(),
            'translation': SimpleTranslation(),
        }

        if include_fuel:
            if not include_resources:
                raise ValueError("include_fuel requires include_resources")
            self.fuel = Metabolite(fuel)
            self.add_components([self.fuel])
        else:
            self.fuel = None

        if include_machinery:
            self.rnap = Protein(rnap)
            self.ribosome = Protein(ribosome)
            self.add_components([self.rnap, self.ribosome])

            default_mechanisms |= {
                'transcription': Transcription_MM(rnap=self.rnap.species),
                'translation': Translation_MM(ribosome=self.ribosome.species),
            }

        if include_resources:
            if not include_machinery:
                raise ValueError(
                    "include_resources requires include_machinery"
                )

            # create default Components to represent cellular machinery
            self.ntps = Metabolite(ntps)
            self.amino_acids = Metabolite(amino_acids)
            self.add_components([self.amino_acids, self.ntps])

            # Create default TX-TL Mechanisms
            default_mechanisms['transcription'] = Energy_Transcription_MM(
                rnap=self.rnap.get_species(),
                fuels=[self.ntps.get_species()],
                wastes=[],
            )
            if include_fuel:
                default_mechanisms['translation'] = Energy_Translation_MM(
                    ribosome=self.ribosome.get_species(),
                    fuels=4 * [self.fuel.get_species()]
                    + [self.amino_acids.get_species()],
                    wastes=[],
                )
            else:
                default_mechanisms['translation'] = Energy_Translation_MM(
                    ribosome=self.ribosome.get_species(),
                    fuels=[self.amino_acids.get_species()],
                    wastes=[],
                )
        self.add_mechanisms(mechanisms=default_mechanisms, overwrite=False)


class BasicPURE(Mixture):
    """PURE cell-free protein synthesis system with energy consumption.

    A mixture that models the PURE (Protein synthesis Using Recombinant
    Elements) reconstituted cell-free transcription-translation system with
    explicit representation of RNA polymerase (RNAP), ribosomes, and
    energy carrier molecules. This extract uses Michaelis-Menten kinetics
    with length-dependent fuel consumption to model realistic TX-TL
    energetics.

    Unlike `EnergyTxTlExtract`, this mixture does not include fuel
    regeneration mechanisms. Energy carriers (ATP, NTPs, amino acids) are
    consumed but not regenerated, making this suitable for modeling
    resource-limited PURE systems. Different amino acids and nucleotides are
    lumped into single meta-species for simplicity.

    Note that fuel (default 'ATP') is modeled as a separate molecule from
    other nucleotides ('NTPs'), allowing independent tracking of energy
    consumption.

    Energy usage for transcription and translation is length-dependent,
    reflecting stoichiometric consumption during biopolymer synthesis.

    Parameters
    ----------
    name : str, default='PURE'
        Name identifier for the mixture.
    rnap : str, default='RNAP'
        Name for the RNA polymerase protein species.
    ribosome : str, default='Ribo'
        Name for the ribosome protein species.
    ntps : str, default='NTPs'
        Name for the nucleotide triphosphate species (lumped NTPs excluding
        ATP).
    ndps : str, default='NDPs'
        Name for the nucleotide diphosphate species (lumped NDPs).
    amino_acids : str, default='AAs'
        Name for the amino acid species (lumped amino acids).
    fuel : str, default='ATP'
        Name for the primary energy carrier species (ATP).
    parameter_file : str, default='mixtures/pure_parameters.tsv'
        Path to file containing default parameter values for the PURE
        system.
    **kwargs
        Additional keyword arguments passed to the parent Mixture class.

    Attributes
    ----------
    rnap : Protein
        RNA polymerase component.
    ribosome : Protein
        Ribosome component.
    ntps : Metabolite
        Nucleotide triphosphate metabolite component (excluding ATP).
    amino_acids : Metabolite
        Amino acid metabolite component.
    fuel : Metabolite
        Fuel metabolite component (ATP).
    name : str
        Name of the mixture.

    See Also
    --------
    EnergyTxTlExtract : TX-TL with fuel regeneration.
    TxTlExtract : TX-TL with machinery but no energy.
    Energy_Transcription_MM : Mechanism for energy-consuming transcription.
    Energy_Translation_MM : Mechanism for energy-consuming translation.
    Mixture : Base class for all mixtures.

    Notes
    -----
    This mixture automatically adds the following components:

    - RNA polymerase (RNAP)
    - Ribosome
    - Amino acids (lumped)
    - NTPs (nucleotide triphosphates excluding ATP, lumped)
    - NDPs (nucleotide diphosphates, lumped)
    - Fuel (ATP for energy)

    Default mechanisms included:

    - 'transcription' : `Energy_Transcription_MM` - Michaelis-Menten
      transcription with length-dependent ATP and NTP consumption
    - 'translation' : `Energy_Translation_MM` - Michaelis-Menten translation
      with length-dependent amino acid and ATP consumption
    - 'catalysis' : `MichaelisMenten` - General Michaelis-Menten enzyme
      catalysis for user-defined enzymatic reactions
    - 'binding' : `One_Step_Binding` - Simple multi-species binding for
      forming complexes

    Key features of this mixture:

    - Explicit modeling of PURE system components
    - Length-dependent energy consumption (realistic stoichiometry)
    - No fuel regeneration mechanisms (finite resource pool)
    - Resource competition effects (genes compete for RNAP and ribosomes)
    - Resource depletion dynamics (ATP, NTPs, amino acids deplete)
    - Enzyme sequestration in complexes
    - Separate tracking of ATP vs other NTPs
    - Suitable for modeling batch-mode PURE reactions

    Energy model details:

    - Transcription: Consumes L NTPs and L ATPs per mRNA of length L
    - Translation: Consumes L amino acids and 4L ATPs per protein of length
      L (4 ATPs per amino acid reflect GTP hydrolysis during elongation)
    - No regeneration: ATP, NTPs, and amino acids are consumed but not
      regenerated
    - Energy depletion: Expression stops when resources are exhausted
    - Length parameter L: Represents gene/protein length in appropriate
      units
    - Lumped species: Different nucleotides lumped into NTPs, different
      amino acids lumped into single species
    - Separate ATP: ATP tracked separately from other NTPs for independent
      energy accounting

    Differences from `EnergyTxTlExtract`:

    - No fuel regeneration pathway (no NTP regeneration from 3PGA or other
      fuel sources)
    - ATP modeled as separate fuel species rather than included in NTPs
    - Default parameter file points to PURE-specific parameters
    - Intended for modeling finite-resource batch reactions
    - More realistic for in vitro PURE systems

    Common applications include:

    - PURE cell-free TX-TL systems
    - Resource-limited gene expression modeling
    - TX-TL system optimization with fixed resource budgets
    - Batch mode TX-TL reactions
    - Energy budget and resource allocation studies
    - Multi-gene expression burden analysis
    - In vitro synthetic biology applications

    Examples
    --------
    Create a PURE mixture for GFP expression:

    >>> gfp_gene = bcp.DNAassembly(
    ...     name='gfp_construct',
    ...     promoter='pconst',
    ...     rbs='bcd2',
    ...     transcript='gfp_mrna',
    ...     protein='GFP'
    ... )
    >>> mixture = bcp.BasicPURE(
    ...     name='pure_mixture',
    ...     components=[gfp_gene],
    ...     parameter_file='mixtures/pure_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        name='PURE',
        rnap='RNAP',
        ribosome='Ribo',
        ntps='NTPs',
        ndps='NDPs',
        amino_acids='AAs',
        fuel='ATP',
        parameter_file='mixtures/pure_parameters.tsv',
        **kwargs,
    ):
        PURE.__init__(
            self,
            name=name,
            parameter_file=parameter_file,
            include_machinery=True,
            include_resources=True,
            include_fuel=True,
            **kwargs,
        )
