# pure.py - mixture model for PURE
# RMM, 20 Sep 2025

from ..components.basic import Metabolite, Protein
from ..core.mixture import Mixture
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import MichaelisMenten
from ..mechanisms.global_mechanisms import Degradation_mRNA_MM
from ..mechanisms.txtl import Energy_Transcription_MM, Energy_Translation_MM


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
        Mixture.__init__(
            self, name=name, parameter_file=parameter_file, **kwargs
        )

        # create default Components to represent cellular machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.ntps = Metabolite(ntps)
        self.amino_acids = Metabolite(amino_acids)
        self.fuel = Metabolite(fuel)

        default_components = [
            self.rnap,
            self.ribosome,
            self.amino_acids,
            self.ntps,
            self.fuel,
        ]
        self.add_components(default_components)

        # Create default TX-TL Mechanisms
        mech_tx = Energy_Transcription_MM(
            rnap=self.rnap.get_species(),
            fuels=[self.ntps.get_species()],
            wastes=[],
        )
        mech_tl = Energy_Translation_MM(
            ribosome=self.ribosome.get_species(),
            fuels=4 * [self.fuel.get_species()]
            + [self.amino_acids.get_species()],
            wastes=[],
        )
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)
