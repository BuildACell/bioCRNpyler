# pure.py - mixture model for PURE
# RMM, 20 Sep 2025

from warnings import warn

from ..components.basic import Metabolite, Protein
from ..core.mixture import Mixture
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
    include_energy : bool, default=True
        Include the energy carrier species and their utilization during
        translation.  Requires that `include_machinery` and
        `include_resources` be True.  When False, the `atp`, `adp`, and
        `fuel` components are set to None.
    rnap : str, default='RNAP'
        Name for the RNA polymerase protein species.
    ribosome : str, default='Ribo'
        Name for the ribosome protein species.
    ntps : str, default='NTPs'
        Name for the nucleotide triphosphate species, lumping the four
        nucleotides that are polymerized during transcription.  ATP is
        counted here in its role as a monomer; `atp` tracks its separate
        role as an energy carrier.
    amino_acids : str, default='AAs'
        Name for the amino acid species (lumped amino acids).
    atp : str, default='ATP'
        Name for the primary energy carrier species, consumed during
        translation.
    adp : str, default='ADP'
        Name for the spent energy carrier species produced from `atp`.
    fuel : str or None, default='Fuel_CP'
        Name for the secondary fuel species that, together with `adp`,
        regenerates `atp`.  If None, no fuel species is created and `atp`
        is consumed without being regenerated.
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
        Nucleotide triphosphate metabolite component.
    amino_acids : Metabolite
        Amino acid metabolite component.
    fuel : Metabolite
        Secondary fuel metabolite component, or None if there is no
        regeneration.
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
    - NTPs (nucleotide triphosphates, lumped)
    - ATP and ADP (the energy carrier and its spent form)
    - Fuel (regenerates ATP, unless `fuel` is None)

    Default mechanisms included:

    - 'transcription' : `Energy_Transcription_MM` - Michaelis-Menten
      transcription with length-dependent ATP and NTP consumption
    - 'translation' : `Energy_Translation_MM` - Michaelis-Menten translation
      with length-dependent amino acid and ATP consumption
    - 'catalysis' : `MichaelisMenten` - General Michaelis-Menten enzyme
      catalysis for user-defined enzymatic reactions
    - 'binding' : `One_Step_Binding` - Simple multi-species binding for
      forming complexes

    Unlike the extract mixtures, no 'rna_degradation' mechanism is included.
    PURE is reconstituted from purified components and carries no
    ribonucleases, and the detailed PURE reaction network of Jurado, Pandey
    and Murray (2023) has no mRNA degradation step either; their MGapt
    measurements show transcript continuing to accumulate past two hours.
    Apparent decay of the MGapt signal in commercial PURE is a property of
    the dye chemistry rather than of the mRNA (Jurado and Murray, 2024).

    A consequence is that mRNA accumulates for the length of the
    simulation and progressively binds ribosomes.  That is intended: once
    ATP is exhausted the ribosomes stall on transcript they can no longer
    translate.

    Key features of this mixture:

    - Explicit modeling of PURE system components
    - Length-dependent energy consumption (realistic stoichiometry)
    - ATP regeneration from a finite pool of fuel, or none at all if
      `fuel` is None
    - Resource competition effects (genes compete for RNAP and ribosomes)
    - Resource depletion dynamics (ATP, NTPs, amino acids deplete)
    - Enzyme sequestration in complexes
    - Separate tracking of ATP as an energy carrier and as a nucleotide
    - Suitable for modeling batch-mode PURE reactions

    Parameters are looked up by species name, so when renaming a species
    the values in the default `parameter_file` will no longer match, and
    replacements have to be supplied for the new name.  For most species
    it is only the 'initial concentration' parameter that is affected,
    but in some cases, such as `atp`, there can also be reaction rates.
    For example::

        PURE(
            ribosome='Ribosome',
            atp='GTP',
            parameters={
                ('initial concentration', None, 'Ribosome'): 3,
                ('initial concentration', None, 'GTP'): 1000,
                ('one_step_pathway', 'GTP_production', 'k'): 0.02,
                ('one_step_pathway', 'GTP_degradation', 'k'): 0.0000177,
            },
        )

    Energy model details:

    - Transcription: Consumes L NTPs per mRNA of length L.  The
      nucleotides are incorporated as monomers, so no energy carrier is
      spent beyond the nucleotides themselves
    - Translation: Consumes L amino acids and 4L ATPs per protein of
      length L.  The four high-energy phosphates per residue are two for
      charging the tRNA and two for the GTP hydrolyzed during elongation,
      all counted here as ATP
    - Energy depletion: Expression stops when resources are exhausted
    - Length parameter L: Represents gene/protein length in appropriate
      units
    - Lumped species: Different nucleotides lumped into NTPs, different
      amino acids lumped into single species
    - Separate ATP: ATP appears in two roles, as one of the nucleotides
      polymerized during transcription, counted in NTPs, and as the
      energy carrier spent during translation, counted in `atp`

    Differences from `EnergyTxTlExtract`:

    - Regenerates ATP from creatine phosphate rather than from 3PGA, and
      only when `fuel` is given
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
        include_energy=True,
        name='PURE',
        rnap='RNAP',
        ribosome='Ribo',
        ntps='NTPs',
        amino_acids='AAs',
        atp='ATP',
        adp='ADP',
        fuel='Fuel_CP',
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

        if include_energy:
            if not include_resources:
                raise ValueError("include_energy requires include_resources")
            self.adp = Metabolite(adp)
            if fuel is None:
                # ATP is consumed but never replenished
                self.fuel = None
                self.atp = Metabolite(atp)
                self.add_components([self.atp, self.adp])
            else:
                self.fuel = Metabolite(fuel)
                self.atp = Metabolite(
                    atp, precursors=[self.fuel, self.adp], products=[self.adp]
                )  # fuel becomes ATP, and ATP is degraded
                self.add_components([self.atp])  # includes ADP, fuel
        else:
            self.atp = None
            self.adp = None
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
            if include_energy:
                default_mechanisms['translation'] = Energy_Translation_MM(
                    ribosome=self.ribosome.get_species(),
                    fuels=4 * [self.atp.get_species()]
                    + [self.amino_acids.get_species()],
                    wastes=4 * [self.adp.get_species()],
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

    .. deprecated:: 1.3.0
        Use `PURE` instead, which covers this configuration and others.
        `BasicPURE` is equivalent to `PURE` with `include_machinery`,
        `include_resources`, and `include_energy` all set to True.

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

    Note that the energy carrier (default 'ATP') is modeled as a separate
    molecule from other nucleotides ('NTPs'), allowing independent tracking
    of energy consumption.

    Energy usage for transcription and translation is length-dependent,
    reflecting stoichiometric consumption during biopolymer synthesis.

    Parameters
    ----------
    name : str, default='PURE'
        Name identifier for the mixture.
    fuel : str or None, optional
        Name for the energy carrier species, retained for compatibility
        with earlier versions of this mixture.  Equivalent to the `atp`
        argument of `PURE`; giving both is an error.
    parameter_file : str, default='mixtures/pure_parameters.tsv'
        Path to file containing default parameter values for the PURE
        system.
    **kwargs
        Additional keyword arguments passed to `PURE`, including the
        species names `rnap`, `ribosome`, `ntps`, `amino_acids`, `atp`,
        and `adp`.

    Attributes
    ----------
    rnap : Protein
        RNA polymerase component.
    ribosome : Protein
        Ribosome component.
    ntps : Metabolite
        Nucleotide triphosphate metabolite component.
    amino_acids : Metabolite
        Amino acid metabolite component.
    atp : Metabolite
        Energy carrier metabolite component.
    name : str
        Name of the mixture.

    See Also
    --------
    EnergyTxTlExtract : TX-TL with fuel regeneration.
    TxTlExtract : TX-TL with machinery but no energy.
    Energy_Transcription_MM : Mechanism for energy-consuming transcription.
    Energy_Translation_MM : Mechanism for energy-consuming translation.
    Mixture : Base class for all mixtures.

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

    Notes
    -----
    The `fuel` argument does not mean the same thing here as it does in
    `PURE`.  In this mixture it names the energy carrier, as it did before
    `PURE` was introduced, and is an alias for the `atp` argument of
    `PURE`.  In `PURE` it names the species that regenerates the carrier,
    which this mixture does not have.

    """

    def __init__(
        self,
        name='PURE',
        fuel=None,
        parameter_file='mixtures/pure_parameters.tsv',
        **kwargs,
    ):
        warn(
            "BasicPURE is deprecated; use PURE instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        # Here fuel names the energy carrier, which PURE calls atp.  Note
        # that fuel=None means "not given" in this signature, but means
        # "no regeneration" in the call to PURE below.
        if fuel is not None:
            if 'atp' in kwargs:
                raise ValueError(
                    "BasicPURE accepts fuel or atp, not both; they name "
                    "the same species"
                )
            kwargs['atp'] = fuel
        PURE.__init__(
            self,
            name=name,
            parameter_file=parameter_file,
            include_machinery=True,
            include_resources=True,
            include_energy=True,
            fuel=None,
            **kwargs,
        )
