# pure.py - mixture model for PURE
# RMM, 20 Sep 2025

from ..components.basic import Metabolite, Protein
from ..core.mixture import Mixture
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import MichaelisMenten
from ..mechanisms.global_mechanisms import Degradation_mRNA_MM
from ..mechanisms.txtl import Energy_Transcription_MM, Energy_Translation_MM


class BasicPURE(Mixture):
    """Reconstituted protein synthesis system with resource limits.

    This model includes energy carrier molecules in the form of NTPs, amino
    acids, and a fuel species (such as ATP) used for transcription,
    translation, and other core mechanisms.  This model is equivalent to
    `EnergyTxTlExtract`, but without a fuel generation mechanism.  Amino
    acids and nucleotides are lumped together into a single meta-species.

    Note that fuel is modeled as a separate molecule so if the default
    'ATP' is used, it is separate from the other nucleotides ('NTPs').

    Energy usage for transcription and translation is length dependent.

    """

    def __init__(
        self,
        name='PURE',
        rnap='RNAP',
        ribosome='Ribo',
        rnaase='RNAase',
        ntps='NTPs',
        ndps='NDPs',
        amino_acids='AAs',
        fuel='ATP',
        parameter_file='mixtures/pure_parameters.tsv',
        **kwargs,
    ):
        """Initialize the PURE mixture.

        :param name: name of the mixture
        :param rnap: name of the RNA polymerase, default: RNAP
        :param ribosome: name of the ribosome, default: Ribo
        :param rnaase: name of the Ribonuclease, default: RNAase
        :param ntps: name of the nucleotide fuel source (eg ATP + GTP etc),
            default: NTP
        :param amino_acids: name of the amino acids species, default:
            amino_acids
        :param fuel: name of the energy carrier species
        :param parameter_file: file containing default parameter values
        :param parameter: dictionary with parameter values
        :param kwargs: keywords passed into the parent Class (Mixture)

        """
        Mixture.__init__(
            self, name=name, parameter_file=parameter_file, **kwargs
        )

        # create default Components to represent cellular machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)
        self.ntps = Metabolite(ntps)
        self.amino_acids = Metabolite(amino_acids)
        self.fuel = Metabolite(fuel)

        default_components = [
            self.rnap,
            self.ribosome,
            self.rnaase,
            self.amino_acids,
            self.ntps,
            self.fuel,
        ]
        self.add_components(default_components)

        # Create default TX-TL Mechanisms
        mech_tx = Energy_Transcription_MM(
            rnap=self.rnap.get_species(),
            fuels=[self.fuel.get_species()]  # TODO: one ATP per bp
            + [self.ntps.get_species()],
            wastes=[],
        )
        mech_tl = Energy_Translation_MM(
            ribosome=self.ribosome.get_species(),
            fuels=4 * [self.fuel.get_species()]  # TODO: why 4 ATP per AA?
            + [self.amino_acids.get_species()],
            wastes=[],
        )
        mech_rna_deg = Degradation_mRNA_MM(
            nuclease=self.rnaase.get_species()
        )  # TODO: add fuel usage?
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms)
