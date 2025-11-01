# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..components.basic import Metabolite, Protein
from ..components.dna.assembly import DNAassembly
from ..core.chemical_reaction_network import ChemicalReactionNetwork
from ..core.mechanism import EmptyMechanism
from ..core.mixture import Mixture
from ..mechanisms.binding import OneStepBinding
from ..mechanisms.enzyme import BasicCatalysis, MichaelisMenten
from ..mechanisms.global_mechanisms import Dilution, RNAdegradation_MM
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
    """Gene expression without any machinery (ribosomes, polymerases, etc.).

    Here transcription and Translation are lumped into one reaction:
    expression.

    """

    def __init__(self, name='', **kwargs):
        """Initializes an ExpressionExtract instance.

        :param name: name of the mixture
        :param kwargs: kwargs passed into the parent Class (Mixture)
        """
        # always call the superlcass Mixture.__init__(...)
        Mixture.__init__(self, name=name, **kwargs)

        # Create default Expression Mechanisms
        dummy_translation = EmptyMechanism(
            name='dummy_translation', mechanism_type='translation'
        )
        mech_expression = OneStepGeneExpression()
        mech_cat = BasicCatalysis()
        mech_bind = OneStepBinding()

        default_mechanisms = {
            mech_expression.mechanism_type: mech_expression,
            dummy_translation.mechanism_type: dummy_translation,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }

        self.add_mechanisms(default_mechanisms, overwrite=None)

    def compile_crn(self, **kwargs) -> ChemicalReactionNetwork:
        """Compile CRN, turning off transcription.

        Overwriting compile_crn to turn off transcription in all
        DNAassemblies.

        :return: compiled CRN instance

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
    """Transcription and translation in extract w/out any machinery.

    Transcriptoin and translation without ribosomes, polymerases,
    etc.  RNA is degraded via a global mechanism.

    """

    def __init__(self, name='', **kwargs):
        """Initializes a SimpleTxTlExtract instance.

        :param name: name of the mixture
        :param kwargs: kwargs passed into the parent Class (Mixture)
        """
        # Always call the superlcass Mixture.__init__(...)
        Mixture.__init__(self, name=name, **kwargs)

        # TxTl Mechanisms
        mech_tx = SimpleTranscription()
        mech_tl = SimpleTranslation()
        mech_cat = BasicCatalysis()
        mech_bind = OneStepBinding()

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
    """Transcription and translation with expression machinery.

    A Model for Transcription and Translation in Cell Extract with
    Ribosomes, Polymerases, and Endonucleases.

    This model does not include any energy.

    """

    def __init__(
        self, name='', rnap='RNAP', ribosome='Ribo', rnaase='RNAase', **kwargs
    ):
        """Initializes a TxTlExtract instance.

        :param name: name of the mixture
        :param rnap: name of the RNA polymerase, default: RNAP
        :param ribosome: name of the ribosome, default: Ribo
        :param rnaase: name of the Ribonuclease, default: RNAase
        :param kwargs: kwargs passed into the parent Class (Mixture)

        """
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
        mech_rna_deg = RNAdegradation_MM(nuclease=self.rnaase.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = OneStepBinding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)


class EnergyTxTlExtract(Mixture):
    """Transcription and translation in extract with machinery, energy.

    This model include energy carrier molcules in the form of NTPs,
    Amino Acids, and a Fuel Species (such as 3PGA) used for NTP
    regeneration. This model is equivalent to TxTl extract, but with
    limited fuel. Note that different amino acids and nucleotides are
    lumped together.

    Energy usage for transcription and translation is length dependent.

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
        """Initailize the TX-TL mixture.

        :param name: name of the mixture
        :param rnap: name of the RNA polymerase, default: RNAP
        :param ribosome: name of the ribosome, default: Ribo
        :param rnaase: name of the Ribonuclease, default: RNAase
        :param ntps: name of the nucleotide fuel source (eg ATP + GTP etc),
            default: NTP
        :param amino_acids: name of the amino acids species, default:
            amino_acids
        :param fuel: name of the fuel species that regenerates ATP
        :param kwargs: keywords passed into the parent Class (Mixture)

        """
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
        mech_rna_deg = RNAdegradation_MM(nuclease=self.rnaase.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = OneStepBinding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
        }
        self.add_mechanisms(default_mechanisms, overwrite=None)
