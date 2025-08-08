# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.chemical_reaction_network import ChemicalReactionNetwork
from ..core.mechanism import EmptyMechanism
from ..core.mixture import Mixture

from ..components.basic import Protein
from ..components.dna.assembly import DNAassembly
from ..mechanisms.global_mechanisms import Degredation_mRNA_MM, Dilution
from ..mechanisms.binding import One_Step_Binding
from ..mechanisms.enzyme import BasicCatalysis, MichaelisMenten
from ..mechanisms.txtl import (
    OneStepGeneExpression,
    SimpleTranscription,
    SimpleTranslation,
    Transcription_MM,
    Translation_MM,
)


class ExpressionDilutionMixture(Mixture):
    """A Model for in-vivo Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.).


    Here transcription and Translation are lumped into one reaction: expression.
    A global mechanism is used to dilute all non-dna species
    """

    def __init__(self, name="", **kwargs):
        """Initializes an ExpressionDilutionMixture instance.

        :param name: name of the mixture
        :param kwargs: keywords passed into the parent Class (Mixture)
        """
        Mixture.__init__(self, name=name, **kwargs)

        # Create default mechanisms for Gene Expression
        dummy_translation = EmptyMechanism(
            name="dummy_translation", mechanism_type="translation"
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
        self.add_mechanisms(default_mechanisms)

        # Create global mechanism for dilution
        dilution_mechanism = Dilution(
            name="dilution", filter_dict={"dna": False}, default_on=True
        )
        global_mechanisms = {"dilution": dilution_mechanism}
        self.add_mechanisms(global_mechanisms)

    def compile_crn(self, **keywords) -> ChemicalReactionNetwork:
        """Overwriting compile_crn to replace transcripts with proteins for all DNA_assemblies.

        Overwriting compile_crn to turn off transcription in all DNAassemblies

        :return: compiled CRN instance
        """
        for component in self.components:
            if isinstance(component, DNAassembly):
                # Only turn off transcription for an Assembly that makes a Protein.
                # Some assemblies might only make RNA!
                if component.protein is not None:
                    # This will turn off transcription and set Promoter.transcript = False
                    # Mechanisms that recieve no transcript but a protein will use the protein instead.
                    component.update_transcript(False)

        # Call the superclass function
        return Mixture.compile_crn(self, **keywords)


class SimpleTxTlDilutionMixture(Mixture):
    """Mixture with continuous dilution for non-DNA species.

    Transcription and Translation are both modeled as catalytic with no cellular machinery.
    mRNA is also degraded via a separate reaction to represent endonucleases
    """

    def __init__(self, name="", **keywords):
        """Initializes a SimpleTxTlDilutionMixture instance.

        :param name: name of the mixture
        :param kwargs: keywords passed into the parent Class (Mixture)
        """
        # Always call the superclass __init__ with **keywords
        Mixture.__init__(self, name=name, **keywords)

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
        self.add_mechanisms(default_mechanisms)

        # Global Dilution Mechanisms
        # By Default Species are diluted S-->0 Unless:
        # They are of type 'dna'
        # They have the attribute 'machinery'
        dilution_mechanism = Dilution(filter_dict={"dna": False}, default_on=True)
        deg_mrna = Dilution(
            name="rna_degredation", filter_dict={"rna": True}, default_on=False
        )

        global_mechanisms = {
            "dilution": dilution_mechanism,
            "rna_degredation": deg_mrna,
        }
        self.add_mechanisms(global_mechanisms)


class TxTlDilutionMixture(Mixture):
    """A Model for Transcription and Translation with Ribosomes, Polymerases, and Endonucleases labelled as Machinery.

    This model includes a background load "cellular processes" which represents innate loading effects in the cell.
    Effects of loading on cell growth are not modelled.
    Unlike TxTlExtract, has global dilution for non-DNA and non-Machinery
    This model does not include any energy
    """

    def __init__(
        self, name="", rnap="RNAP", ribosome="Ribo", rnaase="RNAase", **kwargs
    ):
        """Initializes a TxTlDilutionMixture instance.

        :param name: name of the mixture
        :param rnap: name of the RNA polymerase, default: RNAP
        :param ribosome: name of the ribosome, default: Ribo
        :param rnaase: name of the Ribonuclease, default: RNAase
        :param kwargs: keywords passed into the parent Class (Mixture)
        """
        Mixture.__init__(self, name=name, **kwargs)

        # Create Components for TxTl machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        self.rnap.add_attribute("machinery")
        self.ribosome.add_attribute("machinery")
        self.rnaase.add_attribute("machinery")

        # DNAassmbly represents background processes / loading in a cell
        background_parameters = {
            ("transcription", None, "ku"): 50,
            ("transcription", None, "kb"): 500,
            ("transcription", None, "ktx"): 0.1,
            ("translation", None, "ku"): 5,
            ("translation", None, "kb"): 500,
            ("translation", None, "ktl"): 0.1,
            ("rna_degredation", None, "ku"): 50,
            ("rna_degredation", None, "kb"): 500,
            ("rna_degredation", None, "kdeg"): 0.1,
        }
        BackgroundProcesses = DNAassembly(
            name="cellular_processes",
            promoter="average_promoter",
            rbs="average_rbs",
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
            filter_dict={"dna": False, "machinery": False}, default_on=True
        )
        mech_rna_deg = Degredation_mRNA_MM(nuclease=self.rnaase.get_species())

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            "dilution": dilution_mechanism,
        }

        self.add_mechanisms(default_mechanisms)
