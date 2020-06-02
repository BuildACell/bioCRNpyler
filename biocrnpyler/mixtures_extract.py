# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .mechanism import EmptyMechanism
from .mechanisms_txtl import Transcription_MM, Translation_MM, Degredation_mRNA_MM, OneStepGeneExpression, SimpleTranscription, SimpleTranslation
from .mixture import Mixture
from .chemical_reaction_network import Species
       
#A Model for Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.)
# Here transcription and Translation are lumped into one reaction: expression.
class ExpressionExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        dummy_translation = EmptyMechanism(name = "dummy_translation", mechanism_type = "translation")
        mech_expression = OneStepGeneExpression()

        default_mechanisms = {
            "transcription": mech_expression,
            "translation": dummy_translation
        }

        default_components = []
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

#A Model for Transcription and Translation in an extract any Machinery (eg Ribosomes, Polymerases, etc.)
#RNA is degraded via a global mechanism
class SimpleTxTlExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        mech_tx = SimpleTranscription(name = "simple_transcription", mechanism_type = "transcription")
        mech_tl = SimpleTranscription(name = "simple_translation", mechanism_type = "translation")

        default_mechanisms = {
            "transcription": mech_tx,
            "translation": mech_tl
        }

        mech_rna_deg_global = Dilution(name = "rna_degredation", filter_dict = {"rna":True}, default_on = False)
        global_mechanisms = {"rna_degredation":dilution_mrna}

        default_components = []
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, global_mechanisms= global_mechanisms, **kwargs)

#A Model for Transcription and Translation in Cell Extract with Ribosomes, Polymerases, and Endonucleases.
#This model does not include any energy
class TxTlExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                 rnap = "RNAP", ribosome = "Ribo", rnaase = "RNAase", **kwargs):
        
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        init = kwargs.get('init')
        if init:
            self.rnap.get_species().initial_concentration = init[rep(rnap)]
            self.rnaase.get_species().initial_concentration = init[repr(rnaase)]
            self.ribosome.get_species().initial_concentration = init[repr(ribosome)]

        mech_tx = Transcription_MM(rnap = self.rnap.get_species())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_species())
        mech_rna_deg = Degredation_mRNA_MM(nuclease = self.rnaase.get_species()) 


        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.rnaase]
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)
