# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from warnings import resetwarnings
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .mechanism import EmptyMechanism
from .mechanisms_enzyme import BasicCatalysis, MichalisMenten
from .mechanisms_binding import One_Step_Binding
from .mechanisms_txtl import Transcription_MM, Translation_MM, Degredation_mRNA_MM, OneStepGeneExpression, SimpleTranscription, SimpleTranslation
from .global_mechanism import Dilution
from .mixture import Mixture
from .chemical_reaction_network import Species, ChemicalReactionNetwork
from .dna_assembly import DNAassembly
       
#A Model for Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.)
# Here transcription and Translation are lumped into one reaction: expression.
class ExpressionExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        dummy_translation = EmptyMechanism(name = "dummy_translation", mechanism_type = "translation")
        mech_expression = OneStepGeneExpression()
        mech_cat = BasicCatalysis()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_expression.mechanism_type: mech_expression,
            dummy_translation.mechanism_type: dummy_translation,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind
        }

        default_components = []
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

    #Overwriting compile_crn to replace transcripts with proteins for all DNA_assemblies
    def compile_crn(self) -> ChemicalReactionNetwork:
        """ Creates a chemical reaction network from the species and reactions associated with a mixture object
        :return: ChemicalReactionNetwork
        """
        resetwarnings()#Reset warnings - better to toggle them off manually.
        species = self.update_species()
        reactions = self.update_reactions()

        for comp in self.components:
            if isinstance(comp, DNAassembly):
                if comp.transcript is not None and comp.protein is not None:
                    for i, s in enumerate(species):
                        species[i] = s.replace_species(comp.transcript, comp.protein)
                    for i, r in enumerate(reactions):
                        reactions[i] = r.replace_species(comp.transcript, comp.protein)

        self.crn_species = list(set(species))
        self.crn_reactions = reactions
        #global mechanisms are applied last and only to all the species 
        global_mech_species, global_mech_reactions = self.apply_global_mechanisms()

        species += global_mech_species
        reactions += global_mech_reactions

        species = self.set_initial_condition(species)
        species.sort(key = lambda s:repr(s))
        reactions.sort(key = lambda r:repr(r))
        CRN = ChemicalReactionNetwork(species, reactions)
        return CRN

#A Model for Transcription and Translation in an extract any Machinery (eg Ribosomes, Polymerases, etc.)
#RNA is degraded via a global mechanism
class SimpleTxTlExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        mech_tx = SimpleTranscription()
        mech_tl = SimpleTranslation()
        mech_cat = BasicCatalysis()
        mech_bind = One_Step_Binding()

        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind
        }

        mech_rna_deg_global = Dilution(name = "rna_degredation", filter_dict = {"rna":True}, default_on = False)
        global_mechanisms = {"rna_degredation":mech_rna_deg_global}

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
        mech_cat = MichalisMenten()
        mech_bind = One_Step_Binding()


        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind
        }

        default_components = [self.rnap, self.ribosome, self.rnaase]
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)
