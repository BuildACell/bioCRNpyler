
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from warnings import resetwarnings
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .mechanism import EmptyMechanism
from .mechanisms_enzyme import BasicCatalysis, MichaelisMenten
from .mechanisms_binding import One_Step_Binding
from .mechanisms_txtl import Transcription_MM, Translation_MM, OneStepGeneExpression, SimpleTranscription, SimpleTranslation
from .global_mechanism import Dilution, Degredation_mRNA_MM
from .mixture import Mixture
from .chemical_reaction_network import ChemicalReactionNetwork
from .species import Species
from .dna_assembly import DNAassembly
       

class ExpressionExtract(Mixture):
    """A Model for Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.)

    Here transcription and Translation are lumped into one reaction: expression.
    """
    def __init__(self, name="", **kwargs):
        #always call the superlcass Mixture.__init__(...) 
        Mixture.__init__(self, name=name, **kwargs)

        #Create default Expression Mechanisms
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

        self.add_mechanisms(default_mechanisms)

    #Overwriting compile_crn to turn off transcription in all DNAassemblies
    def compile_crn(self) -> ChemicalReactionNetwork:

        for component in self.components:
            if isinstance(component, DNAassembly):
                #Only turn off transcription for an Assembly that makes a Protein. 
                #Some assemblies might only make RNA!
                if component.protein is not None:
                     #This will turn off transcription and set Promoter.transcript = False
                     #Mechanisms that recieve no transcript but a protein will use the protein instead.
                    component.update_transcript(False)

        #Call the superclass function
        return Mixture.compile_crn(self)


class SimpleTxTlExtract(Mixture):
    """
    A Model for Transcription and Translation in an extract any Machinery (eg Ribosomes, Polymerases, etc.)
    RNA is degraded via a global mechanism.
    """

    def __init__(self, name="", **kwargs):
        #Always call the superlcass Mixture.__init__(...) 
        Mixture.__init__(self, name=name, **kwargs)

        #TxTl Mechanisms
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
        self.add_mechanisms(default_mechanisms)

        #global mechanisms for dilution and rna degredation
        mech_rna_deg_global = Dilution(name = "rna_degredation", filter_dict = {"rna":True}, default_on = False)
        global_mechanisms = {"rna_degredation":mech_rna_deg_global}
        self.add_mechanisms(global_mechanisms)        


class TxTlExtract(Mixture):
    """A Model for Transcription and Translation in Cell Extract with Ribosomes, Polymerases, and Endonucleases.

    This model does not include any energy
    """
    def __init__(self, name="", rnap="RNAP", ribosome="Ribo", rnaase="RNAase", **kwargs):
        #Always call the superlcass Mixture.__init__(...) 
        Mixture.__init__(self, name=name, **kwargs)
        
        #create default Components to represent cellular machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        init = kwargs.get('init')
        if init:
            self.rnap.get_species().initial_concentration = init[repr(rnap)]
            self.rnaase.get_species().initial_concentration = init[repr(rnaase)]
            self.ribosome.get_species().initial_concentration = init[repr(ribosome)]

        default_components = [self.rnap, self.ribosome, self.rnaase]
        self.add_components(default_components)

        #Create default TxTl Mechanisms
        mech_tx = Transcription_MM(rnap = self.rnap.get_species())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_species())
        mech_rna_deg = Degredation_mRNA_MM(nuclease = self.rnaase.get_species()) 
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()


        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type: mech_bind
        }
        self.add_mechanisms(default_mechanisms)
