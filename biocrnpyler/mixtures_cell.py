from warnings import warn
from warnings import resetwarnings
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .mechanism import EmptyMechanism
from .mechanisms_txtl import Transcription_MM, Translation_MM, Degredation_mRNA_MM, OneStepGeneExpression, SimpleTranscription, SimpleTranslation
from .mixture import Mixture
from .chemical_reaction_network import Species, ChemicalReactionNetwork
from .global_mechanism import Dilution
from .dna_assembly import DNAassembly



#A Model for in-vivo Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.)
# Here transcription and Translation are lumped into one reaction: expression.
#A global mechanism is used to dilute all non-dna species
class ExpressionDilutionMixture(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        dummy_translation = EmptyMechanism(name = "dummy_translation", mechanism_type = "translation")
        mech_expression = OneStepGeneExpression()

        default_mechanisms = {
            "transcription": mech_expression,
            "translation": dummy_translation
        }

        dilution_mechanism= Dilution(name = "dilution", filter_dict = {"dna":False}, default_on = True)
        global_mechanisms = {"dilution":dilution_mechanism}

        default_components = []

        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, global_mechanisms = global_mechanisms, **kwargs)

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

#A Mixture with continious dilution for non-DNA species
#mRNA is also degraded via a seperate reaction to represent endonucleases
class SimpleTxTlDilutionMixture(Mixture):
    def __init__(self, name="", **keywords):
        
        simple_transcription = SimpleTranscription() #Transcription will not involve machinery
        simple_translation = SimpleTranslation()
        
        default_mechanisms = {
            "transcription": simple_transcription,
            "translation": simple_translation
        }
    
        #By Default Species are diluted S-->0 Unless:
        # They are of type 'dna'
        # They have the attribute 'machinery'
        dilution_mechanism = Dilution(filter_dict = {"dna":False}, default_on = True)
        dilution_mrna = Dilution(name = "rna_degredation", filter_dict = {"rna":True}, default_on = False)

        global_mechanisms = {"dilution":dilution_mechanism, "rna_degredation":dilution_mrna}
        
        #Always call the superclass __init__ with **keywords
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, global_mechanisms = global_mechanisms, **keywords)

#A Model for Transcription and Translation with Ribosomes, Polymerases, and Endonucleases labelled as Machinery. 
#Unlike TxTlExtract, has global dilution for non-DNA and non-Machinery
#This model does not include any energy

#TODO:
#Include some "internal" gene which provides background loading of all machinery
class TxTlDilutionMixture(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                 rnap = "RNAP", ribosome = "Ribo", rnaase = "RNAase", **kwargs):
        
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        self.rnap.add_attribute("machinery")
        self.ribosome.add_attribute("machinery")
        self.rnaase.add_attribute("machinery")

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

        dilution_mechanism = Dilution(filter_dict = {"dna":False, "machinery":False}, default_on = True)
        global_mechanisms = {"dilution":dilution_mechanism}

        background_parameters = {("transcription", "ku"):50, ("transcription", "kb"):500, ("transcription", "ktx"):0.1, 
              ("translation","ku"):5, ("translation","kb"):500, ("translation", "ktl"):.1,
              ("rna_degredation","ku"):50, ("rna_degredation","kb"):500, ("rna_degredation", "kdeg"):0.1}

        BackgroundProcesses = DNAassembly(name = "cellular_processes", promoter = "average_promoter", rbs = "average_rbs", parameters = background_parameters)

        default_components = [
            self.rnap, self.ribosome, self.rnaase, BackgroundProcesses
        ]

        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, global_mechanisms = global_mechanisms, **kwargs)