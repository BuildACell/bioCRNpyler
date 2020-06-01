from warnings import warn
from .basic_components import DNA, RNA, Protein, ChemicalComplex
from .mechanism import EmptyMechanism
from .txtl_mechanisms import Transcription_MM, Translation_MM, Degredation_mRNA_MM, OneStepGeneExpression, SimpleTranscription, SimpleTranslation
from .mixture import Mixture
from .chemical_reaction_network import Species
from .global_mechanism import Dilution


#A Mixture with continious dilution for non-DNA species
#mRNA is also degraded via a seperate reaction to represent endonucleases
class DilutionMixture(Mixture):
    def __init__(self, name="", **keywords):
        
        simple_transcription = SimpleTranscription() #Transcription will not involve machinery
        simple_translation = SimpleTranslation()
        
        default_mechanisms = {
            "transcription": simple_transcription, #This will be overwritten by the NegativeHillPromotor
            "translation": simple_translation
        }
    
        #By Default Species are diluted S-->0 Unless:
        # They are of type 'dna'
        # They have the attribute 'machinery'
        dilution_mechanism = Dilution(filter_dict = {"dna":False}, default_on = True)
        dilution_mrna = Dilution(name = "rna_degredation", filter_dict = {"rna":True}, default_on = False)

        #Add this mechanism to a dictionary which is passed into the Mixture txtl.TxTlExtract
        global_mechanisms = {"dilution":dilution_mechanism, "rna_degredation":dilution_mrna}
        
        #Always call the superclass __init__ with **keywords
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, global_mechanisms = global_mechanisms, **keywords)