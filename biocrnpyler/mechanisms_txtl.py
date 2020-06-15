from .mechanism import *
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer
from .mechanisms_enzyme import *


class Transcription_MM(MichalisMentenCopyRXN):
    """Michalis Menten Transcription
        G + RNAP <--> G:RNAP --> G+RNAP+mRNA
    """

    def __init__(self, name="transcription_mm", rnap="RNAP", **keywords):
        if isinstance(rnap, Species):
            self.rnap = rnap
        elif isinstance(rnap, str):
            self.rnap = Species(name=rnap, material_type="protein")
        elif isinstance(rnap, Component) and rnap.get_species() != None:
            self.rnap = rnap.get_species()
        else:
            raise ValueError(
                "'rnap' parameter must be a string or a Component with defined "
                "get_species(), or a chemical_reaction_network.Species object")

        MichalisMentenCopyRXN.__init__(self=self, name=name, enzyme=self.rnap,
                                       mechanism_type="transcription")

    def update_species(self, dna, transcript=None, return_rnap=True,
                       **keywords):
        species = [dna]
        if return_rnap:
            species += [self.rnap]

        species += MichalisMentenCopyRXN.update_species(self, dna)
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")
        
        species += [transcript]

        return species

    def update_reactions(self, dna, component, part_id = None, complex=None, transcript=None,
                         **keywords):

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        rxns = []
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")
        rxns += MichalisMentenCopyRXN.update_reactions(self, dna, transcript,
                                                       complex=complex, kb=kb,
                                                       ku=ku, kcat=ktx)

        return rxns


class Translation_MM(MichalisMentenCopyRXN):
    """ Michalis Menten Translation
        mRNA + Rib <--> mRNA:Rib --> mRNA + Rib + Protein
    """

    def __init__(self, name="translation_mm", ribosome="Ribo", **keywords):
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        elif isinstance(ribosome, str):
            self.ribosome = Species(name=ribosome, material_type="ribosome")
        elif isinstance(ribosome, Component) and ribosome.get_species() != None:
            self.ribosome = ribosome.get_species()
        else:
            raise ValueError(
                "'ribosome' parameter must be a string, a Component with defined "
                "get_species, or a chemical_reaction_network.species")
        MichalisMentenCopyRXN.__init__(self=self, name=name,
                                       enzyme=self.ribosome,
                                       mechanism_type="translation")

    def update_species(self, transcript, protein=None,
                       return_ribosome=True, **keywords):
        species = [transcript]
        if return_ribosome:
            species += [self.ribosome]
        if protein is None:
            protein = Species(transcript.name, material_type="protein")
        species += [protein]

        species += MichalisMentenCopyRXN.update_species(self, transcript)

        return species

    def update_reactions(self, transcript, component, part_id = None, complex=None,
                         protein=None, **keywords):
        rxns = []

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktl = component.get_parameter("ktl", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        if protein is None:
            protein = Species(transcript.name, material_type="protein")
        rxns += MichalisMentenCopyRXN.update_reactions(self, transcript,
                                                       protein, complex=complex,
                                                       kb=kb, ku=ku,
                                                       kcat=ktl)
        return rxns


class Degredation_mRNA_MM(MichalisMentenRXN):
    """Michalis Menten mRNA Degredation by Endonucleases
       mRNA + Endo <--> mRNA:Endo --> Endo
    """
    def __init__(self, name="rna_degredation_mm", nuclease="RNAase",
                 **keywords):
        if isinstance(nuclease, Species):
            self.nuclease = nuclease
        elif isinstance(nuclease, str):
            self.nuclease = Species(name=nuclease, material_type="protein")
        else:
            raise ValueError("'nuclease' parameter requires a "
                             "chemical_reaction_network.species or a string")
        MichalisMentenRXN.__init__(self=self, name=name, enzyme=self.nuclease,
                                   mechanism_type="rna_degredation")

    def update_species(self, rna, return_nuclease=True, **keywords):
        species = [rna]
        if return_nuclease:
            species += [self.nuclease]
        species += MichalisMentenRXN.update_species(self, rna)
        return species

    def update_reactions(self, rna, component, part_id = None, complex=None, **keywords):

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        kdeg = component.get_parameter("kdeg", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        rxns = []
        rxns += MichalisMentenRXN.update_reactions(self, rna, Prod=None, complex=complex, kb=kb, ku=ku, kcat=kdeg)
        return rxns

class SimpleTranscription(Mechanism):
    def __init__(self, name = "simple_transcription", mechanism_type = "transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, transcript = None, **keywords):
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")

        return [dna, transcript]

    def update_reactions(self, dna, component = None, ktx = None, part_id = None, transcript = None, **keywords):

        if ktx == None and Component != None:
            ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        elif component == None and ktx == None:
            raise ValueError("Must pass in component or a value for ktx")

        if transcript is None:
            transcript = Species(dna.name, material_type="rna")

        rxns = [Reaction(inputs = [dna], outputs = [dna, transcript], k = ktx)]
        return rxns

class SimpleTranslation(Mechanism):
    def __init__(self, name = "simple_translation", mechanism_type = "translation"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, transcript, protein = None,  **keywords):
        if protein is None:
            protein = Species(transcript.name, material_type="protein")

        return [transcript, protein]

    def update_reactions(self, transcript, component = None, ktl = None, part_id = None, protein = None, **keywords):

        if ktl == None and Component != None:
            ktl = component.get_parameter("ktl", part_id = part_id, mechanism = self)
        elif component == None and ktl == None:
            raise ValueError("Must pass in component or a value for ktl")

        if protein is None:
            protein = Species(transcript.name, material_type="protein")

        rxns = [Reaction(inputs = [transcript], outputs = [transcript, protein], k = ktl)]
        return rxns

class OneStepGeneExpression(Mechanism):
    def __init__(self, name="gene_expression",
                 mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, protein=None, transcript=None, **keywords):
        species = [dna]
        if protein == None:
            protein = Species(dna.name, material_type="protein")

        species += [protein]
        return species

    def update_reactions(self, dna, component = None, kexpress = None,
                         protein=None, transcript = None, part_id = None, **keywords):

        if kexpress == None and Component != None:
            kexpress = component.get_parameter("kexpress", part_id = part_id, mechanism = self)
        elif component == None and kexpress == None:
            raise ValueError("Must pass in component or a value for kexpress")

        if protein is None:
            protein = Species(dna.name, material_type="protein")
        rxns = [Reaction(inputs=[dna], outputs=[dna, protein], k = kexpress)]
        return rxns


class PositiveHillTranscription(Mechanism):
    #Set the name and mechanism_type
    def __init__(self, name="positivehill_transcription", mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)
    
    #Overwrite update_species
    def update_species(self, dna, regulator, transcript = None, **keywords):
        
        if transcript is None: #Species names can be automatically created
            transcript = Species(dna.name, material_type = "rna")
            
        return [dna, transcript, regulator] #it is best to return all species that will be involved in the reactions

    
    #Overwrite update_reactions
    #This always requires the inputs component and part_id to find the relevant parameters
    def update_reactions(self, dna, regulator, component, part_id, transcript = None, **keywords):
        
        if transcript is None: #Species names should be automatically created the same here as above
            transcript = Species(dna.name, material_type = "rna")
        
        ktx = component.get_parameter("k", part_id = part_id, mechanism = self)
        n = component.get_parameter("n", part_id = part_id, mechanism = self)
        K = component.get_parameter("K", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)
        
        params = {"k":ktx, "n":n, "K":K, "s1":regulator, "d":dna}
        
        reaction = Reaction(inputs = [dna], outputs = [dna, transcript], 
                            propensity_type = "proportionalhillpositive", propensity_params = params)
        
        reaction_leak = Reaction(inputs = [dna], outputs = [dna, transcript], k = kleak)
        
        #In this case, we just return one reaction
        return [reaction, reaction_leak]

class NegativeHillTranscription(Mechanism):
    #Set the name and mechanism_type
    def __init__(self, name="negativehill_transcription", mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)
    
    #Overwrite update_species
    def update_species(self, dna, regulator, transcript = None, **keywords):
        
        if transcript is None: #Species names can be automatically created
            transcript = Species(dna.name, material_type = "rna")
            
        return [dna, transcript, regulator] #it is best to return all species that will be involved in the reactions
    
    #Overwrite update_reactions
    #This always requires the inputs component and part_id to find the relevant parameters
    def update_reactions(self, dna, regulator, component, part_id, transcript = None, **keywords):
        
        if transcript is None: #Species names should be automatically created the same here as above
            transcript = Species(dna.name, material_type = "rna")
        
        ktx = component.get_parameter("k", part_id = part_id, mechanism = self)
        n = component.get_parameter("n", part_id = part_id, mechanism = self)
        K = component.get_parameter("K", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)
        
        params = {"k":ktx, "n":n, "K":K, "s1":regulator, "d":dna}
        
        reaction = Reaction(inputs = [dna], outputs = [dna, transcript], 
                            propensity_type = "proportionalhillnegative", propensity_params = params)
        
        reaction_leak = Reaction(inputs = [dna], outputs = [dna, transcript], k = kleak)
        
        #In this case, we just return one reaction
        return [reaction, reaction_leak]
