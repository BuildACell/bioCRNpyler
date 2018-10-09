import component
import chemical_reaction_network as crn
from warnings import warn as pywarn
def warn(txt):
    pywarn(txt)
import mechanism



class Promoter(component.Component):
    def __init__(self, name, assembly = None,
                 transcript = None, length = 0,
                 mechanisms = {}, parameters = {}, **keywords):
        self.assembly = assembly
        self.length = length
        if transcript == None and assembly == None:
            self.transcript = None
        elif transcript == None:
            self.transcript = crn.specie(assembly.name, type = "rna")
        elif isinstance(transcript, str):
            self.transcript = crn.specie(transcript, type = "rna")
        elif isinstance(transcript, crn.specie):
            self.transcript = transcript
        else:
            raise ValueError("Improper value of transcriped passed into promoter. Transcript should be a string or chemical_reaction_network.specie")

        component.Component.__init__(self, name = name, mechanisms = mechanisms, parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        species = [self.transcript]
        species += mech_tx.update_species(self.assembly.dna)
        return species

    def update_reactions(self):
        mech_tx = self.mechanisms["transcription"]
        reactions = []
        ktx = self.get_parameter("ktx")
        kb = self.get_parameter("kb")
        ku = self.get_parameter("ku")

        reactions += mech_tx.update_reactions(self.assembly.dna, ktx = ktx, ku = ku, kb = kb,
                                              complex = None, transcript = self.transcript)
        return reactions


class RegulatedPromoter(Promoter):
    def __init__(self, name, regulators, leak = False, assembly=None, transcript=None, length=0,
                 mechanisms={}, parameters={}, **keywords):

        if not isinstance(regulators, list):
            regulators = [regulators]

        self.regulators = []
        for regulator in regulators:
            if isinstance(regulator, crn.specie):
                self.regulators += [regulator]
            elif isinstance(regulator, str):
                self.regulators += [crn.specie(name = regulator, type = "protein")]
            else:
                raise ValueError("Invalid parameter 'regulators' passed to RegulatedPromoter: valid types are string, chemical_reaction_network.specie and lists of these strings or species")

        self.leak = leak

        self.default_mechanisms = {"cooperative_binding":mechanism.One_Step_Cooperative_Binding()}

        Promoter.__init__(self, name = name, assembly = assembly, transcript = transcript, length = length, mechanisms = mechanisms, parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['cooperative_binding']

        species = []
        self.complexes = []
        if self.leak != False:
            species += mech_tx.update_species(dna=self.assembly.dna)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            coop = self.get_parameter("cooperativity", regulator.name, mechanism = mech_b)
            species_b = mech_b.update_species(binder = regulator, bindee = self.assembly.dna, cooperativity = coop)
            species += species_b
            complex = species_b[0]
            self.complexes += [complex]
            species += mech_tx.update_species(dna = complex)
        return species


    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['cooperative_binding']

        if self.leak != False:
            ktx = self.get_parameter("ktx_leak", mechanism = mech_tx)
            ku = self.get_parameter("ku_leak", mechanism = mech_tx, part_id=self.name)
            kb = self.get_parameter("kb_leak", mechanism = mech_tx, part_id=self.name)
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, ktx = ktx, ku = ku, kb = kb)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex = self.complexes[i]
            ktx = self.get_parameter("ktx", mechanism = mech_tx, part_id = regulator.name)

            ku_tx = self.get_parameter("ku", mechanism = mech_tx, part_id = regulator.name)
            kb_tx = self.get_parameter("kb", mechanism = mech_tx, part_id = regulator.name)

            ku_c = self.get_parameter("ku", mechanism=mech_b, part_id=regulator.name)
            kb_c = self.get_parameter("kb", mechanism=mech_b, part_id=regulator.name)

            coop = self.get_parameter("cooperativity", part_id = regulator.name, mechanism = mech_tx)

            reactions += mech_b.update_reactions(binder = regulator, bindee = self.assembly.dna, ku = ku_c, kb = kb_c, cooperativity = coop)
            reactions += mech_tx.update_reactions(dna = complex, kb = kb_tx, ku = ku_tx, ktx = ktx, transcript = self.transcript)

        return reactions


class RBS(component.Component):
    def __init__(self, name, assembly,
                 transcript=None, protein = None, length=0,
                 mechanisms={}, parameters={}, **keywords):
        self.assembly = assembly
        self.length = length

        component.Component.__init__(self, name=name, mechanisms=mechanisms, parameters=parameters, **keywords)

        if transcript == None and assembly == None:
            self.transcript = None
        elif transcript == None:
            self.transcript = crn.specie(assembly.name, type="rna")
        elif isinstance(transcript, str):
            self.transcript = crn.specie(transcript, type="rna")
        elif isinstance(transcript, crn.specie):
            self.transcript = transcript
        else:
            raise ValueError(
                "Improper value of transcript passed into promoter: transcript should be a string or chemical_reaction_network.specie")

        if protein == None:
            self.protein = crn.specie(assembly.name, type="protein")
        elif isinstance(protein, str):
            self.protein = crn.specie(protein, type="protein")
        elif isinstance(protein, crn.specie):
            self.protein = protein
        else:
            raise ValueError(
                "Improper value of protein passed into RBS: protein should be a string or chemical_reaction_network.specie")

        component.Component.__init__(self, name = name, mechanisms = mechanisms, parameters = parameters, **keywords)

    def update_species(self):
        mech_tl = self.mechanisms['translation']
        species = [self.transcript, self.protein]
        species += mech_tl.update_species(self.transcript)
        return species

    def update_reactions(self):
        mech_tl = self.mechanisms['translation']
        reactions = []

        ktl = self.get_parameter("ktl")
        kb = self.get_parameter("kb")
        ku = self.get_parameter("ku")
        reactions += mech_tl.update_reactions(transcript = self.transcript, protein = self.protein,
                                              ku = ku, kb = kb, ktl = ktl)
        return reactions

class DNAassembly(component.DNA):
    def __init__(self, name, dna = None,
                 promoter = None, transcript = None,
                 rbs=None, protein = None,
                 mechanisms = {}, parameters = {}, **keywords):

        self.name = name
        self.rbs, self.promoter = None, None #Filler until updates

        length = 0
        warn("length not yet implemented for DNA assemblies")
        component.DNA.__init__(self, name, length=length, mechanisms=mechanisms, parameters=parameters, **keywords)

        self.update_dna(dna)
        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_promoter(promoter, transcript = self.transcript)
        self.update_rbs(rbs, transcript=self.transcript, protein=self.protein)



    def update_dna(self, dna):
        if isinstance(dna, crn.specie):
            self.dna = dna
        elif isinstance(dna, str):
            self.dna = crn.specie(dna, type = "dna")
        elif dna == None:
            self.dna = crn.specie(self.name, type = "dna")
        else:
            raise ValueError("Invalid value of 'dna' passed to "+str(self)+": dna must be None, a chemical_reaction_network.specie, or a string.")

    def update_transcript(self, transcript):
        if isinstance(transcript, crn.specie):
            self.transcript = transcript
        elif isinstance(transcript, str):
            self.transcript = crn.specie(transcript, type = "rna")
        elif transcript == None:
            self.transcript = crn.specie(self.name, type = "rna")
        else:
            raise ValueError("Invalid value of transcript passed to "+str(self)+": transcript must be None, a chemical_reaction_network.specie, or a string.")

        if self.promoter != None:
            self.promoter.transcript = self.transcript
        if self.rbs != None:
            self.rbs.transcript = self.transcript

    def update_protein(self, protein):
        if isinstance(protein, crn.specie):
            self.protein = protein
        elif isinstance(protein, str):
            self.protein = crn.specie(protein, type="protein")
        elif protein == None:
            self.protein = crn.specie(self.name, type="protein")
        else:
            raise ValueError("Invalid value of 'protein' passed to DNA Assembly " + self.name + ". 'protein' must be None, a chemical_reaction_network.specie, or a string.")

        if self.rbs != None:
            self.rbs.transcript = self.protein

    def update_promoter(self, promoter, transcript = None):
        if transcript != None:
            self.update_transcript(transcript)

        if isinstance(promoter, str):
            self.promoter = Promoter(assembly= self, name = promoter, transcript = self.transcript, parameters = self.parameters)
        elif isinstance(promoter, Promoter):
            self.promoter = promoter
            self.promoter.assembly = self
            self.promoter.transcript = self.transcript
        elif promoter != None:
            raise ValueError("Improper promoter type recieved by DNAassembly. Expected string or promoter object. Recieved "+repr(promoter))



    def update_rbs(self, rbs, transcript = None, protein = None):
        if protein != None:
            self.update_protein(protein)

        if transcript != None:
            self.update_transcript(transcript)

        if isinstance(rbs, str):
            self.rbs = RBS(assembly=self, name = rbs, protein = self.protein, transcript = self.transcript, parameters = self.parameters)
        elif isinstance(rbs, RBS):
            self.rbs = rbs
            self.rbs.assembly = self
            self.rbs.transcript = self.transcript
            self.rbs.protein = self.protein
        elif rbs != None:
            raise ValueError("Improper rbs type recieved by DNAassembly. Expected string or RBS object. Recieved "+repr(rbs))



    def update_species(self):
        species = []
        species.append(self.dna)
        if self.promoter != None and self.rbs != None:
            #species += [self.transcript, self.protein]
            species += self.promoter.update_species()
            species += self.rbs.update_species()

        elif self.promoter != None and self.rbs == None:
            #species += [self.transcript]
            species+=self.promoter.update_species()

        if "rna_degredation" in self.mechanisms and self.promoter != None:
            deg_mech = self.mechanisms["rna_degredation"]
            species += deg_mech.update_species(self.transcript)

        return list(set(species))

    def update_reactions(self):
        reactions = []
        if self.promoter != None:
            reactions+=self.promoter.update_reactions()

        if self.rbs != None:
            reactions+=self.rbs.update_reactions()

        if "rna_degredation" in self.mechanisms and self.promoter != None:
            deg_mech = self.mechanisms["rna_degredation"]
            ku = self.get_parameter("ku", mechanism =deg_mech)
            kb = self.get_parameter("kb", mechanism=deg_mech)
            kdeg = self.get_parameter("kdeg", mechanism=deg_mech)

            reactions += deg_mech.update_reactions(self.transcript, ku = ku, kb = kb, kdeg = kdeg)

        return reactions

    def update_parameters(self, mixture_parameters = {}, parameters={}):
        component.DNA.update_parameters(self = self, mixture_parameters = mixture_parameters, parameters=parameters)

        if self.promoter != None:
            self.promoter.update_parameters(mixture_parameters=mixture_parameters, parameters=parameters, overwrite_custom_parameters=False)
        if self.rbs != None:
            self.rbs.update_parameters(mixture_parameters=mixture_parameters, parameters=parameters, overwrite_custom_parameters=False)

    def update_mechanisms(self, mixture_mechanisms={}, mechanisms={}):
        component.DNA.update_mechanisms(self = self, mixture_mechanisms = mixture_mechanisms, mechanisms=mechanisms)

        if self.promoter != None and "transcription" in self.mechanisms:
            mech_tx = self.mechanisms["transcription"]
            self.promoter.update_mechanisms(mechanisms={"transcription":mech_tx}, overwrite_custom_mechanisms=False)
        if self.rbs != None and "translation" in self.mechanisms:
            mech_tl = self.mechanisms["translation"]
            self.rbs.update_mechanisms(mechanisms={"translation":mech_tl}, overwrite_custom_mechanisms=False)

    def __str__(self):
        return type(self).__name__+": "+ self.name

    def __repr__(self):
        txt = str(self)
        if self.promoter != None:
            txt+="\n\t"+repr(self.promoter)
            txt+="\n\ttranscript = "+repr(self.transcript)
        if self.rbs != None:
            txt+="\n\t"+repr(self.rbs)
            txt+="\n\tprotein = "+repr(self.protein)
        return txt


