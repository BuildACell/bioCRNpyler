from .component import Component
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .chemical_reaction_network import ComplexSpecies, Species
from .mechanisms_binding import *
from .mechanisms_txtl import *
from warnings import warn as pywarn
import itertools as it
import numpy as np
from .dna_part import DNA_part

class Promoter(DNA_part):
    def __init__(self, name, assembly=None,
                 transcript=None, length=0,
                 mechanisms={}, parameters={}, **keywords):
        self.assembly = assembly
        self.length = length
        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type="rna")
        else:
            self.transcript = self.set_species(transcript, material_type = 'rna')

        DNA_part.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        species = []
        species += mech_tx.update_species(dna = self.assembly.dna, \
            transcript = self.transcript,
            component = self, part_id = self.name)
        return species

    def update_reactions(self):
        mech_tx = self.mechanisms["transcription"]
        reactions = []


        reactions += mech_tx.update_reactions(dna = self.assembly.dna, \
                        component = self, part_id = self.name, complex = None,
                        transcript = self.transcript)
        return reactions

class RegulatedPromoter(Promoter):
    def __init__(self, name, regulators, leak = True, assembly = None,
                 transcript = None, length = 0, mechanisms = {},
                 parameters = {}, **keywords):

        if not isinstance(regulators, list):
            regulators = [regulators]

        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]

        self.leak = leak

        self.default_mechanisms = {"binding": One_Step_Cooperative_Binding()}
        self.complexes = []
        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters,
                          **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']
        species = []
        self.complexes = []
        if self.leak is not False:
            species += mech_tx.update_species(dna = self.assembly.dna, component = self)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]

            species_b = mech_b.update_species(regulator, self.assembly.dna, part_id = self.name+"_"+regulator.name, component = self)
            species += species_b

            #Find complexes containing DNA and the regulator
            for s in species_b:
                if isinstance(s, ComplexSpecies) and self.assembly.dna in s.species and regulator in s.species:
                    self.complexes += [s]

                    species += mech_tx.update_species(dna = s, transcript = self.transcript, part_id = self.name+"_"+regulator.name, component = self)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak != False:
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, component = self, part_id = self.name, \
                                                            transcript = self.transcript)
        if(len(self.complexes)<len(self.regulators)):
            self.update_species()
        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]

            reactions += mech_b.update_reactions(regulator, self.assembly.dna, component = self, \
                                                                    part_id = self.name+"_"+regulator.name)
            reactions += mech_tx.update_reactions(dna = complex_, component = self, part_id = self.name+"_"+regulator.name, \
                                            transcript = self.transcript)

        return reactions

#A class for a promoter which can be activated by a single species, modelled as a positive hill function
class ActivatablePromotor(Promoter):
    def __init__(self, name, transcript, activator, **keywords):
        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        self.activator = self.set_species(activator) 
        
        #Mechanisms are inherited from the Mixture unless set specifically in self.default_mechanisms.
        custom_mechanisms = {"transcription": PositiveHillTranscription()}
        
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms,**keywords)

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.mechanisms["transcription"]
        
        species = [] #A list of species must be returned
        species += mech_tx.update_species(dna = self.assembly.dna, transcript = self.transcript, regulator = self.activator, part_id = self.name+"_"+activator.name, component = self)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = [] #a list of reactions must be returned
        reactions += mech_tx.update_reactions(dna = self.assembly.dna, transcript = self.transcript, regulator = self.repressor, 
                                             component = self, part_id = self.name+"_"+self.activator.name, **keywords)
        return reactions

#A class for a promoter which can be repressed by a single species, modelled as a negative hill function
class RepressablePromotor(Promoter):
    def __init__(self, name, transcript, repressor, **keywords):
        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        self.repressor = self.set_species(repressor) 
        
        #Mechanisms are inherited from the Mixture unless set specifically in self.default_mechanisms.
        custom_mechanisms = {"transcription": NegativeHillTranscription()}
        
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms,**keywords)

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.mechanisms["transcription"]
        
        species = [] #A list of species must be returned
        species += mech_tx.update_species(dna = self.assembly.dna, transcript = self.transcript, regulator = self.repressor, component = self, part_id = self.name+"_"+self.repressor.name)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = [] #a list of reactions must be returned
        reactions += mech_tx.update_reactions(dna = self.assembly.dna, transcript = self.transcript, regulator = self.repressor, 
                                             component = self, part_id = self.name+"_"+self.repressor.name, **keywords)
        return reactions

class CombinatorialPromoter(Promoter):
    def __init__(self, name, regulators, leak = True, assembly = None,
                 transcript = None, length = 0, mechanisms = {},
                 parameters = {},tx_capable_list = None,cooperativity = None, **keywords):
        """
        A combinatorial promoter is something where binding multiple regulators result in
        qualitatively different transcription behaviour. For example, maybe it's an AND
        gate promoter where it only transcribes if two regulators are bound, but not if 
        either one is bound.
        
        =============
        inputs
        =============
        name: the name of the promoter
        regulators: a list of strings or species indicating all the possible regualtors that can bind
        
        leak: if true, then a promoter with nothing bound will transcribe
        
        assembly: a DNA_assembly object that contains this promoter
        
        transcript: the transcript that this promoter makes
        
        length: the length in nt? I don't think this is used for anything at the moment
        
        mechanisms: additional mechanisms. formatted with {"mechanism_type":mechanismObject(),...}
        
        parameters: promoter-specific parameters. Formatted as {("identifier1","identifier2"):value,...}
        
        tx_capable_list: list of which combination of regulators bound will lead to transcription. 
                        formatted as [["regulator1","regulator2"],["regulator1"],...] regulators
                        can be strings or Species
        cooperativity: a dictionary of cooperativity values. For example, {"regulator":2,"regulator2":1,....}
        """
        if not isinstance(regulators, list):
            #you could give one string as a regulator
            regulators = [regulators]
        self.cooperativity = cooperativity
        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]
            
        #after we've sanitized the inputs, then sort
        self.regulators = sorted(self.regulators)
        #now let's work out the tx_capable_list
        if(tx_capable_list == None):
            #if nothing is passed, that means everything transcribes
            allcomb = []
            for r in range(1,len(self.regulators)+1):
                #make all combinations of regulators
                allcomb += [set(a) for a in it.combinations([a.name for a in self.regulators],r)]
            self.tx_capable_list = allcomb
        elif(type(tx_capable_list)==list):
            #if the user passed a list then the user knows what they want
            newlist = []
            #this part converts any species in the tx_capable_list into a string
            for element in tx_capable_list:
                sublist = []
                for specie in element:
                    if(isinstance(specie,Species)):
                        sublist += [specie.name]
                    else:
                        sublist += [specie]
                newlist+=[sublist]
            self.tx_capable_list = [set(a) for a in newlist]

        self.leak = leak
        
        self.default_mechanisms = {"binding": Combinatorial_Cooperative_Binding()}

        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters,
                          **keywords)
        self.complex_combinations = {}
        self.tx_capable_complexes = []
        self.leak_complexes = []
    def update_species(self):

        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']
        #set the tx_capable_complexes to nothing because we havent updated species yet!
        self.tx_capable_complexes = []
        self.leak_complexes = []
        species = [self.assembly.dna]
        self.complexes = []
        bound_species = mech_b.update_species(self.regulators,self.assembly.dna,\
                        component = self,part_id = self.name,cooperativity=self.cooperativity)
        #above is all the species with DNA bound to regulators. Now, we need to extract only the ones which
        #are transcribable

        if self.leak is not False:
            #this part takes care of the promoter not bound to anything
            species += mech_tx.update_species(dna = self.assembly.dna, transcript = self.transcript, protein = self.assembly.protein, component = self, part_id = self.name+"_leak")
            #self.leak_complexes += []

        for bound_complex in bound_species: 
            species_inside = []
            for regulator in self.regulators:
                if(regulator in bound_complex.species):
                    species_inside += [regulator.name] 
            if(set(species_inside) in [set(a) for a in self.tx_capable_list]):
                #only the transcribable complexes in tx_capable_list get transcription reactions
                tx_capable_species = mech_tx.update_species(dna = bound_complex, transcript = self.transcript, protein = self.assembly.protein, component = self, part_id = self.name)
                species +=tx_capable_species[1:]
                self.tx_capable_complexes +=[bound_complex]
            else:
                #in this case there's a combination of regulators which does not feature in tx_capable_list
                #this means: 
                # 1) this complex does nothing so don't add it
                # 2) we said we wanted leak, so then you should add this, but with the "_leak" parameters
                #                                                         (that happens in update_reactions)
                if(self.leak is not False):
                    leak_species = mech_tx.update_species(dna = bound_complex, transcript = self.transcript, protein = self.assembly.protein, component = self, part_id = self.name+"_leak")
                    species += leak_species[1:]
                    self.leak_complexes += [bound_complex]
        species+=bound_species  
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak is not False:
            #once again the DNA not bound to anything gets special treatment
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, component = self, part_id = self.name+"_leak", \
                                                            transcript = self.transcript, protein = self.assembly.protein)
        #the binding reactions happen no matter what
        reactions += mech_b.update_reactions(self.regulators,self.assembly.dna,component = self,\
                                                        part_id = self.name,cooperativity=self.cooperativity)
        if((self.tx_capable_complexes == None) or self.tx_capable_complexes == []):
            #this could mean we haven't run update_species() yet
            species = self.update_species()
            if(self.tx_capable_complexes == []):
                if(self.leak_complexes is None or self.leak_complexes == []):
                    #if it's still zero after running update_species then we could be in trouble
                    warn("nothing can transcribe from combinatorial promoter {}".format(self.name))
            
        if(len(self.tx_capable_complexes)>0):
            for specie in self.tx_capable_complexes:
                tx_partid = self.name
                for part in specie.species_set:
                    #construct the name of the promoter with regulators bound
                    if part.material_type == "dna":
                        #the DNA doesn't matter
                        pass
                    else:
                        #put in the regulators!
                        tx_partid += "_"+part.name
                if(tx_partid[0]=="_"):
                    #this will only happen if the name of the dna is ""
                    tx_partid = tx_partid[1:]
                #if it's bound to RNAP then it transcribes, right?
                tx_partid = tx_partid+"_RNAP"
                reactions += mech_tx.update_reactions(dna = specie, component = self, part_id = tx_partid, \
                                            transcript = self.transcript, protein = self.assembly.protein)
        if(len(self.leak_complexes)>0):
            for specie in self.leak_complexes:
                #in this case every reaction uses the "promoter_leak" partid
                leak_partid = self.name+"_leak"
                reactions += mech_tx.update_reactions(dna = specie, component = self, part_id = leak_partid, \
                                            transcript = self.transcript, protein = self.assembly.protein)


        return reactions