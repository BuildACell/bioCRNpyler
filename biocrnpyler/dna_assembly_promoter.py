from .component import Component
from .components_basic import DNA, RNA, Protein, ChemicalComplex
from .chemical_reaction_network import ComplexSpecies, Species
from .mechanisms_binding import *
from .mechanisms_txtl import *
from warnings import warn as pywarn
import itertools as it
import numpy as np

class Promoter(Component):
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

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        species = []
        species += mech_tx.update_species(dna = self.assembly.dna, \
            transcript = self.transcript, protein = self.assembly.protein,
            component = self, part_id = self.name)
        return species

    def update_reactions(self):
        mech_tx = self.mechanisms["transcription"]
        reactions = []


        reactions += mech_tx.update_reactions(dna = self.assembly.dna, \
                        component = self, part_id = self.name, complex = None,
                        transcript = self.transcript, protein = self.assembly.protein)
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

                    species += mech_tx.update_species(dna = s, transcript = self.transcript, protein = self.assembly.protein, part_id = self.name+"_"+regulator.name, component = self)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak != False:
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, component = self, part_id = self.name, \
                                                            transcript = self.transcript, protein = self.assembly.protein)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]

            reactions += mech_b.update_reactions(regulator, self.assembly.dna, component = self, \
                                                                    part_id = self.name+"_"+regulator.name)
            reactions += mech_tx.update_reactions(dna = complex_, component = self, part_id = self.name+"_"+regulator.name, \
                                            transcript = self.transcript, protein = self.assembly.protein)

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