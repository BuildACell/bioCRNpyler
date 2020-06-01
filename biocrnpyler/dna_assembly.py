#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .component import Component, DNA
from .chemical_reaction_network import ComplexSpecies, Species
from .mechanism import One_Step_Cooperative_Binding
from warnings import warn as pywarn
import itertools as it
import numpy as np

def warn(txt):
    pywarn(txt)

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
            transcript = self.transcript, protein = self.assembly.protein)
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

            species_b = mech_b.update_species(regulator, self.assembly.dna, component = self, part_id = self.name+"_"+regulator.name)
            species += species_b

            #Find complexes containing DNA and the regulator
            for s in species_b:
                if isinstance(s, ComplexSpecies) and self.assembly.dna in s.species and regulator in s.species:
                    self.complexes += [s]

            species += mech_tx.update_species(dna = complex_, transcript = self.transcript, protein = self.assembly.protein, part_id = self.name+"_"+regulator.name)
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

class RBS(Component):
    def __init__(self, name, assembly = None,
                 transcript = None, protein = None, length = 0,
                 mechanisms = {}, parameters = {}, **keywords):
        self.assembly = assembly
        self.length = length

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type = "rna")
        else:
            self.transcript = self.set_species(transcript, material_type = "rna")
 
        if protein is None:
            self.protein = Species(assembly.name, material_type = "protein")
        else:
            self.protein = self.set_species(protein, material_type = "protein")

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tl = self.mechanisms['translation']
        species = []
        species += mech_tl.update_species(transcript = self.transcript, protein = self.protein)
        return species

    def update_reactions(self):
        mech_tl = self.mechanisms['translation']
        reactions = []

        reactions += mech_tl.update_reactions(transcript = self.transcript,
                                              protein = self.protein, component = self, part_id = self.name)
        return reactions


class DNAassembly(DNA):
    def __init__(self, name, dna = None, promoter = None, transcript = None,
                 rbs = None, protein = None, length = None,
                 attributes = [], mechanisms = {}, parameters = {}, initial_conc = None,
                 parameter_warnings = True, **keywords):
        self.promoter = None
        self.rbs = None
        self.transcript = None
        self.initial_concentration = initial_conc
        self.name = name

        DNA.__init__(self, name, length = length, mechanisms = mechanisms,
                     parameters = parameters, initial_conc = initial_conc,
                     parameter_warnings = parameter_warnings,
                     attributes = list(attributes), **keywords)

        self.update_dna(dna, attributes = list(attributes))
        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_promoter(promoter, transcript = self.transcript)
        self.update_rbs(rbs, transcript = self.transcript,
                        protein = self.protein)

        self.set_parameter_warnings(parameter_warnings)

            

    def set_parameter_warnings(self, parameter_warnings):
        self.parameter_warnings = parameter_warnings

        if self.parameter_warnings != None:
            if self.promoter != None:
                self.promoter.set_parameter_warnings(parameter_warnings)
            if self.rbs != None:
                self.rbs.set_parameter_warnings(parameter_warnings)


    def update_dna(self, dna, attributes = None):
        if dna is None:
            self.dna = self.set_species(self.name, material_type = "dna", attributes = attributes)
        else:
            self.dna = self.set_species(dna, material_type = "dna", attributes = attributes)
        

    def update_transcript(self, transcript, attributes = None):
        if transcript is None:
            self.transcript = self.set_species(self.name, material_type = "rna", attributes = attributes)
        else:
            self.transcript = self.set_species(transcript, material_type = "rna", attributes = attributes)

        if self.promoter is not None:
            self.promoter.transcript = self.transcript
        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(self, protein, attributes = None):
        if protein is None:
            self._protein = self.set_species(self.name, material_type = "protein", attributes = attributes)
        else:
            self._protein = self.set_species(protein, material_type = "protein", attributes = attributes)

        if self.rbs is not None:
            self.rbs.transcript = self.protein

    def update_promoter(self, promoter, transcript=None):
        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(promoter, str):
            self.promoter = Promoter(assembly = self, name = promoter,
                                     transcript = self.transcript,
                                     parameters = self.parameters)
        elif isinstance(promoter, Promoter):
            self.promoter = promoter
            self.promoter.assembly = self
            self.promoter.transcript = self.transcript
        elif promoter is not None:
            raise ValueError("Improper promoter type recieved by DNAassembly. "
                             "Expected string or promoter object. "
                             f"Recieved {repr(promoter)}.")
        if promoter is not None:
            self.promoter.update_parameters(
                                        mixture_parameters = self.parameters,
                                        overwrite_custom_parameters = False)

    def update_rbs(self, rbs, transcript = None, protein = None):
        if protein is not None:
            self.update_protein(protein)

        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(rbs, str):
            self.rbs = RBS(assembly = self, name = rbs, protein = self._protein,
                           transcript = self.transcript,
                           parameters = self.parameters)
        elif isinstance(rbs, RBS):
            self.rbs = rbs
            self.rbs.assembly = self
            self.rbs.transcript = self.transcript
            self.rbs.protein = self._protein
        elif rbs is not None:
            raise ValueError("Improper rbs type recieved by DNAassembly. "
                             "Expected string or RBS object. Recieved "
                            f"{repr(rbs)}.")

        if rbs is not None:
            self.rbs.update_parameters(mixture_parameters = self.parameters,
                                       overwrite_custom_parameters = False)

    @property
    def protein(self):
        return self._protein

    def update_species(self):
        species = []
        species.append(self.dna)
        if self.promoter is not None and self.rbs is not None:
            species += self.promoter.update_species()
            species += self.rbs.update_species()

        elif self.promoter is not None and self.rbs is None:
            species += self.promoter.update_species()

        if "rna_degredation" in self.mechanisms and self.promoter is not None:
            deg_mech = self.mechanisms["rna_degredation"]
            species += deg_mech.update_species(rna = self.transcript)

        # TODO raise a warning if there were duplicate species
        return list(set(species))

    def update_reactions(self):
        reactions = []
        if self.promoter is not None:
            self.promoter.parameter_warnings = self.parameter_warnings
            reactions += self.promoter.update_reactions()

        if self.rbs is not None:
            self.rbs.parameter_warnings = self.parameter_warnings
            reactions += self.rbs.update_reactions()

        if "rna_degredation" in self.mechanisms and self.promoter is not None:
            deg_mech = self.mechanisms["rna_degredation"]


            reactions += deg_mech.update_reactions(rna = self.transcript, component = self.promoter)
        # TODO check that the reaction list is unique
        return reactions

    def update_parameters(self, mixture_parameters = {}, parameters = {},
                          overwrite_custom_parameters = True):
        DNA.update_parameters(self = self,
                              mixture_parameters = mixture_parameters,
                              parameters = parameters)

        if self.promoter is not None:
            self.promoter.update_parameters(
                                    mixture_parameters = mixture_parameters,
                                    parameters = parameters,
                                    overwrite_custom_parameters = False)
        if self.rbs is not None:
            self.rbs.update_parameters(mixture_parameters = mixture_parameters,
                                       parameters = parameters,
                                       overwrite_custom_parameters = False)

    def update_mechanisms(self, mixture_mechanisms = {}, mechanisms = {},
                          overwrite_custom_parameters = True):
        DNA.update_mechanisms(self = self,
                              mixture_mechanisms = mixture_mechanisms,
                              mechanisms = mechanisms)

        if self.promoter is not None and "transcription" in self.mechanisms:
            mech_tx = self.mechanisms["transcription"]
            mechs = {"transcription": mech_tx}
            self.promoter.update_mechanisms(mechanisms = mechs,
                                            overwrite_custom_mechanisms = False)
        if self.rbs is not None and "translation" in self.mechanisms:
            mech_tl = self.mechanisms["translation"]
            mechs = {"translation": mech_tl}
            self.rbs.update_mechanisms(mechanisms = mechs,
                                       overwrite_custom_mechanisms = False)

    def __str__(self):
        return type(self).__name__ + ": " + self.name

    def __repr__(self):
        txt = str(self)
        if self.promoter is not None:
            txt += "\n\t" + repr(self.promoter)
            txt += "\n\ttranscript = " + repr(self.transcript)
        if self.rbs is not None:
            txt += "\n\t" + repr(self.rbs)
            txt += "\n\tprotein = " + repr(self._protein)
        return txt