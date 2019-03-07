# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from .component import Component, DNA
from .chemical_reaction_network import Species
from .mechanism import One_Step_Cooperative_Binding
from warnings import warn as pywarn


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
        elif isinstance(transcript, str):
            self.transcript = Species(transcript, material_type="rna")
        elif isinstance(transcript, Species):
            self.transcript = transcript
        else:
            raise ValueError(
                "Improper value of transcriped passed into promoter. "
                "Transcript should be a string or "
                "chemical_reaction_network.species")

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        species = [self.transcript]
        species += mech_tx.update_species(self.assembly.dna)
        return species

    def update_reactions(self):
        mech_tx = self.mechanisms["transcription"]
        reactions = []
        ktx = self.get_parameter("ktx", part_id = self.name, mechanism = mech_tx)
        kb = self.get_parameter("kb", part_id = self.name, mechanism = mech_tx)
        ku = self.get_parameter("ku", part_id = self.name, mechanism = mech_tx)

        reactions += mech_tx.update_reactions(self.assembly.dna, ktx = ktx,
                                              ku = ku, kb = kb, complex = None,
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
            if isinstance(regulator, Species):
                self.regulators += [regulator]
            elif isinstance(regulator, str):
                self.regulators += [Species(name = regulator,
                                            material_type = "protein")]
            else:
                raise ValueError(
                    "Invalid parameter 'regulators' passed to "
                    "RegulatedPromoter: valid types are string, "
                    "chemical_reaction_network.species and lists of these "
                    "strings or species")

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
            species += mech_tx.update_species(dna = self.assembly.dna)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            coop = self.get_parameter(param_name = "cooperativity",
                                      part_id = f"{self.name}_{regulator.name}",
                                      mechanism = mech_b)
            species_b = mech_b.update_species(regulator, self.assembly.dna,
                                              cooperativity = coop)
            species += species_b
            complex_ = species_b[0]
            self.complexes += [complex_]
            species += mech_tx.update_species(dna = complex_)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak != False:
            ktx = self.get_parameter("ktx", mechanism = mech_tx,
                                     part_id = self.name)
            ku = self.get_parameter("ku", mechanism = mech_tx,
                                    part_id = self.name)
            kb = self.get_parameter("kb", mechanism = mech_tx,
                                    part_id = self.name)
            reactions += mech_tx.update_reactions(dna = self.assembly.dna,
                                                  ktx = ktx, ku = ku, kb = kb)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]
            ktx = self.get_parameter("ktx", mechanism = mech_tx,
                                     part_id = self.name+"_"+regulator.name)

            ku_tx = self.get_parameter("ku", mechanism = mech_tx,
                                       part_id = self.name+"_"+regulator.name)
            kb_tx = self.get_parameter("kb", mechanism = mech_tx,
                                       part_id = self.name+"_"+regulator.name)

            ku_c = self.get_parameter("ku", mechanism = mech_b,
                                      part_id = self.name+"_"+regulator.name)
            kb_c = self.get_parameter("kb", mechanism = mech_b,
                                      part_id = self.name+"_"+regulator.name)

            coop = self.get_parameter("cooperativity",
                                      part_id = f"{self.name}_{regulator.name}",
                                      mechanism = mech_b)

            reactions += mech_b.update_reactions(regulator, self.assembly.dna,
                                                 ku = ku_c, kb=kb_c,
                                                 cooperativity = coop)
            reactions += mech_tx.update_reactions(dna = complex_, kb = kb_tx,
                                                 ku = ku_tx, ktx = ktx,
                                                 transcript = self.transcript)

        return reactions


class RBS(Component):
    def __init__(self, name, assembly,
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
        elif isinstance(transcript, str):
            self.transcript = Species(transcript, material_type = "rna")
        elif isinstance(transcript, Species):
            self.transcript = transcript
        else:
            raise ValueError(
                "Improper value of transcript passed into promoter: transcript "
                "should be a string or chemical_reaction_network.species")

        if protein is None:
            self.protein = Species(assembly.name, material_type = "protein")
        elif isinstance(protein, str):
            self.protein = Species(protein, material_type = "protein")
        elif isinstance(protein, Species):
            self.protein = protein
        else:
            raise ValueError(
                "Improper value of protein passed into RBS: protein should be "
                "a string or chemical_reaction_network.species")

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tl = self.mechanisms['translation']
        species = [self.transcript, self.protein]
        species += mech_tl.update_species(self.transcript)
        return species

    def update_reactions(self):
        mech_tl = self.mechanisms['translation']
        reactions = []

        ktl = self.get_parameter("ktl", part_id = self.name,
                                 mechanism = mech_tl)
        kb = self.get_parameter("kb", part_id = self.name, mechanism = mech_tl)
        ku = self.get_parameter("ku", part_id = self.name, mechanism = mech_tl)
        reactions += mech_tl.update_reactions(transcript = self.transcript,
                                              protein = self.protein,
                                              ku = ku, kb = kb, ktl = ktl)
        return reactions


class DNAassembly(DNA):
    def __init__(self, name, dna = None,
                 promoter = None, transcript = None,
                 rbs = None, protein = None, length = None,
                 attributes = [],
                 mechanisms = {}, parameters = {}, initial_conc = None,
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
        if isinstance(dna, Species):
            self.dna = dna
        elif dna is None:
            self.dna = Species(self.name, material_type = "dna",
                            initial_concentration = self.initial_concentration,
                            attributes = attributes)
        elif isinstance(dna, str):
            self.dna = Species(dna, material_type = "dna",
                            initial_concentration = self.initial_concentration,
                            attributes = attributes)
        else:
            raise ValueError(f"Invalid value of 'dna' passed to {self}: "
                             "dna must be None, a "
                             "chemical_reaction_network.species, or a string.")

    def update_transcript(self, transcript):
        if isinstance(transcript, Species):
            self.transcript = transcript
        elif isinstance(transcript, str):
            self.transcript = Species(transcript, material_type = "rna")
        elif transcript is None:
            self.transcript = Species(self.name, material_type = "rna")
        elif isinstance(transcript, Component) \
              and transcript.get_species() is not None:
            self.transcript = transcript.get_species()
        else:
            raise ValueError(f"Invalid value of transcript passed to {self}: "
                "transcript must be None, a chemical_reaction_network.species, "
                "or a string.")

        if self.promoter is not None:
            self.promoter.transcript = self.transcript
        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(self, protein):
        if isinstance(protein, Species):
            self._protein = protein
        elif isinstance(protein, str):
            self._protein = Species(protein, material_type="protein")
        elif protein is None:
            self._protein = Species(self.name, material_type="protein")
        else:
            raise ValueError(
                f"Invalid value of 'protein' passed to DNA Assembly "
                "{self.name}. 'protein' must be None, a "
                "chemical_reaction_network.species, or a string.")

        if self.rbs is not None:
            self.rbs.transcript = self._protein

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
                            f"{repr(rbs}.")

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
            species += deg_mech.update_species(self.transcript)

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
            ku = self.get_parameter("ku", mechanism=deg_mech)
            kb = self.get_parameter("kb", mechanism=deg_mech)
            kdeg = self.get_parameter("kdeg", mechanism=deg_mech)

            reactions += deg_mech.update_reactions(self.transcript, ku = ku,
                                                    kb = kb, kdeg = kdeg)
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
