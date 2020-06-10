#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .component import Component
from .components_basic import DNA, RNA, Protein
from .mechanisms_binding import Reversible_Bimolecular_Binding
from .chemical_reaction_network import Species


class guideRNA(RNA):
    def __init__(self, guide_name, dCas9 = "dCas9", **keywords):

        if isinstance(dCas9, Species):
            self.dCas = dCas9
        elif isinstance(dCas9, str):
            self.dCas = Species(dCas9, material_type ="protein")
        elif isinstance(dCas9, Component) and dCas9.get_species() is not None:
            self.dCas = dCas9.get_species()
        else:
            raise ValueError("dCas9 parameter must be a "
                             "chemical_reaction_network.species, Component "
                             "with get_species(), or a string")


        self.default_mechanisms = {
            "dCas9_binding" :
                Reversible_Bimolecular_Binding(name = "dCas9_binding",
                                        mechanism_type = "bimolecular binding")
        }

        RNA.__init__(self, name = guide_name, **keywords)
        self.gRNA = self.get_species()

    def get_dCasComplex(self):
        binding_species = \
           self.mechanisms['dCas9_binding'].update_species(self.gRNA, self.dCas)
        if len(binding_species) > 1:
            raise ValueError("dCas9_binding mechanisms "
                            f"{self.mechanisms['dCas9_binding'].name} returned "
                            "multiple complexes. Unclear which is active." )
        else:
            return binding_species[0]

    def update_species(self):
        species = [self.gRNA, self.dCas]
        species += self.mechanisms['dCas9_binding'].update_species(self.gRNA, self.dCas, component = self, part_id = self.name)
        return species

    def update_reactions(self):
        rxns = self.mechanisms['dCas9_binding'].update_reactions(self.gRNA, self.dCas, component = self, part_id = self.name)
        return rxns
