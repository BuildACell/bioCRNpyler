# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ...core.chemical_reaction_network import Species
from ...core.component import Component
from .construct import DNA_part


class CDS(DNA_part):
    def __init__(self, name, protein=None, no_stop_codons=[], **kwargs):
        """a CDS is a sequence of DNA that codes for a protein"""
        self.name = name
        DNA_part.__init__(self, name, no_stop_codons=no_stop_codons, **kwargs)
        # TODO use set_species()
        if protein is None:
            self.protein = Species(name, material_type='protein')
        elif isinstance(protein, str):
            self.protein = Species(protein, material_type='protein')
        elif isinstance(protein, Component):
            self.protein = protein.get_species()

    def update_species(self):
        return [self.protein]

    def update_reactions(self):
        return []

    def get_species(self):
        return self.protein
