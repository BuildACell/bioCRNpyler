
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from .chemical_reaction_network import Species
from .component import Component
from .dna_part import DNA_part


class CDS(DNA_part):
    def __init__(self,name,protein,no_stop_codons=None, **keywords):
        """a CDS is a sequence of DNA that codes for a protein"""
        self.name = name
        DNA_part.__init__(self,name,no_stop_codons=no_stop_codons,protein=None, **keywords)
        #TODO make this contain a species and not a component
        #TODO use set_species()
        if(protein is None):
            self.protein = Species(name,material_type="protein")
        elif(isinstance(protein,str)):
            self.protein = Species(protein,material_type="protein")
        elif(isinstance(protein,Component)):
            self.protein=protein.get_species()
    def update_species(self):
        return [self.protein]
    def update_reactions(self):
        return []
    def get_species(self):
        return self.protein
