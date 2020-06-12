from .components_basic import Protein
from .dna_part import DNA_part
class CDS(DNA_part):
    def __init__(self,name,protein,no_stop_codons=[], **keywords):
        """a CDS is a sequence of DNA that codes for a protein"""
        self.name = name
        DNA_part.__init__(self,name,no_stop_codons=no_stop_codons, **keywords)
        self.protein = Protein(name)
    def update_species(self):
        return [self.protein]
    def update_reactions(self):
        return []
    def get_species(self):
        return self.protein.get_species()