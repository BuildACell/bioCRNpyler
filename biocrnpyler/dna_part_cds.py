from .component import Component
from .chemical_reaction_network import Species
from .components_basic import Protein
from .dna_part import DNA_part
class CDS(DNA_part):
    def __init__(self,name,protein,no_stop_codons=[], **keywords):
        """a CDS is a sequence of DNA that codes for a protein"""
        self.name = name
        DNA_part.__init__(self,name,no_stop_codons=no_stop_codons,protein=None, **keywords)
        #TODO make this contain a species and not a component
        #TODO use set_species()
        if(protein is None):
            self.protein = Protein(name)
        
        elif(isinstance(protein,str)):
            self.protein = Protein(name)
        elif(isinstance(protein,Component)):
            self.protein=protein
        elif(isinstance(protein,Species)):
            raise ValueError("CDS got Species as 'protein', but it should be a Component")
    def update_species(self):
        return [self.protein]
    def update_reactions(self):
        return []
    def get_species(self):
        return self.protein.get_species()