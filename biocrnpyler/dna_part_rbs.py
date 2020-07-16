from .component import Component
from .components_basic import DNA, RNA, Protein
from .chemical_reaction_network import ComplexSpecies, Species
from .mechanisms_binding import *
from .mechanisms_txtl import *
from .dna_part import DNA_part
import copy
class RBS(DNA_part):
    def __init__(self, name: str, assembly = None,
                 transcript = None, protein = None, length = 0,
                 mechanisms = {}, parameters = {}, **keywords):
        self.assembly = assembly
        self.length = length

        DNA_part.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type = "rna")
        else:
            self.transcript = self.set_species(transcript, material_type = "rna")
 
        if protein is None and assembly is None:
            self.protein = None
        elif protein is None:
            self.protein = Species(assembly.name, material_type = "protein")
        else:
            self.protein = self.set_species(protein, material_type = "protein")

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tl = self.mechanisms['translation']
        species = []
        species += mech_tl.update_species(transcript = self.transcript, protein = self.protein, component = self, part_id = self.name)
        return species
    def get_species(self):
        return None
    def update_reactions(self):
        mech_tl = self.mechanisms['translation']
        reactions = []

        reactions += mech_tl.update_reactions(transcript = self.transcript,
                        protein = self.protein, component = self, part_id = self.name)
        return reactions
    def update_component(self,dna=None,rnas=None,proteins=None,mypos=None):
        """returns a copy of this component, except with the proper fields updated"""
        if(proteins is None):
            ValueError("cannot update rbs {} when proteins is None".format(self))
        my_rna = None
        if(len(proteins)==1):
            my_rna = list(proteins.keys())[0]
        else:
            #if this happens, it means we are being called from DNA
            if(self.parent.get_species().material_type is not "dna"):
                raise ValueError("something went wrong")
            return None
            
            #for rna in proteins.keys():
            #    if(self in rna):
            #        my_rna = rna
            #if(my_rna is None):
            #    raise KeyError(str(self) + " couldn't be found in "+str(list(proteins.keys())))
        print("my rna is "+str(my_rna))
        out_component = None
        if(self in proteins[my_rna]):
            my_proteins = proteins[my_rna][self]
            out_component = copy.deepcopy(self)
            protein_species = []
            for CDS in my_proteins:
                protein_species += CDS.update_species()

            out_component.protein = protein_species
            if(mypos is not None):
                out_component.transcript = dna[mypos]
            else:
                out_component.transcript = dna
            #my_rna.get_species()
        return out_component
