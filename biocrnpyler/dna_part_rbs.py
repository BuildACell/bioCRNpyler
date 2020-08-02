from .component import Component
from .species import Species
from .dna_part import DNA_part
import copy
class RBS(DNA_part):
    """
    A simple RBS class with no regulation. Must be included in a DNAconstruct or DNAassembly to do anything.
    """
    def __init__(self, name: str, assembly=None,
                 transcript=None, protein=None, length=0,
                 mechanisms=None, parameters=None, **keywords):
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
        
        if protein is None and assembly is not None:
            self.protein = Species(assembly.name, material_type = "protein")
        elif protein is None and assembly is None:
            self.protein = None
        else:
            self.protein = self.set_species(protein, material_type = "protein")

    def update_species(self):
        mech_tl = self.get_mechanism('translation')
        species = []
        if self.protein is not None:
            species += mech_tl.update_species(transcript = self.transcript, protein = self.protein, component = self, part_id = self.name)
        return species

    def update_reactions(self):
        mech_tl = self.get_mechanism('translation')
        reactions = []

        if self.protein is not None:
            reactions += mech_tl.update_reactions(transcript = self.transcript, protein = self.protein, component = self, part_id = self.name)
        return reactions
    def update_component(self,dna=None,rnas=None,proteins=None,mypos=None):
        """returns a copy of this component, except with the proper fields updated"""
        if(proteins is None):
            raise ValueError("cannot update rbs {} when proteins is None".format(self))
        my_rna = None
        if(len(proteins)==1):
            my_rna = list(proteins.keys())[0]
        else:
            #if this happens, it means we are being called from DNA
            if(self.parent.get_species().material_type is not "dna"):
                raise ValueError("something went wrong")
            return None
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
