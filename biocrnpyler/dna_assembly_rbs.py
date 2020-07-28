from .component import Component
from .species import Species

class RBS(Component):
    """
    A simple RBS class with no regulation. Must be included in a DNAconstruct or DNAassembly to do anything.
    """
    def __init__(self, name: str, assembly=None,
                 transcript=None, protein=None, length=0,
                 mechanisms=None, parameters=None, **keywords):
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
