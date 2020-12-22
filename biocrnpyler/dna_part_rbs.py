
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy

from .dna_part import DNA_part
from .species import Species



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
    def update_component(self,internal_species=None,**keywords):
        """returns a copy of this component, except with the proper fields updated"""
        from .dna_construct import RNA_construct,DNA_construct
        if(isinstance(self.parent,DNA_construct)):
            return None
        elif(isinstance(self.parent,RNA_construct)):
            if(self.direction=="forward"):
                out_component = copy.copy(self)
                out_component.transcript = internal_species
                return out_component
            elif(self.direction=="reverse"):
                return None
            else:
                raise AttributeError(f"Unknown direction {self.direction} encountered in {self}")
        else:
            return None
        return out_component

    @classmethod
    def from_rbs(cls, name, assembly, transcript, protein):
        """Helper function to initialize a rbs instance from another rbs or str.

        :param name: either string or an other rbs instance
        :param assembly:
        :param transcript:
        :param protein:
        :return: RBS instance
        """
        if isinstance(name, RBS):
            rbs_instance = copy.deepcopy(name)
            rbs_instance.assembly = assembly
            rbs_instance.transcript = transcript
            rbs_instance.protein = protein
        elif isinstance(name, str):
            rbs_instance = cls(name=name, assembly=assembly,
                               transcript=transcript, protein=protein)
        else:
            raise TypeError(f'RBS can be initialized from string or another RBS! We got {type(name)}')
        return  rbs_instance
