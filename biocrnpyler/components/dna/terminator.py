
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.from .dna_part import DNA_part

from .misc import DNA_part


class Terminator(DNA_part):
    def __init__(self,name, **keywords):
        DNA_part.__init__(self,name, **keywords)
        self.name = name
    def update_species(self):
        return []
    def update_reactions(self):
        return []
