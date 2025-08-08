
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex, Species


class One_Step_Reversible_Conformation_Change(Mechanism):
    """A reaction where n binders (A) bind to 1 bindee (B) in one step
       n A + B <--> nA:B
    """

    def __init__(self, name="one_step_conformation_change",
                 mechanism_type="conformation_change"):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, s0, sf, additional_species=None, component=None, part_id=None, **kwords):

        if additional_species is None:
            additional_species = []

        return [s0, sf] + additional_species

    def update_reactions(self, s0, sf, additional_species=None, component=None, part_id=None, **kwords):

        if part_id is None:
            repr(s0)+"-"+repr(sf)

        if additional_species is None:
            additional_species = []

        kf = component.get_parameter("kf", part_id=part_id, mechanism=self)
        kr = component.get_parameter("kr", part_id=part_id, mechanism=self)

        return [Reaction.from_massaction(inputs=[s0]+additional_species, outputs=[sf], k_forward=kf, k_reverse=kr)]
