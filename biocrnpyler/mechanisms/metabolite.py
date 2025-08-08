from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex
from ..core.propensities import HillPositive


# [precursors] --> [products] using Massaction (None OK)
class OneStepPathway(Mechanism):
    def __init__(self, name="one_step_pathway", mechanism_type="metabolic_pathway"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, precursor, product, **keywords):
        species = []
        if precursor is not None:
            species += precursor
        if product is not None:
            species += product
        return species

    def update_reactions(self, precursor, product, component=None, part_id=None, k=None, **keywords):
        if precursor is None:
            inputs = []
        else:
            inputs = precursor

        if product is None:
            outputs = []
        else:
            outputs = product

        if component is None and k is None:
            raise ValueError("Must pass in a component or a rate k.")
        elif k is None:
            k = component.get_parameter("k", part_id=part_id, mechanism=self)

        r = Reaction.from_massaction(inputs, outputs, k_forward=k)
        return [r]
