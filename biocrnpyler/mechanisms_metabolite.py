from .mechanism import Mechanism
from .reaction import Reaction
from .species import Complex
from .propensities import HillPositive


# precursor --> product using Massaction (None OK)
class OneStepPathway(Mechanism):
    def __init__(self, name = "one_step_pathway", mechanism_type = "metabolic_pathway"):
        Mechanism.__init__(self, name = name, mechanism_type = mechanism_type)


    def update_species(self, precursor, product, **keywords):
        species = []
        if precursor is not None:
            species += [precursor]
        if product is not None:
            species += [product]
        return species

    def update_reactions(self, precursor, product, component = None, part_id = None, k = None, **keywords):
        if precursor is None:
            inputs = []
        else:
            inputs = [precursor]

        if product is None:
            outputs = []
        else:
            outputs = [product]

        if component is None and k is None:
            raise ValueError("Must pass in a component or a rate k.")
        elif k is None:
            k = component.get_parameter("k", part_id = part_id, mechanism = self)

        r = reaction.from_massaction(inputs, outputs, k_forward = k)
        return [r]


# precursor --> product using PositiveHillFunction of X (Y = 0 OK. If precursor is None, defaults to a massaction propensity)
class OneStepHillPathway(Mechanism):
    def __init__(self, name = "one_step_hill_pathway", mechanism_type = "metabolic_pathway"):
        Mechanism.__init__(self, name = name, mechanism_type = mechanism_type)

    def update_species(self, precursor, product, **keywords):
        species = []
        if precursor is not None:
            species += [precursor]
        if product is not None:
            species += [product]
        return species

    def update_reactions(self, precursor, product, component = None, part_id = None, vmax = None, K = None, n = 1, **keywords):
        
        if (vmax is None or K is None) and component is None:
            raise ValuError("Must recieve component keyword or vmax and K keywords.")

        if vmax is None:
            vmax = component.get_parameter("vmax", part_id = part_id, mechanism = self)

        if K is None:
            K = component.get_parameter("K", part_id = part_id, mechanism = self)

        if n is None:
            n = component.get_parameter("n", part_id = part_id, mechanism = self)

        if product is None:
            outputs = []
        else:
            outputs = [product]

        if precursor is None:
            inputs = []
            #default to massaction if there is no input
            r = Reaction.from_massaction(inputs, outputs, k_forward = vmax)
        else:
            inputs = [precursor]
            prop = HillPositive(k = vmax, K = K, n = n, s1 = precursor)
            r = Reaction(inputs, outputs, propensity_type = prop)

        return [r]



