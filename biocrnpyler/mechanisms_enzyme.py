from warnings import warn
from .mechanism import *
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer


class BasicCatalysis(Mechanism):
    """
    Mechanism for the schema S + C --> P + C
    """
    def __init__(self, name = "basic_catalysis", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod = None, **keywords):
        if Prod is None:
            return [Enzyme, Sub]
        else:
            return [Enzyme, Sub, Prod]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, kcat = None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)

        return [Reaction([Enzyme, Sub], [Enzyme, Prod], kcat)]

class BasicProduction(Mechanism):
    """
    Mechanism for the schema C --> P + C
    """
    def __init__(self, name = "basic_production", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub = None, Prod = None, **keywords):
        species = [Enzyme]
        if Prod is not None:
            species += [Prod]
        if Sub is not None:
            species += [Sub]

        return species

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, kcat = None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)

        inputs = [Enzyme]
        outputs = [Enzyme]
        if Prod is not None:
            outputs += [Prod]
        if Sub is not None:
            inputs += [Sub]

        return [Reaction(inputs, outputs, kcat)]

class MichaelisMenten(Mechanism):
    """Mechanism to automatically generate Michaelis-Menten Type Reactions
       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz+Prod
    """

    def __init__(self, name = "michalis_menten", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod = None, complex=None, **keywords):
        if complex is None:
            complexS = ComplexSpecies([Sub, Enzyme])
        else:
            complexS = complex
        if Prod is None:
            return [Enzyme, Sub, complexS]
        else:
            return [Enzyme, Sub, Prod, complexS]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku is None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        

        if complex is None:
            complexS = ComplexSpecies([Sub, Enzyme])
        else:
            complexS = complex

        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction(inputs=[Sub, Enzyme], outputs=[complexS],
                               k=kb, k_rev=ku)
        if Prod is not None:
            # Sub:Enz --> Enz + Prod
            cat_rxn = Reaction(inputs=[complexS],
                               outputs=[Prod, Enzyme], k=kcat)
        else:  # Degradation Reaction
            # Sub:Enz --> Enz
            cat_rxn = Reaction(inputs=[complexS], outputs=[Enzyme],
                               k=kcat)
        return [binding_rxn, cat_rxn]


class MichaelisMentenReversible(Mechanism):
    """Mechanism to automatically generate Michaelis-Menten Type Reactions with products that can bind to enzymes
       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz:Prod <--> Enz + Prod
    """

    def __init__(self, name = "michalis_menten_reverse_binding", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod, complex=None, complex2 = None, **keywords):
        if complex is None:
            complex1 = ComplexSpecies([Sub, Enzyme])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = ComplexSpecies([Prod, Enzyme])
        else:
            complex2 = complex2
        return [Enzyme, Sub, Prod, complex1, complex2]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, complex2 = None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
            kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
        else:
            kb1, kb2 = kb
        if ku is None:
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
        else:
            ku1, ku2 = ku
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
            kcat_rev = component.get_parameter("kcat_rev", part_id = part_id, mechanism = self)
        else:
            kcat, kcat_rev = kcat
        

        if complex is None:
            complex1 = ComplexSpecies([Sub, Enzyme])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = ComplexSpecies([Prod, Enzyme])

        # Sub + Enz <--> Sub:Enz
        binding_rxn1 = Reaction(inputs=[Sub, Enzyme], outputs=[complex1],
                               k=kb1, k_rev=ku1)

        binding_rxn2 = Reaction(inputs=[Prod, Enzyme], outputs=[complex2],
                               k=kb2, k_rev=ku2)

        # Sub:Enz --> Enz:Prod
        cat_rxn = Reaction(inputs=[complex1], outputs=[complex2], k=kcat, k_rev = kcat_rev)
        
        return [binding_rxn1, binding_rxn2, cat_rxn]


class MichaelisMentenCopy(Mechanism):
    """In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Sub+Enz+Prod
    """
    def __init__(self, name = "michalis_menten_copy", mechanism_type = "copy", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, complex=None, Prod = None, **keywords):
        if complex is None:
            complexS = ComplexSpecies([Sub, Enzyme])
        else:
            complexS = complex
            
        if Prod is None:
            return [Enzyme, Sub, complexS]
        else:
            return [Enzyme, Sub, Prod, complexS]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        if complex is None:
            complexS = ComplexSpecies([Sub, Enzyme])
        else:
            complexS = complex

        #Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if kb is None and component is not None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku is None and component is not None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if kcat is None and component is not None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction(inputs=[Sub, Enzyme], outputs=[complexS],
                               k=kb, k_rev=ku)

        # Sub:Enz --> Enz + Prod + Sub
        cat_rxn = Reaction(inputs=[complexS], outputs=[Sub, Prod, Enzyme],
                           k=kcat)

        return [binding_rxn, cat_rxn]

