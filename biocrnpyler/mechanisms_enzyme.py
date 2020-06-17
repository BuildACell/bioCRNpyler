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
        return [Enzyme, Sub, Prod]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, kcat = None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)

        return [Reaction([Enzyme, Sub], [Enzyme, Prod], kcat)]

class MichalisMenten(Mechanism):
    """Mechanism to automatically generate Michalis-Menten Type Reactions
       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz+Prod
    """

    def __init__(self, name = "michalis_menten", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod = None, **keywords):
        complex = ComplexSpecies([Sub, Enzyme])
        return [Enzyme, Sub, complex]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        if component == None and (kb == None or ku == None or kcat == None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku is None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if kcat is None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        

        if complex == None:
            complex = ComplexSpecies([Sub, Enzyme])

        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction(inputs=[Sub, Enzyme], outputs=[complex],
                               k=kb, k_rev=ku)
        if Prod is not None:
            # Sub:Enz --> Enz + Prod
            cat_rxn = Reaction(inputs=[complex],
                               outputs=[Prod, Enzyme], k=kcat)
        else:  # Degradation Reaction
            # Sub:Enz --> Enz
            cat_rxn = Reaction(inputs=[complex], outputs=[Enzyme],
                               k=kcat)
        return [binding_rxn, cat_rxn]


class MichalisMentenReversible(Mechanism):
    """Mechanism to automatically generate Michalis-Menten Type Reactions with products that can bind to enzymes
       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz:Prod <--> Enz + Prod
    """

    def __init__(self, name = "michalis_menten_reverse_binding", mechanism_type = "catalysis", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod, **keywords):
        complex1 = ComplexSpecies([Sub, Enzyme])
        complex2 = ComplexSpecies([Prod, Enzyme])
        return [Enzyme, Sub, Prod, complex1, complex2]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, complex2 = None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        if component == None and (kb == None or ku == None or kcat == None):
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
        

        if complex == None:
            complex = ComplexSpecies([Sub, Enzyme])
        if complex2 == None:
            complex2 = ComplexSpecies([Prod, Enzyme])

        # Sub + Enz <--> Sub:Enz
        binding_rxn1 = Reaction(inputs=[Sub, Enzyme], outputs=[complex],
                               k=kb1, k_rev=ku1)

        binding_rxn2 = Reaction(inputs=[Prod, Enzyme], outputs=[complex2],
                               k=kb2, k_rev=ku2)

        # Sub:Enz --> Enz:Prod
        cat_rxn = Reaction(inputs=[complex], outputs=[complex2], k=kcat, k_rev = kcat_rev)
        
        return [binding_rxn1, binding_rxn2, cat_rxn]


class MichalisMentenCopy(Mechanism):
    """In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Sub+Enz+Prod
    """
    def __init__(self, name = "michalis_menten_copy", mechanism_type = "copy", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod = None, **keywords):
        complex = ComplexSpecies([Sub, Enzyme])
        return [complex]

    def update_reactions(self, Enzyme, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        if complex == None:
            complex = ComplexSpecies([Sub, Enzyme])

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        if kb == None and component != None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku == None and component != None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if kcat == None and component != None:
            kcat = component.get_parameter("kcat", part_id = part_id, mechanism = self)
        if component == None and (kb == None or ku == None or kcat == None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")
        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction(inputs=[Sub, Enzyme], outputs=[complex],
                               k=kb, k_rev=ku)

        # Sub:Enz --> Enz + Prod + Sub
        cat_rxn = Reaction(inputs=[complex], outputs=[Sub, Prod, Enzyme],
                           k=kcat)

        return [binding_rxn, cat_rxn]

