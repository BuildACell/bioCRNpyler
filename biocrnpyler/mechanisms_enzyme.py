from warnings import warn
from .mechanism import *
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer

class MichalisMentenRXN(Mechanism):
    """Helper class to automatically generate Michalis-Menten Type Reactions
       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz+Prod
    """

    def __init__(self, name, enzyme, mechanism_type = "catalysis", **keywords):
        if isinstance(enzyme, Species):
            self.Enzyme = enzyme
        else:
            raise ValueError("MichalisMentenRXN takes a species object for its "
                             "enzyme argument.")

        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Sub, **keywords):
        complex = ComplexSpecies([Sub, self.Enzyme])
        return [complex]

    def update_reactions(self, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        if kb == None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku == None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if kcat == None:
            kcat = component.get_parameter("kcat", part_id = component.name, mechanism = self)
        if component == None and (kb == None or ku == None or kcat == None):
            raise ValueError("Must pass in a Component or values for kb, ku, and kcat.")

        if complex == None:
            complex = ComplexSpecies([Sub, self.Enzyme])

        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction(inputs=[Sub, self.Enzyme], outputs=[complex],
                               k=kb, k_rev=ku)
        if Prod is not None:
            # Sub:Enz --> Enz + Prod
            cat_rxn = Reaction(inputs=[complex],
                               outputs=[Prod, self.Enzyme], k=kcat)
        else:  # Degradation Reaction
            # Sub:Enz --> Enz
            cat_rxn = Reaction(inputs=[complex], outputs=[self.Enzyme],
                               k=kcat)
        return [binding_rxn, cat_rxn]


class MichalisMentenCopyRXN(Mechanism):
    """In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Sub+Enz+Prod
    """
    def __init__(self, name, enzyme, mechanism_type = "copy", **keywords):
        if isinstance(enzyme, Species):
            self.Enzyme = enzyme
        else:
            raise ValueError("MichalisMentenCopyRXN takes a species object "
                             "for its enzyme argument")

        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Sub, **keywords):
        complex = ComplexSpecies([Sub, self.Enzyme])
        return [complex]

    def update_reactions(self, Sub, Prod, component = None, part_id = None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        if complex == None:
            complex = ComplexSpecies([Sub, self.Enzyme])

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
        binding_rxn = Reaction(inputs=[Sub, self.Enzyme], outputs=[complex],
                               k=kb, k_rev=ku)

        # Sub:Enz --> Enz + Prod + Sub
        cat_rxn = Reaction(inputs=[complex], outputs=[Sub, Prod, self.Enzyme],
                           k=kcat)

        return [binding_rxn, cat_rxn]