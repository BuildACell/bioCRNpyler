# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex


class BasicCatalysis(Mechanism):
    """Mechanism for the schema S + C --> P + C."""

    def __init__(self, name: str = "basic_catalysis", mechanism_type: str = "catalysis", **keywords):
        """Initializes a BasicCatalysis instance.

        :param name: name of the Mechanism, default: basic_catalysis
        :param mechanism_type: type of the Mechanism, default: catalysis
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod=None, **keywords):
        if Prod is None:
            return [Enzyme, Sub]
        else:
            return [Enzyme, Sub, Prod]

    def update_reactions(self, Enzyme, Sub, Prod, component=None, part_id=None, kcat=None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)

        return [Reaction.from_massaction(inputs=[Enzyme, Sub], outputs=[Enzyme, Prod], k_forward=kcat)]


class BasicProduction(Mechanism):
    """Mechanism for the schema C --> P + C."""

    def __init__(self, name="basic_production", mechanism_type="catalysis", **keywords):
        """Initializes a BasicProduction instance.

        :param name: name of the Mechanism, default: basic_production
        :param mechanism_type: type of the Mechanism, default: catalysis
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub=None, Prod=None, **keywords):
        species = [Enzyme]
        if Prod is not None:
            species += [Prod]
        if Sub is not None:
            species += [Sub]

        return species

    def update_reactions(self, Enzyme, Sub, Prod, component=None, part_id=None, kcat=None, **keywords):
        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)

        inputs = [Enzyme]
        outputs = [Enzyme]
        if Prod is not None:
            outputs += [Prod]
        if Sub is not None:
            inputs += [Sub]

        return [Reaction.from_massaction(inputs=inputs, outputs=outputs, k_forward=kcat)]


class MichaelisMenten(Mechanism):
    """Mechanism to automatically generate Michaelis-Menten Type Reactions.

       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz --> Enz+Prod
    """

    def __init__(self, name="michalis_menten", mechanism_type="catalysis", **keywords):
        """Initializes a MichaelisMenten instance.

        :param name: name of the Mechanism, default: michalis_menten
        :param mechanism_type: type of the Mechanism, default: catalysis
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod=None, complex=None, **keywords):
        if complex is None:
            complexS = Complex([Sub, Enzyme])
        else:
            complexS = complex
        if Prod is None:
            return [Enzyme, Sub, complexS]
        else:
            return [Enzyme, Sub, Prod, complexS]

    def update_reactions(self, Enzyme, Sub, Prod, component=None, part_id=None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        if ku is None:
            ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        if kcat is None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)

        if complex is None:
            complexS = Complex([Sub, Enzyme])
        else:
            complexS = complex

        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction.from_massaction(inputs=[Sub, Enzyme],
                                               outputs=[complexS],
                                               k_forward=kb,
                                               k_reverse=ku)
        if Prod is not None:
            # Sub:Enz --> Enz + Prod
            cat_rxn = Reaction.from_massaction(inputs=[complexS],
                                               outputs=[Prod, Enzyme],
                                               k_forward=kcat)
        else:  # Degradation Reaction
            # Sub:Enz --> Enz
            cat_rxn = Reaction.from_massaction(inputs=[complexS],
                                               outputs=[Enzyme],
                                               k_forward=kcat)
        return [binding_rxn, cat_rxn]


class MichaelisMentenReversible(Mechanism):
    """Mechanism to automatically generate Michaelis-Menten Type Reactions with products that can bind to enzymes.

       In the Copy RXN version, the Substrate is not Consumed
       Sub+Enz <--> Sub:Enz <--> Enz:Prod <--> Enz + Prod
    """

    def __init__(self, name="michalis_menten_reverse_binding", mechanism_type="catalysis", **keywords):
        """Initializes a MichaelisMentenReversible instance.

        :param name: name of the Mechanism, default: michalis_menten_reverse_binding
        :param mechanism_type: type of the Mechanism, default: catalysis
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, Prod, complex=None, complex2=None, **keywords):
        if complex is None:
            complex1 = Complex([Sub, Enzyme])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([Prod, Enzyme])
        else:
            complex2 = complex2
        return [Enzyme, Sub, Prod, complex1, complex2]

    def update_reactions(self, Enzyme, Sub, Prod, component=None, part_id=None, complex=None, complex2=None, kb=None, ku=None,
                         kcat=None, **keywords):
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter(
                "kb1", part_id=part_id, mechanism=self)
            kb2 = component.get_parameter(
                "kb2", part_id=part_id, mechanism=self)
        else:
            kb1, kb2 = kb
        if ku is None:
            ku1 = component.get_parameter(
                "ku1", part_id=part_id, mechanism=self)
            ku2 = component.get_parameter(
                "ku2", part_id=part_id, mechanism=self)
        else:
            ku1, ku2 = ku
        if kcat is None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)
            kcat_rev = component.get_parameter(
                "kcat_rev", part_id=part_id, mechanism=self)
        else:
            kcat, kcat_rev = kcat

        if complex is None:
            complex1 = Complex([Sub, Enzyme])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([Prod, Enzyme])

        # Sub + Enz <--> Sub:Enz
        binding_rxn1 = Reaction.from_massaction(inputs=[Sub, Enzyme],
                                                outputs=[complex1],
                                                k_forward=kb1,
                                                k_reverse=ku1)

        binding_rxn2 = Reaction.from_massaction(inputs=[Prod, Enzyme],
                                                outputs=[complex2],
                                                k_forward=kb2,
                                                k_reverse=ku2)

        # Sub:Enz --> Enz:Prod
        cat_rxn = Reaction.from_massaction(inputs=[complex1],
                                           outputs=[complex2],
                                           k_forward=kcat,
                                           k_reverse=kcat_rev)

        return [binding_rxn1, binding_rxn2, cat_rxn]


class MichaelisMentenCopy(Mechanism):
    """In the Copy RXN version, the Substrate is not Consumed.

       Sub+Enz <--> Sub:Enz --> Sub+Enz+Prod
    """

    def __init__(self, name="michalis_menten_copy", mechanism_type="copy", **keywords):
        """Initializes a MichaelisMentenCopy instance.

        :param name: name of the Mechanism, default: michalis_menten_copy
        :param mechanism_type: type of the Mechanism, default: copy
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, Enzyme, Sub, complex=None, Prod=None, **keywords):
        if complex is None:
            complexS = Complex([Sub, Enzyme])
        else:
            complexS = complex

        if Prod is None:
            return [Enzyme, Sub, complexS]
        elif (type(Prod) == list):
            return [Enzyme, Sub, complexS]+Prod
        else:
            return [Enzyme, Sub, Prod, complexS]

    def update_reactions(self, Enzyme, Sub, Prod, component=None, part_id=None, complex=None, kb=None, ku=None,
                         kcat=None, **keywords):
        if complex is None:
            complexS = Complex([Sub, Enzyme])
        else:
            complexS = complex

        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if kb is None and component is not None:
            kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        if ku is None and component is not None:
            ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        if kcat is None and component is not None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)
        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat.")
        # Sub + Enz <--> Sub:Enz
        binding_rxn = Reaction.from_massaction(inputs=[Sub, Enzyme],
                                               outputs=[complexS],
                                               k_forward=kb,
                                               k_reverse=ku)

        # Sub:Enz --> Enz + Prod + Sub
        cat_rxn = Reaction.from_massaction(inputs=[complexS],
                                           outputs=[Sub, Prod, Enzyme],
                                           k_forward=kcat)

        return [binding_rxn, cat_rxn]
