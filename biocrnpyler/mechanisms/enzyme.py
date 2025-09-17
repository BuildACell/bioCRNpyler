# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex


class BasicCatalysis(Mechanism):
    """Mechanism for the schema S + C --> P + C."""

    def __init__(
            self, name: str='basic_catalysis',
            mechanism_type: str='catalysis'):
        """Initializes a BasicCatalysis instance.

        :param name: name of the Mechanism, default: 'basic_catalysis'
        :param mechanism_type: type of the Mechanism, default: 'catalysis'
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, enzyme, substrate, product=None):
        if product is None:
            return [enzyme, substrate]
        else:
            return [enzyme, substrate, product]

    def update_reactions(
            self, enzyme, substrate, product, component=None,
            part_id=None, kcat=None):

        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self)

        if product is None:
            product = []

        return [Reaction.from_massaction(
            inputs=[enzyme, substrate], outputs=[enzyme, product],
            k_forward=kcat)]


class BasicProduction(Mechanism):
    """Mechanism for the schema C --> P + C."""

    def __init__(
            self, name='basic_production', mechanism_type='catalysis'):
        """Initializes a BasicProduction instance.

        :param name: name of the Mechanism, default: basic_production
        :param mechanism_type: type of the Mechanism, default: catalysis
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
            self, enzyme, substrate=None, product=None):

        species = [enzyme]
        if product is not None:
            species += [product]
        if substrate is not None:
            species += [substrate]

        return species

    def update_reactions(
            self, enzyme, substrate, product, component=None,
            part_id=None, kcat=None):

        if part_id is None and component is not None:
            part_id = component.name

        if kcat is None and component is None:
            raise ValueError("Must pass in either a component or kcat.")
        elif kcat is None:
            kcat = component.get_parameter(
                "kcat", part_id=part_id, mechanism=self)

        inputs = [enzyme]
        outputs = [enzyme]
        if product is not None:
            outputs += [product]
        if substrate is not None:
            inputs += [substrate]

        return [Reaction.from_massaction(
            inputs=inputs, outputs=outputs, k_forward=kcat)]


class MichaelisMenten(Mechanism):
    """Mechanism to automatically generate Michaelis-Menten Type Reactions.

       In the Copy RXN version, the substrates is not Consumed:
       sub + enz <--> sub:enz --> enz + prod
    """

    def __init__(
            self, name='michaelis_menten', mechanism_type='catalysis'):
        """Initializes a MichaelisMenten instance.

        :param name: name of the Mechanism, default: 'michaelis_menten'
        :param mechanism_type: type of the Mechanism, default: 'catalysis'
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
            self, enzyme, substrate, product=None, complex=None):

        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex
        if product is None:
            return [enzyme, substrate, complexS]
        else:
            return [enzyme, substrate, product, complexS]

    def update_reactions(
            self, enzyme, substrate, product, component=None, part_id=None,
            complex=None, kb=None, ku=None, kcat=None):

        # Get parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        if ku is None:
            ku = component.get_parameter('ku', part_id=part_id, mechanism=self)
        if kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self)

        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex

        # substrate + Enz <--> substrate:Enz
        binding_rxn = Reaction.from_massaction(
            inputs=[substrate, enzyme],
            outputs=[complexS], k_forward=kb, k_reverse=ku)
        if product is not None:
            # substrate:Enz --> Enz + product
            cat_rxn = Reaction.from_massaction(
                inputs=[complexS], outputs=[product, enzyme],
                k_forward=kcat)
        else:  # degredation Reaction
            # substrate:Enz --> Enz
            cat_rxn = Reaction.from_massaction(
                inputs=[complexS], outputs=[enzyme], k_forward=kcat)
        return [binding_rxn, cat_rxn]


class MichaelisMentenReversible(Mechanism):
    """Michaelis-Menten reactions with product that can bind to enzymes.

       In the Copy RXN version, the substrate is not Consumed
       sub + enz <--> sub:enz <--> enz:prod <--> enz + prod
    """

    def __init__(
            self, name='michaelis_menten_reverse_binding',
            mechanism_type='catalysis'):
        """Initializes a MichaelisMentenReversible instance.

        :param name: name of the Mechanism, default:
            'michaelis_menten_reverse_binding'
        :param mechanism_type: type of the Mechanism, default: 'catalysis'
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
            self, enzyme, substrate, product, complex=None, complex2=None):

        if complex is None:
            complex1 = Complex([substrate, enzyme])
        else:
            complex1 = complex
        if complex2 is None:
            complex2 = Complex([product, enzyme])
        else:
            complex2 = complex2
        return [enzyme, substrate, product, complex1, complex2]

    def update_reactions(
            self, enzyme, substrate, product, component=None,
            part_id=None, complex=None, complex2=None, kb=None, ku=None,
            kcat=None):

        # Get parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (kb is None or ku is None or kcat is None):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and kcat.")
        if kb is None:
            kb1 = component.get_parameter(
                'kb1', part_id=part_id, mechanism=self)
            kb2 = component.get_parameter(
                'kb2', part_id=part_id, mechanism=self)
        else:
            kb1, kb2 = kb
        if ku is None:
            ku1 = component.get_parameter(
                'ku1', part_id=part_id, mechanism=self)
            ku2 = component.get_parameter(
                'ku2', part_id=part_id, mechanism=self)
        else:
            ku1, ku2 = ku
        if kcat is None:
            kcat = component.get_parameter(
                'kcat', part_id=part_id, mechanism=self)
            kcat_rev = component.get_parameter(
                'kcat_rev', part_id=part_id, mechanism=self)
        else:
            kcat, kcat_rev = kcat

        if complex is None:
            complex1 = Complex([substrate, enzyme])
        else:
            complex1 = complex
        if complex2 == None:
            complex2 = Complex([product, enzyme])

        # substrate + Enz <--> substrate:Enz
        binding_rxn1 = Reaction.from_massaction(
            inputs=[substrate, enzyme], outputs=[complex1],
            k_forward=kb1, k_reverse=ku1)

        binding_rxn2 = Reaction.from_massaction(
            inputs=[product, enzyme], outputs=[complex2],
            k_forward=kb2, k_reverse=ku2)

        # substrate:Enz --> Enz:product
        cat_rxn = Reaction.from_massaction(
            inputs=[complex1], outputs=[complex2],
            k_forward=kcat, k_reverse=kcat_rev)

        return [binding_rxn1, binding_rxn2, cat_rxn]


class MichaelisMentenCopy(Mechanism):
    """In the Copy RXN version, the substrate is not consumed.

       substrate + Enz <--> substrate:Enz --> substrate + Enz + product
    """

    def __init__(
            self, name='michaelis_menten_copy', mechanism_type='copy'):
        """Initializes a MichaelisMentenCopy instance.

        :param name: name of the Mechanism, default: 'michaelis_menten_copy'
        :param mechanism_type: type of the Mechanism, default: 'copy'
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
            self, enzyme, substrate, complex=None, product=None):

        if complex is None:
            complexS = Complex([substrate, enzyme])
        else:
            complexS = complex

        if product is None:
            return [enzyme, substrate, complexS]
        elif (type(product) == list):
            return [enzyme, substrate, complexS]+product
        else:
            return [enzyme, substrate, product, complexS]

    def update_reactions(
            self, enzyme, substrate, product, component=None, part_id=None,
            complex=None, kb=None, ku=None, kcat=None):

        if complex is None:
            complexS = Complex([substrate, enzyme])
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
        # substrate + Enz <--> substrate:Enz
        binding_rxn = Reaction.from_massaction(
            inputs=[substrate, enzyme], outputs=[complexS],
            k_forward=kb, k_reverse=ku)

        # substrate:Enz --> Enz + product + substrate
        cat_rxn = Reaction.from_massaction(
            inputs=[complexS], outputs=[substrate, product, enzyme],
            k_forward=kcat)

        return [binding_rxn, cat_rxn]
