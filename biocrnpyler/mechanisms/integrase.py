# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Species

# from ..core.component import Component  # Commented out to avoid circular import


class BasicIntegration(Mechanism):
    """Mechanism for the schema DNA1 + DNA2 --> DNA3 + DNA4."""

    def __init__(
        self,
        name: str = "basic_integration",
        mechanism_type: str = "integration",
        **keywords,
    ):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs=None, **keywords):
        # this doesn't make any species because I use a Binding mechanism for that
        # maybe if we do the tetramerization mechanism then this would do something
        return []

    def update_reactions(
        self,
        DNA_inputs,
        DNA_outputs,
        component=None,
        part_id=None,
        kint=None,
        **keywords,
    ):
        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")
        elif kint is None:
            kint = component.get_parameter("kint", part_id=part_id, mechanism=self)

        return [
            Reaction.from_massaction(
                inputs=DNA_inputs, outputs=DNA_outputs, k_forward=kint
            )
        ]


class EnzymeIntegration(Mechanism):
    """Mechanism for the schema integrase+DNA1 + DNA2 --> integrase+DNA3 + DNA4."""

    def __init__(
        self,
        name: str = "enzyme_integration",
        mechanism_type: str = "integration",
        integrase="Int1",
        **keywords,
    ):
        """Initializes a BasicIntegration instance.

        :param name: name of the Mechanism, default: basic_integration
        :param mechanism_type: type of the Mechanism, default: integration
        :param keywords:
        """
        # TODO ZAT: remove unused keywords argument
        self.integrase = Species(name=integrase, material_type="protein")
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, DNA_inputs, DNA_outputs=None, **keywords):
        # this doesn't make any species because I use a Binding mechanism for that
        # maybe if we do the tetramerization mechanism then this would do something
        return []

    def update_reactions(
        self,
        DNA_inputs,
        DNA_outputs,
        component=None,
        part_id=None,
        kint=None,
        **keywords,
    ):

        if part_id is None and component is not None:
            part_id = component.name

        if kint is None and component is None:
            raise ValueError("Must pass in either a component or kint.")

        elif kint is None:
            kint = component.get_parameter("kint", part_id=part_id, mechanism=self)

        return [
            Reaction.from_massaction(
                inputs=[self.integrase] * 4 + DNA_inputs,
                outputs=[self.integrase] * 4 + DNA_outputs,
                k_forward=kint,
            )
        ]
