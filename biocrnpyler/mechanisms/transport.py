# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex
from ..core.propensities import ProportionalHillNegative, GeneralPropensity
from ..core.parameter import Parameter, ParameterEntry


class Simple_Diffusion(Mechanism):
    """A mechanism to model the diffusion of a substrate through a membrane channel.
    Does not require energy and follows diffusion rules.
    Reaction schema: substrate <-> product
    """

    def __init__(self, name="simple_diffusion", mechanism_type="diffusion", **keywords):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, substrate, product, **keywords):
        return [substrate, product]

    def update_reactions(
        self, substrate, product, component=None, part_id=None, k_diff=None, **keywords
    ):
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_diff is None):
            raise ValueError("Must pass in a Component or values for k_diff.")
        if k_diff is None:
            k_diff = component.get_parameter("k_diff", part_id=part_id, mechanism=self)
        else:
            k_diff = k_diff

        # Simple diffusion
        # Sub (Internal) <--> Product (External)
        diffusion_rxn = Reaction.from_massaction(
            inputs=[substrate], outputs=[product], k_forward=k_diff, k_reverse=k_diff
        )
        return [diffusion_rxn]


class Membrane_Protein_Integration(Mechanism):
    """A simple mehanism to integrate into the membrane protein in the membrane.
    Reaction schema for monomers: monomer -> intergral membrane protein
    Reaction schema for oligomer: monomer*[size] -> oligomer -> intergral membrane protein
    """

    def __init__(
        self,
        name="membrane_protein_integration",
        mechanism_type="membrane_insertion",
        **keywords,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self, integral_membrane_protein, product, complex=None, **keywords
    ):

        if complex is None:
            size = integral_membrane_protein.size
            if size > 1:
                complex1 = Complex([integral_membrane_protein] * size,
                                   compartment=integral_membrane_protein.compartment)
            else:
                complex1 = complex
        else:
            complex1 = complex

        return [integral_membrane_protein, product, complex1]

    def update_reactions(
        self,
        integral_membrane_protein,
        product,
        complex=None,
        component=None,
        part_id=None,
        **keywords,
    ):
        """This always requires the inputs component and part_id to find the relevant parameters"""

        # Get Parameters
        kb_oligomer = component.get_parameter(
            "kb_oligomer", part_id=part_id, mechanism=self
        )
        ku_oligomer = component.get_parameter(
            "ku_oligomer", part_id=part_id, mechanism=self
        )
        kex = component.get_parameter("kex", part_id=part_id, mechanism=self)
        kcat = component.get_parameter("kcat", part_id=part_id, mechanism=self)

        size = integral_membrane_protein.size

        if complex is None:
            if size > 1:
                complex1 = Complex([integral_membrane_protein] * size,
                                   compartment=integral_membrane_protein.compartment)
            else:
                complex1 = complex
        else:
            complex1 = complex

        # Membrane protein integration
        # Integration steps based on if protein is monomer or oligomer
        if size > 1:
            # homo: monomer --> oligomer
            binding_rxn1 = Reaction.from_massaction(
                inputs=[integral_membrane_protein] * size,
                outputs=[complex1],
                k_forward=kb_oligomer,
                k_reverse=ku_oligomer,
            )

            # oligomer-->integrated
            prophill_negative = ProportionalHillNegative(
                k=kex, d=complex1, K=kcat, n=4, s1=product
            )
            integration_rxn1 = Reaction(
                [complex1], [product], propensity_type=prophill_negative
            )
        else:
            # monomer-->integrated
            prophill_negative = ProportionalHillNegative(
                k=kex, d=integral_membrane_protein, K=kcat, n=4, s1=product
            )
            integration_rxn1 = Reaction(
                [integral_membrane_protein],
                [product],
                propensity_type=prophill_negative,
            )

        if size > 1:
            return [binding_rxn1, integration_rxn1]
        else:
            return [integration_rxn1]


class Simple_Transport(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane channel.
    Does not require energy and has unidirectional transport, following diffusion rules.
    Reaction schema: membrane_channel + substrate <-> membrane_channel + product
    """

    def __init__(
        self,
        name="simple_membrane_protein_transport",
        mechanism_type="transport",
        **keywords,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, membrane_channel, substrate, product, **keywords):
        if membrane_channel.attributes[0] != "Passive":
            raise ValueError(
                "Protein is not classified as a channel with passive transport of small molecules. Use mechanism Facilitated_Passive_Transport instead."
            )

        return [membrane_channel, substrate, product]

    def update_reactions(
        self,
        membrane_channel,
        substrate,
        product,
        component=None,
        part_id=None,
        k_trnsp=None,
        **keywords,
    ):
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_trnsp is None):
            raise ValueError("Must pass in a Component or values for k_trnsp.")
        if k_trnsp is None:
            k_trnsp = component.get_parameter(
                "k_trnsp", part_id=part_id, mechanism=self
            )
        else:
            k_trnsp = k_trnsp

        # Simple membrane protein transport
        # Sub (Internal) <--> Product (External)
        SimpleTransport_rxn = Reaction.from_massaction(
            inputs=[substrate, membrane_channel],
            outputs=[product, membrane_channel],
            k_forward=k_trnsp,
            k_reverse=k_trnsp,
        )
        return [SimpleTransport_rxn]


class Facilitated_Transport_MM(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane carrier.
    Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.
    Mechanism for the schema: Sub+MC <--> Sub:MC --> Prod:MC --> Prod + MC
    """

    def __init__(
        self,
        name="facilitated_membrane_protein_transport",
        mechanism_type="transport",
        **keywords,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self, membrane_carrier, substrate, product, complex_dict=None, **keywords
    ):

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict["sub:MC"] = Complex([substrate, membrane_carrier],
                                             compartment=substrate.compartment)
            # Complex2
            complex_dict["prod:MC"] = Complex([product, membrane_carrier],
                                              compartment=product.compartment)

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [membrane_carrier, substrate, product, complex_array]

    def update_reactions(
        self,
        membrane_carrier,
        substrate,
        product,
        complex_dict=None,
        component=None,
        part_id=None,
        **keywords,
    ):
        """This always requires the inputs component and part_id to find the relevant parameters"""

        # Get Parameters
        kb_subMC = component.get_parameter("kb_subMC", part_id=part_id, mechanism=self)
        ku_subMC = component.get_parameter("ku_subMC", part_id=part_id, mechanism=self)
        k_trnspMC = component.get_parameter(
            "k_trnspMC", part_id=part_id, mechanism=self
        )
        ku_prodMC = component.get_parameter(
            "ku_prodMC", part_id=part_id, mechanism=self
        )

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict["sub:MC"] = Complex([substrate, membrane_carrier],
                                             compartment=substrate.compartment)
            # Complex2
            complex_dict["prod:MC"] = Complex([product, membrane_carrier],
                                              compartment=product.compartment)

        # Facilitated membrane protein transport
        # Sub + MC --> Sub:MC
        prop_subMC = GeneralPropensity(
            f"kb_subMC*{substrate}*{membrane_carrier}*Heaviside({substrate}-{product})-kb_subMC*{product}*{membrane_carrier}*Heaviside({substrate}-{product})",
            propensity_species=[product, substrate, membrane_carrier],
            propensity_parameters=[kb_subMC],
        )
        binding_rxn1 = Reaction(
            [substrate, membrane_carrier],
            [complex_dict["sub:MC"]],
            propensity_type=prop_subMC,
        )

        # Sub:MC --> Sub + MC
        unbinding_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict["sub:MC"]],
            outputs=[membrane_carrier, substrate],
            k_forward=ku_subMC,
        )

        # Sub:MC --> Prod:MC
        transport_rxn = Reaction.from_massaction(
            inputs=[complex_dict["sub:MC"]],
            outputs=[complex_dict["prod:MC"]],
            k_forward=k_trnspMC,
        )

        # MC:Prod --> MC + Prod
        unbinding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict["prod:MC"]],
            outputs=[product, membrane_carrier],
            k_forward=ku_prodMC,
        )

        return [binding_rxn1, unbinding_rxn1, transport_rxn, unbinding_rxn2]


class Primary_Active_Transport_MM(Mechanism):
    """A mechanism to model the transport of a substrate through a membrane carrier.
    Mechanism follows Michaelis-Menten Type Reactions with products that can bind to membrane carriers.
    Mechanism for the schema: Sub+MP <--> Sub:MP + E --> Sub:MP:E --> MP:Prod:E --> Prod + MP:W --> Prod + MP+ W
    """

    def __init__(
        self,
        name="active_membrane_protein_transport",
        mechanism_type="transport",
        **keywords,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self,
        membrane_pump,
        substrate,
        product,
        energy,
        waste,
        complex_dict=None,
        **keywords,
    ):

        nATP = membrane_pump.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict["Pump:Sub"] = Complex([substrate, membrane_pump],
                                               compartment=substrate.compartment)
            # Complex2
            complex_dict["Pump:Sub:ATP"] = Complex([nATP * [energy],
                                                    complex_dict["Pump:Sub"]],
                                                   compartment=substrate.compartment)
            # Complex3
            complex_dict["Pump:Prod:ATP"] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=product.compartment
            )
            # Complex4
            complex_dict["Pump:ADP"] = Complex([nATP * [waste], membrane_pump],
                                               compartment=membrane_pump.compartment)

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [membrane_pump, substrate, product, energy, waste, complex_array]

    def update_reactions(
        self,
        membrane_pump,
        substrate,
        product,
        energy,
        waste,
        complex_dict=None,
        component=None,
        part_id=None,
        **keywords,
    ):
        """This always requires the inputs component and part_id to find the relevant parameters"""

        # Get Parameters
        kb_subMP = component.get_parameter("kb_subMP", part_id=part_id, mechanism=self)
        ku_subMP = component.get_parameter("ku_subMP", part_id=part_id, mechanism=self)
        kb_subMPnATP = component.get_parameter(
            "kb_subMPnATP", part_id=part_id, mechanism=self
        )
        ku_subMPnATP = component.get_parameter(
            "ku_subMPnATP", part_id=part_id, mechanism=self
        )
        k_trnspMP = component.get_parameter(
            "k_trnspMP", part_id=part_id, mechanism=self
        )
        ku_prodMP = component.get_parameter(
            "ku_prodMP", part_id=part_id, mechanism=self
        )
        ku_MP = component.get_parameter("ku_MP", part_id=part_id, mechanism=self)

        nATP = membrane_pump.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}

            # Complex1
            complex_dict["Pump:Sub"] = Complex([substrate, membrane_pump],
                                               compartment=substrate.compartment)
            complex1 = complex_dict["Pump:Sub"]

            # Complex2
            complex_dict["Pump:Sub:ATP"] = Complex(
                [nATP * [energy], complex_dict["Pump:Sub"]],
                compartment=substrate.compartment,
            )
            # Complex3
            complex_dict["Pump:Prod:ATP"] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=product.compartment,
            )

            # Complex4
            complex_dict["Pump:ADP"] = Complex([nATP * [waste], membrane_pump],
                                               compartment=membrane_pump.compartment)

        # Active membrane protein transport
        # Sub + MP<--> Sub:MP
        prop_subMP = GeneralPropensity(
            f"kb_subMP*{substrate}*{membrane_pump}*Heaviside({membrane_pump})",
            propensity_species=[substrate, membrane_pump],
            propensity_parameters=[kb_subMP],
        )
        binding_rxn1 = Reaction(
            [substrate, membrane_pump],
            [complex_dict["Pump:Sub"]],
            propensity_type=prop_subMP,
        )

        unbinding_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict["Pump:Sub"]],
            outputs=[substrate, membrane_pump],
            k_forward=ku_subMP,
        )

        # Sub:MP + E <--> Sub:MP:E

        prop_subMPnATP = GeneralPropensity(
            f"kb_subMPnATP*{complex1}*{energy}*Heaviside({complex1})",
            propensity_species=[complex1, energy],
            propensity_parameters=[kb_subMPnATP],
        )
        binding_rxn2 = Reaction(
            [complex1, nATP * [energy]],
            [complex_dict["Pump:Sub:ATP"]],
            propensity_type=prop_subMPnATP,
        )

        unbinding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict["Pump:Sub:ATP"]],
            outputs=[complex_dict["Pump:Sub"], nATP * [energy]],
            k_forward=ku_subMPnATP,
        )

        # Sub:MP:E --> Prod:MP:E
        transport_rxn = Reaction.from_massaction(
            inputs=[complex_dict["Pump:Sub:ATP"]],
            outputs=[complex_dict["Pump:Prod:ATP"]],
            k_forward=k_trnspMP,
        )
        # Prod:MP:E--> Prod+MP:W
        unbinding_rxn3 = Reaction.from_massaction(
            inputs=[complex_dict["Pump:Prod:ATP"]],
            outputs=[complex_dict["Pump:ADP"], product],
            k_forward=ku_prodMP,
        )
        # MP:W --> MP+W
        unbinding_rxn4 = Reaction.from_massaction(
            inputs=[complex_dict["Pump:ADP"]],
            outputs=[nATP * [waste], membrane_pump],
            k_forward=ku_MP,
        )

        return [
            binding_rxn1,
            unbinding_rxn1,
            binding_rxn2,
            unbinding_rxn2,
            transport_rxn,
            unbinding_rxn3,
            unbinding_rxn4,
        ]
