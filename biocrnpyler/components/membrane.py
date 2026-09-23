#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import warnings
from typing import List, Union

from ..core.compartment import Compartment
from ..core.component import Component
from ..mechanisms import (
    Diffusion_Simple, Integration_MembraneProtein, Diffusion_Facilitated_Channel,
    Diffusion_Facilitated_Carrier, Transport_PrimaryActive_ABCexporter,
    Sensor_TwoComponentSystem)
from ..core.species import Species


class DiffusibleMolecule(Component):
    r"""Molecule that diffuses passively through a membrane.

    A `DiffusibleMolecule` component represents a molecule that is membrane
    permeable, thus can diffuse across a membrane separating two compartments.
    The component uses the 'diffusion' mechanism to generate bidirectional
    diffusion reactions based on concentration gradients.

    Parameters
    ----------
    substrate : Species, str, or Component
        The diffusible molecule species. Can be a or a List of `Species`
        object, string name, or `Component` with an associated species.
    internal_compartment : str or Compartment, default='internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='external'
        The external compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    attributes : list of str, optional
        List of attribute tags to associate with the substrate species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    substrate : Species
        The substrate species in the internal compartment.
    product : Species
        The same substrate species in the external compartment (diffusion
        product).

    See Also
    --------
    MembraneChannel : Membrane channel for substrate diffusion.
    MembraneCarrier : Membrane carrier for substrate diffusion or transport.
    MembranePump : Membrane pump for ATP-dependent substrate transport.
    Component : Base class for biomolecular components.

    Notes
    -----
    Passive diffusion follows concentration gradients and does not require
    energy. The diffusion mechanism generates bidirectional reactions:

    - Forward: substrate_internal --> substrate_external

    - Reverse: substrate_external --> substrate_internal

    If not specified using the `name` keyword, the component name is
    automatically generated as: '<substrate_name>_<internal_compartment_name>'

    Examples
    --------
    Create a simple diffusible molecule:

    >>> glucose = bcp.DiffusibleMolecule(
    ...     substrate='Glucose',
    ...     internal_compartment='cytoplasm',
    ...     external_compartment='extracellular'
    ... )

    Use with a mixture and diffusion mechanism:

    >>> mixture = bcp.Mixture(
    ...     components=[glucose],
    ...     mechanisms={'diffusion': bcp.Diffusion_Simple()},
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        substrate: Union[List[Union[Species, str, Component]],
                         Union[Species, str, Component]],
        internal_compartment: Union[str, Compartment] = 'internal',
        external_compartment: Union[str, Compartment] = 'external',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Catch single item
        if not isinstance(substrate, list):
            substrate = [substrate]

        # Populate the internal and external species lists
        self.substrate = []
        self.product = []

        for sub in substrate:
            if isinstance(sub, Species):
                if sub.compartment.name != internal_compartment:
                    sub.compartment = internal_compartment
                sub_species = sub
            elif isinstance(sub, str):
                sub_species = self.set_species(
                    sub,
                    compartment=internal_compartment,
                    attributes=attributes,
                )

            prod_species = self.set_species(
                sub_species.name,
                material_type=sub_species.material_type,
                compartment=external_compartment,
                attributes=sub_species.attributes,
            )

            self.substrate.append(sub_species)
            self.product.append(prod_species)

        # Name the component dynamically
        if (name := kwargs.pop('name', None)) is None:
            if len(self.substrate) == 1:
                name = (
                    self.substrate[0].name
                    + '_'
                    + self.substrate[0].compartment.name)
            else:
                name = f"MultiDiffusible_{internal_compartment.name}"

        Component.__init__(
            self=self, name=name, attributes=attributes,
            mechanisms={'diffusion':Diffusion_Simple()}, **kwargs
        )

    def get_species(self):
        """Get the substrate species in the internal compartment.

        Returns
        -------
        Species
            The substrate species in the internal compartment.

        """
        return self.substrate

    def update_species(self):
        """Uses the 'diffusion' mechanism to generate diffusion species.

        Uses the 'diffusion' mechanism to generate species in both
        compartments.

        Returns
        -------
        list of Species
            List of species in internal and external compartments generated
            by the diffusion mechanism.

        """
        mech_diff = self.get_mechanism('diffusion')
        species_list = []

        for sub, prod in zip(self.substrate, self.product):
            species_list.extend(mech_diff.update_species(sub, prod))

        # Return a unique list of species
        return list(set(species_list))

    def update_reactions(self):
        """Use 'diffusion' mechanism to generate diffusion reactions.

        Uses the 'diffusion' mechanism to generate reactions for
        diffusion between compartments.

        Returns
        -------
        list of Reaction
            List of diffusion reactions (forward and reverse) between
            internal and external compartments.

        """
        mech_diff = self.get_mechanism('diffusion')
        reactions_list = []

        for sub, prod in zip(self.substrate, self.product):
            reactions_list.extend(
                mech_diff.update_reactions(
                    sub, prod,
                    component=self,
                    part_id=sub.name)
                )

        return reactions_list


class IntegralMembraneProtein(Component):
    """A protein that integrates into the membrane.

    An `IntegralMembraneProtein` component represents a protein
    that integrates into a membrane compartment. The component uses the
    'membrane_integration' mechanism to generate reactions for protein
    integration into the membrane. The size parameter allows modeling of
    oligomeric channels (dimers, trimers, etc.).

    Parameters
    ----------
    membrane_protein : Species, str, or Component
        The membrane protein species before insertion. Can be a `Species`
        object, string name, or `Component` with an associated species.
    product : Species, str, or Component
        The integrated membrane protein species. Can be a `Species` object,
        string name, or `Component`.
    size : int, optional
        Number of monomers needed to form the functional channel. Used to
        model oligomeric channels (e.g., size=2 for dimers, size=3 for
        trimers). Default is 1.
    compartment : str or Compartment, default='internal'
        The compartment containing the membrane protein before insertion.
        Can be a string name or `Compartment` object.
    membrane_compartment : str or Compartment, default='membrane'
        The membrane compartment where the protein integrates. Can be a
        string name or `Compartment` object.
    attributes : list of str, optional
        List of attribute tags to associate with the membrane protein.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    membrane_protein : Species
        The membrane protein species before insertion.
    product : Species
        The integrated transmembrane protein species in the membrane
        compartment.

    See Also
    --------
    MembraneChannel : Membrane channel for substrate diffusion.
    MembraneCarrier : Membrane carrier for substrate diffusion or transport.
    MembranePump : Membrane pump for ATP-dependent substrate transport.
    Component : Base class for biomolecular components.

    Notes
    -----
    The membrane_integration mechanism generates reactions for protein
    integration into the membrane. For oligomeric channels, the size
    parameter determines the stoichiometry:

    - size=1: Monomer insertion
    - size=2: Dimer formation (2 proteins --> 1 channel)
    - size=3: Trimer formation (3 proteins --> 1 channel)

    The component name is automatically generated as:
    '<membrane_protein_name>_<compartment_name>'

    Examples
    --------
    Create a simple membrane protein:

    >>> channel = bcp.IntegralMembraneProtein(
    ...     membrane_protein='ChannelProtein',
    ...     product='ChannelProtein_membrane',
    ... )

    Create a dimeric channel protein:

    >>> dimer = bcp.IntegralMembraneProtein(
    ...     membrane_protein='Aquaporin',
    ...     product='Aquaporin_channel',
    ...     size=2,
    ... )

    """

    def __init__(
        self,
        membrane_protein: Union[Species, str],
        product: Union[Species, str, Component],
        size: int = None,
        compartment: Union[str, Compartment] = 'internal',
        membrane_compartment: Union[str, Compartment] = 'membrane',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(compartment, str):
            compartment = Compartment(name=compartment)
        if isinstance(membrane_compartment, str):
            membrane_compartment = Compartment(name=membrane_compartment)

        # PROTEIN
        if isinstance(membrane_protein, Species):
            self.membrane_protein = membrane_protein
        elif isinstance(membrane_protein, str):
            membrane_protein_name = membrane_protein
            self.membrane_protein = self.set_species(
                membrane_protein_name,
                material_type='protein',
                compartment=compartment,
                attributes=attributes,
            )
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(membrane_protein)}")

        # Logic for prioritizing compartments
        if self.membrane_protein.compartment.name == 'default':
            self.membrane_protein.compartment = compartment
        elif (
            self.membrane_protein.compartment.name != compartment.name
            and compartment.name == 'internal'
        ):
            warnings.warn(
                "Inconsistent compartments, prioritizing membrane protein "
                "compartment.",
                UserWarning,
            )
            compartment = self.membrane_protein.compartment
        elif (self.membrane_protein.compartment.name != compartment.name
            and compartment.name != 'internal'
        ):
            warnings.warn(
                "Inconsistent compartments, prioritizing integral membrane "
                "protein compartment.",
                UserWarning,
            )
            self.membrane_protein.compartment = compartment

        # PRODUCT is an integrated membrane protein
        if product is None:
            product_name = self.membrane_protein.name + '_IMP'
        elif isinstance(product, Component):
            product_name = product.get_species()[0].name + '_IMP'
        elif isinstance(product, Species):
            product_name = product.name
        elif isinstance(product, str):
            product_name = product
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(product)}")

        self.product = self.set_species(
            product_name,
            material_type='protein',
            compartment=membrane_compartment,
        )

        # Indicates the number of monomers that compose the membrane protein,
        # will be used in Integration_MembraneProtein(Mechanism)
        if size is None:
            self.membrane_protein.size = 1
        else:
            self.membrane_protein.size = size

        # Name the component
        name = (
            self.membrane_protein.name
            + '_'
            + self.membrane_protein.compartment.name
        )

        Component.__init__(
            self=self,
            name=name,
            mechanisms={
                'membrane_integration': Integration_MembraneProtein()
            },
            **kwargs)

    def get_species(self):
        """Get the membrane protein species before integration.

        Returns
        -------
        Species
            The membrane protein species in the compartment before
            integration into the membrane.

        """
        return self.membrane_protein

    def update_species(self):
        """Uses the 'membrane_integration' to generate integration species.

        Uses the 'membrane_integration' mechanism to generate species for
        the protein before and after integration.

        Returns
        -------
        list of Species
            List of species generated by the membrane_integration mechanism,
            including the protein and integrated product.

        """
        mech_ins = self.get_mechanism('membrane_integration')
        return mech_ins.update_species(self.membrane_protein, self.product)

    def update_reactions(self):
        """Use 'membrane_integration' to generate integration reactions.

        Uses the 'membrane_integration' mechanism to generate reactions for
        protein integration into the membrane.

        Returns
        -------
        list of Reaction
            List of reactions for protein integration into the membrane.

        """
        mech_ins = self.get_mechanism('membrane_integration')
        return mech_ins.update_reactions(
            self.membrane_protein,
            self.product,
            component=self,
            part_id=self.name,
        )


class MembraneChannel(Component):
    """Membrane proteins that facilitate solute diffusion across membranes.

    A `MembraneChannel` component represents a membrane channel or
    pore that facilitates substrate movement across a membrane
    following concentration gradients. The component uses the 'diffusion'
    mechanism to generate diffusion reactions.

    Parameters
    ----------
    membrane_channel : Species, str, or Component
        The integral membrane protein that forms the channel. Can be a
        `Species` object, string name, or `Component`. If a string,
        automatically creates a protein species with appropriate direction
        attribute.
    substrate : Species, str, or Component
        The substrate to be diffused through the channel. Can be a
        `Species` object, string name, or `Component`.
    internal_compartment : str or Compartment, default='internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='external'
        The external compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    attributes : list of str, optional
        List of attribute tags to associate with substrate species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    membrane_channel : Species
        The membrane channel protein species.
    substrate : Species
        The substrate species in the source compartment (depends on
        direction).
    product : Species
        The same substrate in the destination compartment.

    See Also
    --------
    IntegralMembraneProtein : Protein insertion into membranes.
    MembranePump : Membrane pump for ATP-dependent substrate transport.
    DiffusibleMolecule : Passive diffusion without channels.
    Component : Base class for biomolecular components.

    Notes
    -----
    The component name is automatically generated as:
    '<membrane_channel_name>_<compartment_name>'

    Examples
    --------
    Create a passive channel:

    >>> channel = bcp.MembraneChannel(
    ...     membrane_channel='WaterChannel',
    ...     substrate='Water',
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[channel],
    ...     mechanisms={'diffusion': bcp.Diffusion_Facilitated_Channel()},
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_channel: Union[Species, str, Component],
        substrate: Union[List[Union[Species, str, Component]],
                         Union[Species, str, Component]],
        internal_compartment: Union[str, Compartment] = 'internal',
        external_compartment: Union[str, Compartment] = 'external',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Channel
        if isinstance(membrane_channel, Component):
            self.membrane_channel = membrane_channel.product
        elif isinstance(membrane_channel, Species):
            self.membrane_channel = membrane_channel
        elif isinstance(membrane_channel, str):
            self.membrane_channel = self.set_species(
                membrane_channel,
                material_type='protein',
                attributes=attributes,
            )
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(membrane_channel)}")

        # Substrate and product assignments.
        # In the case of membrane components, the substrate is the
        # substance on which the channel acts without distinction
        # of compartment. The substrate and product are the same substance
        # and the substance does not change as a result except for the
        # compartment. The substrate and the product are explicitly given
        # the `internal_compartment` and `external_compartment`, respectively.

        if substrate is None:
            substrate = []
        elif not isinstance(substrate, list):
            substrate = [substrate]

        # Initialize lists for substrate and products
        self.substrate_in = []
        self.substrate_out = []

        # Iterate over each substrate
        for sub in substrate:
            if isinstance(sub, Species):
                sub_species = sub
                if sub.compartment.name != internal_compartment.name:
                    sub.compartment = internal_compartment
            elif isinstance(sub, Component):
                sub_name = sub.get_species()[0].name
                sub_species = self.set_species(
                    sub_name,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
            elif isinstance(sub, str):
                sub_name = sub
                sub_species = self.set_species(
                    sub,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
            else:
                raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(sub)}")

            prod_species = self.set_species(
                sub_species.name,
                material_type=sub_species.material_type,
                compartment=external_compartment,
                attributes=sub_species.attributes,
            )

            self.substrate_in.append(sub_species)
            self.substrate_out.append(prod_species)

        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = (
                self.membrane_channel.name
                + '_'
                + self.membrane_channel.compartment.name
            )

        Component.__init__(
            self=self, name=name,
            mechanisms={'diffusion':Diffusion_Facilitated_Channel()},
            **kwargs)

    def get_species(self):
        """Get the membrane channel species.

        Returns
        -------
        Species
            The integral membrane protein species that
              forms the channel or pore.

        """
        return self.membrane_channel

    def update_species(self):
        """Uses the 'diffusion' mechanism to generate channel species.

        Uses the 'diffusion' mechanism to generate species including the
        channel protein, substrate, and product.

        Returns
        -------
        list of Species
            List of species generated by the diffusion mechanism.

        """
        mech_tra = self.get_mechanism('diffusion')
        species_list = []

        for sub_in, sub_out in zip(self.substrate_in, self.substrate_out):
            species_list.extend(
                mech_tra.update_species(
                    self.membrane_channel, sub_in, sub_out
                )
            )

        return list(set(species_list))

    def update_reactions(self):
        """Use 'diffusion' mechanism to generate channel-mediated reactions.

        Uses the 'diffusion' mechanism to generate reactions for substrate
        diffusion through the channel.

        Returns
        -------
        list of Reaction
            List of diffusion reactions through the membrane channel.

        """
        mech_tra = self.get_mechanism('diffusion')
        reactions_list = []

        for sub_in, sub_out in zip(self.substrate_in, self.substrate_out):
            reactions_list.extend(
                mech_tra.update_reactions(
                    self.membrane_channel,
                    sub_in,
                    sub_out,
                    component=self,
                    part_id=f"{self.name}_{sub_in.name}"
                )
            )
        return reactions_list


class MembraneCarrier(Component):
    """Membrane proteins that facilitate solute diffusion or transport across.

    A `MembraneCarrier` component represents a membrane carrier facilitates
    substrate movement across a membrane following concentration gradients of
    the substrate or driving ion. The component uses the 'diffusion' or
    'transport' mechanism to generate diffusion or transport reactions.

    Parameters
    ----------
    membrane_carrier : Species, str, or Component
        The integral membrane protein that forms the carrier. Can be a
        `Species` object, string name, or `Component`. If a string,
        automatically creates a protein species with appropriate direction
        attribute.
    substrate : Species, str, or Component
        The substrate to be transported through the channel. Can be a
        `Species` object, string name, or `Component`.
    internal_compartment : str or Compartment, default='internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='external'
        The external compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    attributes : list of str, optional
        List of attribute tags to associate with substrate species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    membrane_carrier : Species
        The membrane channel protein species.
    substrate : Species
        The substrate species in the source compartment (depends on
        direction).
    product : Species
        The same substrate in the destination compartment.

    See Also
    --------
    IntegralMembraneProtein : Protein integration into membranes.
    MembranePump : Membrane pump for ATP-dependent substrate transport.
    DiffusibleMolecule : Passive diffusion without channels.
    Component : Base class for biomolecular components.

    Notes
    -----
    The component name is automatically generated as:
    '<membrane_carrier_name>_<compartment_name>'

    Examples
    --------
    Create a passive carrier:

    >>> carrier = bcp.MembraneCarrier(
    ...     membrane_carrier='SubCarrier',
    ...     substrate='S1',
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[carrier],
    ...     mechanisms={'diffusion': bcp.Diffusion_Facilitated_Carrier()},
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_carrier: Union[Species, str, Component],
        substrate: Union[List[Union[Species, str, Component]],
                         Union[Species, str, Component]],
        internal_compartment: Union[str, Compartment] = 'internal',
        external_compartment: Union[str, Compartment] = 'external',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Carrier
        if isinstance(membrane_carrier, Component):
            self.membrane_carrier = membrane_carrier.product
        elif isinstance(membrane_carrier, Species):
            self.membrane_carrier = membrane_carrier
        elif isinstance(membrane_carrier, str):
            self.membrane_carrier = self.set_species(
                membrane_carrier,
                material_type='protein',
                attributes=attributes,
            )
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(membrane_carrier)}")

        # Substrate and product assignments.
        #
        # In the case of membrane components, the substrate is the
        # substance on which the carrier acts without distinction
        # of compartment. The substrate and product are the same substance
        # and the substance does not change as a result except for the
        # compartment. The substrate and the product are explicitly given the
        # `internal_compartment` and `external_compartment`, respectively.
        if substrate is None:
                    substrate = []
        elif not isinstance(substrate, list):
            substrate = [substrate]

        self.substrate_in = []
        self.substrate_out = []
        # Substrate
        # Iterate over each substrate
        for sub in substrate:
            if isinstance(sub, Species):
                sub_species = sub

                if sub.compartment.name != internal_compartment:
                    sub.compartment = internal_compartment
            elif isinstance(sub, Component):
                sub_name = sub.get_species()[0].name

                sub_species = self.set_species(
                    sub_name,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
            elif isinstance(sub, str):
                sub_name = sub
                sub_species = self.set_species(
                    sub,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
            else:
                raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(sub)}")

            prod_species = self.set_species(
                sub_species.name,
                material_type=sub_species.material_type,
                compartment=external_compartment,
                attributes=sub_species.attributes,
            )

            self.substrate_in.append(sub_species)
            self.substrate_out.append(prod_species)

        # Name the component
        name = (
            self.membrane_carrier.name
            + '_'
            + self.membrane_carrier.compartment.name
        )

        Component.__init__(
            self=self, name=name,
            mechanisms={'diffusion':Diffusion_Facilitated_Carrier()},
            **kwargs)

    def get_species(self):
        """Get the membrane carrier species.

        Returns
        -------
        Species
            The integral membrane protein species that
              forms the carrier.

        """
        return self.membrane_carrier

    def update_species(self):
        """Uses the 'diffusion' or 'transport' mechanism to generate species.

        Uses the 'diffusion' or 'transport' mechanism to generate species
        including the carrier protein, substrate, and product.

        Returns
        -------
        list of Species
            List of species generated by the diffusion or transport mechanism.

        """
        mech_tra = self.get_mechanism('transport', optional_mechanism=True)

        if mech_tra is None:
            mech_tra = self.get_mechanism('diffusion', optional_mechanism=True)

        if mech_tra is None:
            raise KeyError(
                f"Unable to find mechanism of type diffusion or transport in "
                f"Component {self}."
            )

        species_list = []
        for sub_in, sub_out in zip(self.substrate_in, self.substrate_out):
            try:
                species_list.extend(
                    mech_tra.update_species(
                        self.membrane_carrier,
                        sub_in,
                        sub_out,
                        driving_ion=None
                    )
                )
            except TypeError:
                species_list.extend(
                    mech_tra.update_species(
                        self.membrane_carrier,
                        sub_in,
                        sub_out
                    )
                )

        return list(set(species_list))

    def update_reactions(self):
        """Use 'diffusion' or 'transport' mechanism to generate reactions.

        Uses the 'diffusion' or 'transport' mechanism to generate reactions
        for substrate diffusion or transport through the carrier.

        Returns
        -------
        list of Reaction
            List of diffusion reactions through the membrane channel.

        """
        mech_tra = self.get_mechanism('transport', optional_mechanism=True)

        if mech_tra is None:
            mech_tra = self.get_mechanism('diffusion', optional_mechanism=True)

        if mech_tra is None:
            raise KeyError(
                f"Unable to find mechanism of type diffusion or transport in "
                f"Component {self}."
            )

        if mech_tra is None:
            raise KeyError(
                f"Unable to find mechanism of type diffusion or transport in "
                f"Component {self}."
            )

        reactions_list = []
        for sub_in, sub_out in zip(self.substrate_in, self.substrate_out):
            try:
                reactions_list.extend(
                    mech_tra.update_reactions(
                        self.membrane_carrier,
                        sub_in,
                        sub_out,
                        driving_ion=None,
                        component=self,
                        part_id=self.name,
                    )
                )
            except TypeError:
                reactions_list.extend(
                    mech_tra.update_reactions(
                        self.membrane_carrier,
                        sub_in,
                        sub_out,
                        component=self,
                        part_id=f"{self.name}_{sub_in.name}"
                    )
                )

        return reactions_list


class MembranePump(Component):
    """An ATP-dependent membrane protein that enables active transport.

    A `MembranePump` component represents an active transporter or pump
    that uses ATP to transport substrates across membranes against
    concentration gradients. The pump operates unidirectionally and requires
    energy in the form of ATP. The component uses the 'transport' mechanism
    to generate ATP-dependent transport reactions.

    Parameters
    ----------
    membrane_pump : Species, str, or Component
        The membrane pump protein species. Can be a `Species` object,
        string name, or `Component`. If a string, automatically creates a
        protein species with appropriate direction attribute.
    substrate : Species, str, or Component
        The substrate to be transported by the pump. Can be a `Species`
        object, string name, or `Component`.
    direction : str, optional
        Direction of active transport. Common values: 'importer'
        (external --> internal), 'exporter' (internal --> external),
        Affects substrate and ATP compartment placement.
    internal_compartment : str or Compartment, default='internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='external'
        The external compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    ATP : int, optional
        Number of ATP molecules required per transport cycle. Default is 1.
    attributes : list of str, optional
        List of attribute tags to associate with substrate species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    membrane_pump : Species
        The membrane pump protein species.
    substrate : Species
        The substrate species in the source compartment.
    product : Species
        The same substrate in the destination compartment.
    energy : Species
        ATP species used for energy (compartment depends on direction).
    waste : Species
        ADP species produced (compartment depends on direction).

    See Also
    --------
    MembraneChannel : Membrane channel for substrate diffusion.
    MembraneCarrier : Membrane carrier for substrate diffusion or transport.
    DiffusibleMolecule : Passive diffusion without channels.
    Component : Base class for biomolecular components.

    Notes
    -----
    Active transport requires ATP hydrolysis and can move substrates
    against concentration gradients. The typical reaction scheme is:

    - Exporter: substrate_internal + ATP + pump -->
      substrate_external + ADP + pump

    - Importer: substrate_external + ATP + pump -->
      substrate_internal + ADP + pump

    The ATP parameter controls the stoichiometry of ATP consumption per
    transport event.

    The component name is automatically generated as:
    '<membrane_pump_name>_<compartment_name>'

    Examples
    --------
    Create a simple ATP-dependent exporter:

    >>> pump = bcp.MembranePump(
    ...     membrane_pump='CalciumPump',
    ...     substrate='Calcium',
    ...     direction='importer',
    ...     ATP=2
    ... )

    Create an ABC transporter (importer):

    >>> abc = bcp.MembranePump(
    ...     membrane_pump='ABC_Transporter',
    ...     substrate='Maltose',
    ...     ATP=1
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[pump],
    ...     mechanisms={
    ...         'transport': bcp.Transport_PrimaryActive_ABCexporter()
    ...     },
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_pump: Union[Species, str, Component],
        substrate: Union[List[Union[Species, str, Component]],
                         Union[Species, str, Component]],
        direction: str = 'exporter',
        internal_compartment: Union[str, Compartment] = 'internal',
        external_compartment: Union[str, Compartment] = 'external',
        ATP: int = None,
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Pump
        if isinstance(membrane_pump, Component):
            self.membrane_pump = membrane_pump.product
        elif isinstance(membrane_pump, Species):
            self.membrane_pump = membrane_pump
        elif isinstance(membrane_pump, str):
            self.membrane_pump = self.set_species(
                membrane_pump,
                material_type='protein',
                attributes=attributes,
            )
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(membrane_pump)}")

        self.membrane_pump.direction = direction

        # ENERGY and WASTE
        self.energy = self.set_species(
            'ATP',
            material_type='small_molecule',
            compartment=internal_compartment,
            attributes=attributes,
        )
        self.waste = self.set_species(
            'ADP',
            material_type='small_molecule',
            compartment=internal_compartment,
            attributes=attributes,
        )

        if ATP is None:
            self.membrane_pump.ATP = 1
        else:
            self.membrane_pump.ATP = ATP


        if substrate is None:
            substrate = []
        elif not isinstance(substrate, list):
            substrate = [substrate]

        # Initialize lists for substrate and products
        self.substrate = []
        self.product = []

        if direction == 'importer':
            sub_compartment = external_compartment
            prod_compartment = internal_compartment
        elif direction == 'exporter':
            sub_compartment = internal_compartment
            prod_compartment = external_compartment
        else:
            raise TypeError(
                f"Direction of pump must be defined as 'exporter' or" \
                f"'importer'. Got: {direction}.",
                UserWarning,
            )

        # Iterate over each substrate
        for sub in substrate:
            if isinstance(sub, Species):
                # sub_name = sub.name
                sub_species = sub
                if sub.compartment.name != sub_compartment:
                    sub.compartment = sub_compartment

            elif isinstance(sub, Component):
                sub_name = sub.get_species()[0].name

                sub_species = self.set_species(
                    sub_name,
                    compartment=sub_compartment,
                    attributes=attributes,
                )

            elif isinstance(sub, str):
                sub_name = sub
                sub_species = self.set_species(
                    sub,
                    compartment=sub_compartment,
                    attributes=attributes,
                )

            else:
                raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(sub)}")

            prod_species = self.set_species(
                sub_species.name,
                material_type=sub_species.material_type,
                compartment=prod_compartment,
                attributes=sub_species.attributes,
            )

            self.substrate.append(sub_species)
            self.product.append(prod_species)

        # Name the component
        name = (
            self.membrane_pump.name
            + '_'
            + self.membrane_pump.compartment.name
        )

        Component.__init__(
            self=self, name=name,
            mechanisms={'transport':Transport_PrimaryActive_ABCexporter()},
            **kwargs)

    def get_species(self):
        """Get the membrane pump protein species.

        Returns
        -------
        Species
            The membrane pump protein species.

        """
        return self.membrane_pump

    def update_species(self):
        r"""Uses the 'transport' mechanism to generate pump-mediated species.

        Uses the 'transport' mechanism to generate species including the
        pump protein, substrate, product, ATP, and ADP.

        Returns
        -------
        list of Species
            List of species generated by the transport mechanism,
            including pump, substrate, product, energy, and waste.

        """
        mech_tra = self.get_mechanism('transport')
        species_list = []

        for sub, prod in zip(self.substrate, self.product):
            species_list.extend(
                mech_tra.update_species(
                    self.membrane_pump,
                    sub,
                    prod,
                    self.energy,
                    self.waste,
                )
            )
        return list(set(species_list))

    def update_reactions(self):
        """Uses the 'transport' mechanism to generate pump-mediated reactions.

        Uses the 'transport' mechanism to generate reactions for active
        transport coupled to ATP hydrolysis.

        Returns
        -------
        list of Reaction
            List of ATP-dependent transport reactions.

        """
        mech_tra = self.get_mechanism('transport')
        reactions_list = []

        for sub, prod in zip(self.substrate, self.product):
            reactions_list.extend(
                mech_tra.update_reactions(
                    self.membrane_pump,
                    sub,
                    prod,
                    self.energy,
                    self.waste,
                    component=self,
                    part_id=self.name,
                )
            )
        return reactions_list


class MembraneSensor(Component):
    r"""Two-component system (TCS) membrane sensor.

    A `MembraneSensor` component represents a membrane sensor protein in a
    two-component system. The sensor detects external signal substrates
    and catalyzes the transfer of a chemical group (typically phosphate)
    to a response protein, activating it. The component uses the
    'membrane_sensor' mechanism to generate signal transduction reactions.

    Parameters
    ----------
    membrane_sensor : Species, str, or Component
        The membrane sensor protein (histidine kinase) that detects the
        signal. Can be a `Species` object, string name, or `Component`.
    response_protein : Species, str, or Component
        The cytoplasmic response regulator protein that receives the
        signal. Can be a `Species` object, string name, or `Component`.
    assigned_substrate : Species, str, or Component
        The chemical group to be transferred (typically phosphate). Can be
        a `Species` object, string name, or `Component`.
    signal_substrate : Species, str, or Component
        The external signal molecule that activates the sensor. Can be a
        `Species` object, string name, or `Component`.
    product : Species, str, or Component, optional
        The activated response protein product. If None, automatically
        named as '<response_protein>active'.
    internal_compartment : str or Compartment, default='internal'
        The internal compartment containing response protein. Can be a
        string name (creates new Compartment) or an existing `Compartment`
        object.
    external_compartment : str or Compartment, default='external'
        The external compartment containing signal. Can be a string name
        (creates new Compartment) or an existing `Compartment` object.
    ATP : int, default=2
        Number of ATP molecules required for the signaling process.
    attributes : list of str, optional
        List of attribute tags to associate with species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    membrane_sensor : Species
        The membrane sensor protein species.
    response_protein : Species
        The response regulator protein species.
    assigned_substrate : Species
        The substrate to be transferred (e.g., phosphate).
    signal_substrate : Species
        The external signal molecule species.
    product : Species
        The activated response protein species.
    energy : Species
        ATP species used for energy.
    waste : Species
        ADP species produced.

    See Also
    --------
    MembraneChannel : Membrane channel for substrate diffusion.
    MembraneCarrier : Membrane carrier for substrate diffusion or transport.
    Component : Base class for biomolecular components.

    Notes
    -----
    Two-component systems (TCS) are common bacterial signal transduction
    pathways. The typical mechanism involves:

    1. Signal detection by membrane sensor (histidine kinase)
    2. Autophosphorylation of sensor using ATP
    3. Phosphotransfer to response regulator
    4. Activated response regulator regulates gene expression

    The general reaction scheme:
    $$
        & 'signal' + 'sensor' + 'ATP' + 'response_protein' \\
        &     --> 'signal' + 'sensor' + 'ADP' + 'response_protein-P'
    $$

    The component name is automatically generated as:
    '<membrane_sensor_name>_<compartment_name>'

    Examples
    --------
    Create a simple two-component system:

    >>> tcs = bcp.MembraneSensor(
    ...     membrane_sensor='EnvZ',
    ...     response_protein='OmpR',
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Osmolarity',
    ...     ATP=2
    ... )

    Create a chemotaxis receptor:

    >>> chemoreceptor = bcp.MembraneSensor(
    ...     membrane_sensor='CheA',
    ...     response_protein='CheY',
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Aspartate',
    ...     product='CheY_P'
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[tcs],
    ...     mechanisms={
    ...         'membrane_sensor': bcp.Membrane_Signaling_Pathway_MM()},
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_sensor: Union[Species, str, Component],
        response_protein: Union[Species, str, Component],
        assigned_substrate: Union[Species, str, Component],
        signal_substrate: Union[Species, str, Component],
        product: Union[Species, str, Component] = None,
        internal_compartment: Union[str, Compartment] = 'internal',
        external_compartment: Union[str, Compartment] = 'external',
        ATP: int = 2,
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # RESPONSE PROTEIN
        if response_protein is None:
            self.response_protein = None
        elif isinstance(response_protein, Component):
            self.response_protein = self.set_species(
                response_protein.get_species()[0].name,
                compartment=internal_compartment,
                attributes=attributes,
            )
        else:
            self.response_protein = self.set_species(
                response_protein,
                compartment=internal_compartment,
                attributes=attributes,
            )

        # PRODUCT PROTEIN
        if product is None:
            self.product = self.set_species(
                self.response_protein.name + 'active',
                compartment=internal_compartment,
                attributes=attributes,
            )
        else:
            self.product = self.set_species(
                product,
                compartment=internal_compartment,
                attributes=attributes,
            )

        # ASSIGNED SUBSTRATE
        if assigned_substrate is None:
            self.assigned_substrate = None
        elif isinstance(assigned_substrate, Component):
            self.assigned_substrate = self.set_species(
                assigned_substrate.get_species()[0].name,
                compartment=internal_compartment,
                attributes=attributes,
            )
        else:
            self.assigned_substrate = self.set_species(
                assigned_substrate,
                compartment=internal_compartment,
                attributes=attributes,
            )
        # SIGNAL SUBSTRATE
        if signal_substrate is None:
            self.signal_substrate = None
        elif isinstance(signal_substrate, Component):
            self.signal_substrate = self.set_species(
                signal_substrate.get_species()[0].name,
                compartment=internal_compartment,
                attributes=attributes,
            )
        else:
            self.signal_substrate = self.set_species(
                signal_substrate,
                compartment=internal_compartment,
                attributes=attributes,
            )

        # PROTEIN
        if isinstance(membrane_sensor, Component):
            self.membrane_sensor = membrane_sensor.product
        elif isinstance(membrane_sensor, Species):
            self.membrane_sensor = membrane_sensor
        elif isinstance(membrane_sensor, str):
            self.membrane_sensor = self.set_species(
                membrane_sensor,
                material_type='protein',
                attributes=attributes,
            )
        else:
            raise TypeError(f"Expected Species, Component, or str. "
                            f"Got: {type(membrane_sensor)}")

        # ENERGY: ATP
        if ATP is None:
            self.membrane_sensor.ATP = 1
        else:
            self.membrane_sensor.ATP = ATP

        self.energy = self.set_species(
            'ATP',
            material_type='small_molecule',
            compartment=internal_compartment,
            attributes=attributes,
        )
        self.waste = self.set_species(
            'ADP',
            material_type='small_molecule',
            compartment=internal_compartment,
            attributes=attributes,
        )
        # Name the component
        name = (
            self.membrane_sensor.name
            + '_'
            + self.membrane_sensor.compartment.name
        )

        Component.__init__(
            self=self, name=name,
            mechanisms={'membrane_sensor':Sensor_TwoComponentSystem()},
            **kwargs)

    def get_species(self):
        """Get the membrane sensor protein species.

        Returns
        -------
        Species
            The membrane sensor protein (histidine kinase) species.

        """
        return self.membrane_sensor

    def update_species(self):
        """Use 'membrane_sensor' to generate species signaling species.

        Uses the 'membrane_sensor' mechanism to generate all species
        involved in the signaling pathway including sensor, response
        protein, substrates, signal, product, ATP, and ADP.

        Returns
        -------
        list of Species
            List of species generated by the membrane_sensor mechanism.

        """
        mech_sen = self.get_mechanism('membrane_sensor')
        return mech_sen.update_species(
            self.membrane_sensor,
            self.response_protein,
            self.assigned_substrate,
            self.signal_substrate,
            self.product,
            self.energy,
            self.waste,
        )

    def update_reactions(self):
        """Use 'membrane_sensor' to generate species signaling reactions.

        Uses the 'membrane_sensor' mechanism to generate reactions for
        signal detection, ATP-dependent phosphorylation, and
        phosphotransfer to the response regulator.

        Returns
        -------
        list of Reaction
            List of signal transduction reactions including sensing,
            autophosphorylation, and phosphotransfer.

        """
        mech_sen = self.get_mechanism('membrane_sensor')
        return mech_sen.update_reactions(
            self.membrane_sensor,
            self.response_protein,
            self.assigned_substrate,
            self.signal_substrate,
            self.product,
            self.energy,
            self.waste,
            component=self,
            part_id=self.name,
        )
