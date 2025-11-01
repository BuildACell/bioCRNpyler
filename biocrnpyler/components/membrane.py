#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import warnings
from typing import Union

from ..core.compartment import Compartment
from ..core.component import Component
from ..core.species import Species


class DiffusibleMolecule(Component):
    """Molecule that diffuses passively through a membrane.

    A `DiffusibleMolecule` component represents a molecule that undergoes
    passive diffusion across a membrane between two compartments. The
    component uses a 'diffusion' mechanism to generate bidirectional
    diffusion reactions based on concentration gradients.

    Parameters
    ----------
    substrate : Species, str, or Component
        The diffusible molecule species. Can be a `Species` object, string
        name, or `Component` with an associated species.
    internal_compartment : str or Compartment, default='Internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='External'
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
    MembraneChannel : Active transport through membrane channels.
    MembranePump : ATP-dependent active transport.
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
    ...     internal_compartment='Cytoplasm',
    ...     external_compartment='Extracellular'
    ... )

    Use with a mixture and diffusion mechanism:

    >>> mixture = bcp.Mixture(
    ...     components=[glucose],
    ...     mechanisms={'diffusion': bcp.SimpleDiffusion()},
    ...     parameters={'k_diff': 0.01}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        substrate: Union[Species, str, Component],
        internal_compartment: Union[str, Compartment] = 'Internal',
        external_compartment: Union[str, Compartment] = 'External',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Substrate
        self.substrate = self.set_species(
            substrate, compartment=internal_compartment
        )
        self.product = self.set_species(
            substrate, compartment=external_compartment
        )

        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = self.substrate.name + '_' + self.substrate.compartment.name

        Component.__init__(
            self=self, name=name, attributes=attributes, **kwargs
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
        """Generate species for diffusion reactions.

        Uses the 'diffusion' mechanism to generate species in both
        compartments.

        Returns
        -------
        list of Species
            List of species in internal and external compartments generated
            by the diffusion mechanism.

        """
        mech_diff = self.get_mechanism('diffusion')
        return mech_diff.update_species(self.substrate, self.product)

    def update_reactions(self):
        """Generate bidirectional diffusion reactions.

        Uses the 'diffusion' mechanism to generate reactions for passive
        diffusion between compartments.

        Returns
        -------
        list of Reaction
            List of diffusion reactions (forward and reverse) between
            internal and external compartments.

        """
        mech_diff = self.get_mechanism('diffusion')
        return mech_diff.update_reactions(
            self.substrate, self.product, component=self, part_id=self.name
        )


class IntegralMembraneProtein(Component):
    """Transmembrane protein that integrates into the membrane.

    An `IntegralMembraneProtein` component represents a membrane protein
    that integrates into a membrane compartment. The component uses a
    'membrane_insertion' mechanism to generate reactions for protein
    insertion into the membrane. The size parameter allows modeling of
    oligomeric channels (dimers, trimers, etc.).

    Parameters
    ----------
    membrane_protein : Species, str, or Component
        The membrane protein species before insertion. Can be a `Species`
        object, string name, or `Component` with an associated species.
    product : Species, str, or Component
        The integrated membrane protein species. Can be a `Species` object,
        string name, or `Component`.
    direction : str, optional
        Transport direction attribute for the integrated protein.
        Default is 'passive'. Common values: 'passive', 'importer',
        'exporter'.
    size : int, optional
        Number of monomers needed to form the functional channel. Used to
        model oligomeric channels (e.g., size=2 for dimers, size=3 for
        trimers). Default is 1.
    compartment : str or Compartment, default='Internal'
        The compartment containing the membrane protein before insertion.
        Can be a string name or `Compartment` object.
    membrane_compartment : str or Compartment, default='Membrane'
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
    MembraneChannel : Membrane channel for substrate transport.
    Component : Base class for biomolecular components.

    Notes
    -----
    The membrane_insertion mechanism generates reactions for protein
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
    ...     direction='passive'
    ... )

    Create a dimeric channel protein:

    >>> dimer = bcp.IntegralMembraneProtein(
    ...     membrane_protein='Aquaporin',
    ...     product='Aquaporin_channel',
    ...     size=2,
    ...     direction='passive'
    ... )

    """

    def __init__(
        self,
        membrane_protein: Union[Species, str, Component],
        product: Union[Species, str, Component],
        direction: str = None,
        size: int = None,
        compartment: Union[str, Compartment] = 'Internal',
        membrane_compartment: Union[str, Compartment] = 'Membrane',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(compartment, str):
            compartment = Compartment(name=compartment)
        if isinstance(membrane_compartment, str):
            membrane_compartment = Compartment(name=membrane_compartment)

        # Logic for prioritizing compartments
        self.membrane_protein = self.set_species(membrane_protein)
        if self.membrane_protein.compartment.name == 'default':
            self.membrane_protein.compartment = compartment
        elif (
            self.membrane_protein.compartment.name != compartment.name
            and compartment.name == 'Internal'
        ):
            warnings.warn(
                "Inconsistent compartments, prioritizing membrane protein "
                "compartment.",
                UserWarning,
            )
            compartment = self.membrane_protein.compartment
        else:
            warnings.warn(
                "Inconsistent compartments, prioritizing integral membrane "
                "protein compartment.",
                UserWarning,
            )
            self.membrane_protein.compartment = compartment

        # PROTEIN
        self.membrane_protein = self.set_species(
            membrane_protein,
            material_type='protein',
            compartment=compartment,
            attributes=attributes,
        )

        # PRODUCT is an integrated membrane protein (transmembrane_protein)
        if product is None:
            if direction is None:
                self.product = self.set_species(
                    product,
                    material_type='protein',
                    compartment=membrane_compartment,
                    attributes=['passive'],
                )
            else:
                self.product = self.set_species(
                    product,
                    material_type='protein',
                    compartment=membrane_compartment,
                    attributes=[direction],
                )
        else:
            if direction is None:
                self.product = self.set_species(
                    product,
                    material_type='protein',
                    compartment=membrane_compartment,
                    attributes=['passive'],
                )
            else:
                self.product = self.set_species(
                    product,
                    material_type='protein',
                    compartment=membrane_compartment,
                    attributes=[direction],
                )

        # Indicates the number of monomers that compose the channel,
        # will be used in Membrane_Protein_Integration(Mechanism)
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

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self):
        """Get the membrane protein species before insertion.

        Returns
        -------
        Species
            The membrane protein species in the compartment before
            integration into the membrane.

        """
        return self.membrane_protein

    def update_species(self):
        """Generate species for membrane protein insertion.

        Uses the 'membrane_insertion' mechanism to generate species for
        the protein before and after insertion.

        Returns
        -------
        list of Species
            List of species generated by the membrane_insertion mechanism,
            including the protein and integrated product.

        """
        mech_ins = self.get_mechanism('membrane_insertion')
        return mech_ins.update_species(self.membrane_protein, self.product)

    def update_reactions(self):
        """Generate reactions for membrane protein insertion.

        Uses the 'membrane_insertion' mechanism to generate reactions for
        protein integration into the membrane.

        Returns
        -------
        list of Reaction
            List of reactions for protein insertion into the membrane.

        """
        mech_ins = self.get_mechanism('membrane_insertion')
        return mech_ins.update_reactions(
            self.membrane_protein,
            self.product,
            component=self,
            part_id=self.name,
        )


class MembraneChannel(Component):
    """Membrane channel for facilitated transport across membranes.

    A `MembraneChannel` component represents a membrane channel or
    transporter that facilitates substrate movement across a membrane
    following concentration gradients. The direction of transport depends
    on the specific transporter type. The component uses a 'transport'
    mechanism to generate transport reactions.

    Parameters
    ----------
    integral_membrane_protein : Species, str, or Component
        The integral membrane protein that forms the channel. Can be a
        `Species` object, string name, or `Component`. If a string,
        automatically creates a protein species with appropriate direction
        attribute.
    substrate : Species, str, or Component
        The substrate to be transported through the channel. Can be a
        `Species` object, string name, or `Component`.
    direction : str, optional
        Direction of transport. If None, extracted from
        integral_membrane_protein attributes. Common values: 'importer'
        (external --> internal), 'exporter' (internal --> external),
        'passive' (bidirectional).
    internal_compartment : str or Compartment, default='Internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='External'
        The external compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    attributes : list of str, optional
        List of attribute tags to associate with substrate species.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    integral_membrane_protein : Species
        The membrane channel protein species.
    substrate : Species
        The substrate species in the source compartment (depends on
        direction).
    product : Species
        The same substrate in the destination compartment.

    See Also
    --------
    IntegralMembraneProtein : Protein insertion into membranes.
    MembranePump : ATP-dependent active transport.
    DiffusibleMolecule : Passive diffusion without channels.
    Component : Base class for biomolecular components.

    Notes
    -----
    The transport mechanism generates reactions based on the direction:

    - 'importer': substrate_external + channel
          --> substrate_internal + channel
    - 'exporter': substrate_internal + channel
          --> substrate_external + channel
    - 'passive': bidirectional transport following gradients

    The component name is automatically generated as:
    '<integral_membrane_protein_name>_<compartment_name>'

    Examples
    --------
    Create a glucose importer:

    >>> importer = bcp.MembraneChannel(
    ...     integral_membrane_protein='GlucoseTransporter',
    ...     substrate='Glucose',
    ...     direction='importer'
    ... )

    Create a passive channel:

    >>> channel = bcp.MembraneChannel(
    ...     integral_membrane_protein='WaterChannel',
    ...     substrate='Water',
    ...     direction='passive'
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[importer],
    ...     mechanisms={'transport': bcp.FacilitatedTransport_MM()},
    ...     parameter_file='mechanisms/transport_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        integral_membrane_protein: Union[Species, str, Component],
        substrate: Union[Species, str, Component],
        direction: str = None,
        internal_compartment: Union[str, Compartment] = 'Internal',
        external_compartment: Union[str, Compartment] = 'External',
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # Set up the integral membrane protein
        # TODO: allow integral_membrane_protein to be a Component
        if isinstance(integral_membrane_protein, str):
            integral_membrane_protein = self.set_species(
                integral_membrane_protein,
                material_type='protein',
                attributes='passive' if direction is None else direction,
            )
        self.integral_membrane_protein = integral_membrane_protein

        # Get the direction from the integral_membrane_protein, if not given
        # TODO: need more complete check for conflicting information
        if direction is None:
            if 'importer' in integral_membrane_protein.attributes:
                direction = 'importer'
            elif 'exporter' in integral_membrane_protein.attributes:
                direction = 'exporter'

        # Substrate and product assignments.
        #
        # In the case of membrane components, the substrate is the
        # substance on which the transporter/channel acts without
        # distinction of compartment.  The substrate and product are the
        # same substance and the substance does not change as a result
        # except for the compartment.  Therefore, the product here is based
        # on the action of the transporter.
        #
        # TODO: if the substrate is not passed as a string, we need more
        # sophisticated logic to set up the product, since the `set_species`
        # method will just return the existing species, without changing
        # the compartment.
        #
        # TODO: we should think about allowing a list of substrates to
        # be supplied and/or some sort of attribute based import (eg,
        # based on size).

        if substrate is None:
            self.substrate = None

        else:
            if direction == 'importer':
                self.substrate = self.set_species(
                    substrate,
                    compartment=external_compartment,
                    attributes=attributes,
                )
                self.product = self.set_species(
                    self.substrate.name,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
            else:
                self.substrate = self.set_species(
                    substrate,
                    compartment=internal_compartment,
                    attributes=attributes,
                )
                self.product = self.set_species(
                    self.substrate.name,
                    compartment=external_compartment,
                    attributes=attributes,
                )

        # Name the component
        name = (
            self.integral_membrane_protein.name
            + '_'
            + self.integral_membrane_protein.compartment.name
        )

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self):
        """Get the integral membrane protein species.

        Returns
        -------
        Species
            The integral membrane protein species that forms the channel.

        """
        return self.integral_membrane_protein

    def update_species(self):
        """Generate species for channel-mediated transport.

        Uses the 'transport' mechanism to generate species including the
        channel protein, substrate, and product.

        Returns
        -------
        list of Species
            List of species generated by the transport mechanism.

        """
        mech_tra = self.get_mechanism('transport')
        return mech_tra.update_species(
            self.integral_membrane_protein, self.substrate, self.product
        )

    def update_reactions(self):
        """Generate reactions for channel-mediated transport.

        Uses the 'transport' mechanism to generate reactions for substrate
        transport through the channel.

        Returns
        -------
        list of Reaction
            List of transport reactions through the membrane channel.

        """
        mech_tra = self.get_mechanism('transport')
        return mech_tra.update_reactions(
            self.integral_membrane_protein,
            self.substrate,
            self.product,
            component=self,
            part_id=self.name,
        )


class MembranePump(Component):
    """ATP-dependent membrane pump for active transport.

    A `MembranePump` component represents an active transporter or pump
    that uses ATP to transport substrates across membranes against
    concentration gradients. The pump operates unidirectionally and requires
    energy in the form of ATP. The component uses a 'transport' mechanism
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
        'passive' (default). Affects substrate and ATP compartment
        placement.
    internal_compartment : str or Compartment, default='Internal'
        The internal compartment. Can be a string name (creates new
        Compartment) or an existing `Compartment` object.
    external_compartment : str or Compartment, default='External'
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
    MembraneChannel : Facilitated transport without ATP.
    DiffusibleMolecule : Passive diffusion.
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
    ...     direction='exporter',
    ...     ATP=2
    ... )

    Create an ABC transporter (importer):

    >>> abc = bcp.MembranePump(
    ...     membrane_pump='ABC_Transporter',
    ...     substrate='Maltose',
    ...     direction='importer',
    ...     ATP=1
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[pump],
    ...     mechanisms={'transport': bcp.PrimaryActiveTransport_MM()},
    ...     parameter_file='mechanisms/transport_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_pump: Union[Species, str, Component],
        substrate: Union[Species, str, Component],
        direction: str = None,
        internal_compartment: Union[str, Compartment] = 'Internal',
        external_compartment: Union[str, Compartment] = 'External',
        ATP: int = None,
        attributes=None,
        **kwargs,
    ):
        # Creates compartment object if compartment is a str
        if isinstance(internal_compartment, str):
            internal_compartment = Compartment(name=internal_compartment)
        if isinstance(external_compartment, str):
            external_compartment = Compartment(name=external_compartment)

        # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            product = substrate
            self.substrate = self.set_species(
                substrate,
                compartment=internal_compartment,
                attributes=attributes,
            )
            self.product = self.set_species(
                product,
                compartment=external_compartment,
                attributes=attributes,
            )

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

        # PROTEIN
        if isinstance(membrane_pump, str):
            if ATP is None:
                ATP = 1
            else:
                ATP = ATP

            if direction is None:
                self.membrane_pump = self.set_species(
                    membrane_pump,
                    material_type='protein',
                    attributes='passive',
                )
                self.membrane_pump.ATP = ATP
            else:
                self.membrane_pump = self.set_species(
                    membrane_pump,
                    material_type='protein',
                    attributes=direction,
                )
                self.membrane_pump.ATP = ATP
                if direction == 'importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product = substrate
                        self.substrate = self.set_species(
                            substrate,
                            compartment=external_compartment,
                            attributes=attributes,
                        )
                        self.product = self.set_species(
                            product,
                            compartment=internal_compartment,
                            attributes=attributes,
                        )
        else:
            if membrane_pump.attributes[0] == 'passive':
                self.integral_membrane_protein = self.set_species(
                    membrane_pump,
                    material_type='protein',
                    attributes='passive',
                )
            elif membrane_pump.attributes[0] == 'exporter':
                self.membrane_pump = self.set_species(
                    membrane_pump,
                    material_type='protein',
                    attributes='exporter',
                )
            elif membrane_pump.attributes[0] == 'importer':
                self.membrane_pump = self.set_species(
                    membrane_pump,
                    material_type='protein',
                    attributes='importer',
                )
                self.energy = self.set_species(
                    'ATP',
                    material_type='small_molecule',
                    compartment=external_compartment,
                    attributes=attributes,
                )
                self.waste = self.set_species(
                    'ADP',
                    material_type='small_molecule',
                    compartment=external_compartment,
                    attributes=attributes,
                )
                if substrate is None:
                    self.substrate = None
                else:
                    product = substrate
                    self.substrate = self.set_species(
                        substrate,
                        compartment=external_compartment,
                        attributes=attributes,
                    )
                    self.product = self.set_species(
                        product,
                        compartment=internal_compartment,
                        attributes=attributes,
                    )
            else:
                print('Membrane channel direction not found.')

            if ATP is None:
                self.membrane_pump.ATP = 1
            else:
                self.membrane_pump.ATP = ATP

        # Name the component
        name = (
            self.membrane_pump.name
            + '_'
            + self.membrane_pump.compartment.name
        )

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self):
        """Get the membrane pump protein species.

        Returns
        -------
        Species
            The membrane pump protein species.

        """
        return self.membrane_pump

    def update_species(self):
        """Generate species for ATP-dependent transport.

        Uses the 'transport' mechanism to generate species including the
        pump protein, substrate, product, ATP, and ADP.

        Returns
        -------
        list of Species
            List of species generated by the transport mechanism,
            including pump, substrate, product, energy, and waste.

        """
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_species(
            self.membrane_pump,
            self.substrate,
            self.product,
            self.energy,
            self.waste,
        )

    def update_reactions(self):
        """Generate reactions for ATP-dependent transport.

        Uses the 'transport' mechanism to generate reactions for active
        transport coupled to ATP hydrolysis.

        Returns
        -------
        list of Reaction
            List of ATP-dependent transport reactions.

        """
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_reactions(
            self.membrane_pump,
            self.substrate,
            self.product,
            self.energy,
            self.waste,
            component=self,
            part_id=self.name,
        )


class MembraneSensor(Component):
    """Two-component system (TCS) membrane sensor protein.

    A `MembraneSensor` component represents a membrane sensor protein in a
    two-component signaling system. The sensor detects external signal
    substrates and catalyzes the transfer of a chemical group (typically
    phosphate) to a response protein, activating it. The component uses a
    'membrane_sensor' mechanism to generate signal transduction reactions.

    Parameters
    ----------
    membrane_sensor_protein : Species, str, or Component
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
    internal_compartment : str or Compartment, default='Internal'
        The internal compartment containing response protein. Can be a
        string name (creates new Compartment) or an existing `Compartment`
        object.
    external_compartment : str or Compartment, default='External'
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
    membrane_sensor_protein : Species
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

        signal + sensor + ATP + response_protein -->
        signal + sensor + ADP + response_protein-P

    The component name is automatically generated as:
    '<membrane_sensor_protein_name>_<compartment_name>'

    Examples
    --------
    Create a simple two-component system:

    >>> tcs = bcp.MembraneSensor(
    ...     membrane_sensor_protein='EnvZ',
    ...     response_protein='OmpR',
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Osmolarity',
    ...     ATP=2
    ... )

    Create a chemotaxis receptor:

    >>> chemoreceptor = bcp.MembraneSensor(
    ...     membrane_sensor_protein='CheA',
    ...     response_protein='CheY',
    ...     assigned_substrate='Phosphate',
    ...     signal_substrate='Aspartate',
    ...     product='CheY_P'
    ... )

    Use with a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[tcs],
    ...     mechanisms={'membrane_sensor': bcp.MembraneSignalingPathway_MM()},
    ...     parameter_file='mechanisms/transport_parameters.tsv'
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        membrane_sensor_protein: Union[Species, str, Component],
        response_protein: Union[Species, str, Component],
        assigned_substrate: Union[Species, str, Component],
        signal_substrate: Union[Species, str, Component],
        product: Union[Species, str, Component] = None,
        internal_compartment: Union[str, Compartment] = 'Internal',
        external_compartment: Union[str, Compartment] = 'External',
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
        else:
            self.response_protein = self.set_species(
                response_protein,
                compartment=internal_compartment,
                attributes=attributes,
            )

        # PRODUCT PROTEIN
        if product is None:
            self.product = self.set_species(
                str(response_protein) + 'active',
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
        else:
            self.assigned_substrate = self.set_species(
                assigned_substrate,
                compartment=internal_compartment,
                attributes=attributes,
            )
        # SIGNAL SUBSTRATE
        if signal_substrate is None:
            self.signal_substrate = None
        else:
            self.signal_substrate = self.set_species(
                signal_substrate,
                compartment=internal_compartment,
                attributes=attributes,
            )

        # PROTEIN
        if membrane_sensor_protein is None:
            self.membrane_sensor_protein = None
        else:
            self.membrane_sensor_protein = self.set_species(
                membrane_sensor_protein,
                material_type='protein',
                attributes=attributes,
            )
        # ENERGY: ATP
        if ATP is None:
            self.membrane_sensor_protein.ATP = 1
        else:
            self.membrane_sensor_protein.ATP = ATP

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
            self.membrane_sensor_protein.name
            + '_'
            + self.membrane_sensor_protein.compartment.name
        )

        Component.__init__(self=self, name=name, **kwargs)

    def get_species(self):
        """Get the membrane sensor protein species.

        Returns
        -------
        Species
            The membrane sensor protein (histidine kinase) species.

        """
        return self.membrane_sensor_protein

    def update_species(self):
        """Generate species for two-component signal transduction.

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
            self.membrane_sensor_protein,
            self.response_protein,
            self.assigned_substrate,
            self.signal_substrate,
            self.product,
            self.energy,
            self.waste,
        )

    def update_reactions(self):
        """Generate reactions for two-component signal transduction.

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
            self.membrane_sensor_protein,
            self.response_protein,
            self.assigned_substrate,
            self.signal_substrate,
            self.product,
            self.energy,
            self.waste,
            component=self,
            part_id=self.name,
        )
