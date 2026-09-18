# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import warnings

from ..core.mechanism import Mechanism
from ..core.propensities import GeneralPropensity, ProportionalHillNegative
from ..core.reaction import Reaction
from ..core.species import Complex, Species


class Diffusion_Simple(Mechanism):
    """Passive diffusion mechanism for substrate transport across membranes.

    A 'diffusion' mechanism that models simple passive diffusion of
    substrates through a membrane without requiring membrane proteins or
    energy. The transport is bidirectional and follows Fick's law of
    diffusion with equal forward and reverse rate constants.

    The reaction follows the schema:
    $$
        'substrate' <--> 'product'
    $$

    where substrate and product represent the same species on opposite sides
    of the membrane.

    Parameters
    ----------
    name : str, default='diffusion_simple'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='diffusion'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('diffusion').

    See Also
    --------
    Diffusion_Facilitated_Channel : Passive transport via membrane channels.
    Diffusion_Facilitated_Carrier : Facilitated diffusion with carriers.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    Simple diffusion models the movement of small, lipophilic molecules
    across lipid bilayers without the assistance of membrane proteins. This
    process is driven purely by concentration gradients and does not require
    cellular energy.

    Common examples include:

    - Diffusion of gases (O2, CO2) across cell membranes
    - Transport of small nonpolar molecules
    - Movement of lipid-soluble substances

    The mechanism generates a single reversible mass-action reaction with
    equal forward and reverse rate constants, reflecting the thermodynamic
    equilibrium of passive diffusion.

    Required parameters for this mechanism:

    - 'k_diff' : Diffusion rate constant (same for both directions)

    Examples
    --------
    Model oxygen diffusion across a membrane:

    >>> O2 = bcp.DiffusibleMolecule('O2')
    >>> mechanism = bcp.Diffusion_Simple()
    >>> mixture = bcp.Mixture(
    ...     components=[O2],
    ...     mechanisms={'diffusion': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='diffusion_simple',
        mechanism_type='diffusion',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(self, substrate, product, **kwargs):
        """Generate species for simple diffusion.

        Returns the substrate and product species involved in the diffusion
        reaction.

        Parameters
        ----------
        substrate : Species
            The substrate species on one side of the membrane (typically
            the intracellular side).
        product : Species
            The product species on the other side of the membrane (typically
            the extracellular side). Usually the same molecular species as
            substrate but in a different compartment.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [substrate, product].

        """
        return [substrate, product]

    def update_reactions(
        self,
        substrate,
        product,
        component=None,
        part_id=None,
        k_diff=None,
        **kwargs,
    ):
        """Generate reactions for simple diffusion.

        Creates a single reversible mass-action reaction representing
        passive diffusion across a membrane with equal forward and reverse
        rate constants.

        Parameters
        ----------
        substrate : Species
            The substrate species on one side of the membrane.
        product : Species
            The product species on the other side of the membrane.
        component : Component, optional
            Component containing parameter values. Required if k_diff is not
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is
            provided, defaults to component.name.
        k_diff : Parameter or float, optional
            Diffusion rate constant. If None, retrieved from component
            parameters. Used as both forward and reverse rate constant.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single reversible mass-action reaction for
            diffusion.

        Raises
        ------
        ValueError
            If component is None and k_diff is not provided.

        Notes
        -----
        The reaction has equal forward and reverse rate constants, reflecting
        the thermodynamic equilibrium of passive diffusion:

        1. substrate <--> product (rates: 'k_diff')

        This method requires both component and part_id parameters to
        retrieve rate constants from the component's parameter database.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_diff is None):
            raise ValueError("Must pass in a Component or values for k_diff.")
        if k_diff is None:
            k_diff = component.get_parameter(
                'k_diff', part_id=part_id, mechanism=self
            )
        else:
            k_diff = k_diff

        # Simple diffusion
        # Sub (Internal) <--> Product (External)
        diffusion_rxn = Reaction.from_massaction(
            inputs=[substrate],
            outputs=[product],
            k_forward=k_diff,
            k_reverse=k_diff,
        )
        return [diffusion_rxn]


class Integration_MembraneProtein(Mechanism):
    """Membrane protein integration mechanism for protein insertion.

    A 'membrane_integration' mechanism that models the integration of newly
    synthesized proteins into cellular membranes. Supports both monomeric
    and oligomeric membrane proteins, handling oligomerization before
    membrane insertion when required.

    The reaction schema depends on protein oligomeric state:

    For monomers (size = 1):
    $$
        'monomer' --> 'integral membrane protein'
    $$

    For oligomers (size > 1):
    $$
        'monomer' * 'size' <--> 'oligomer' --> 'integral membrane protein'
    $$

    Parameters
    ----------
    name : str, default='integration_membraneprotein'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='membrane_integration'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('membrane_integration').

    See Also
    --------
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models the process by which proteins become embedded in
    cellular membranes. For oligomeric proteins, multiple monomers must
    first associate into a complex before integration can occur. The
    integration step uses a `ProportionalHillNegative` propensity function to
    model saturation kinetics and product inhibition.

    The mechanism requires the integral membrane protein to have a size
    attribute (integral_membrane_protein.size) that specifies the number of
    monomers in the functional unit.

    Common examples include:

    - Integration of ion channels (often oligomeric)
    - Insertion of receptor proteins (can be monomeric or oligomeric)
    - Assembly and insertion of transporter complexes

    Required parameters for this mechanism:

    - 'kb_oligomer' : Forward oligomerization rate constant (for size > 1)
    - 'ku_oligomer' : Reverse oligomerization rate constant (for size > 1)
    - 'kex' : Maximum integration rate constant
    - 'kcat' : Michaelis constant for integration

    Examples
    --------
    Model integration of a tetrameric channel:

    >>> channel = bcp.IntegralMembraneProtein(
    ...     membrane_protein='Aquaporin',
    ...     product='Aquaporin_channel',
    ...     size=2,
    ... )
    >>> mechanism = bcp.Integration_MembraneProtein()
    >>> mixture = bcp.Mixture(
    ...     components=[channel],
    ...     mechanisms={'membrane_integration': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='integration_membraneprotein',
        mechanism_type='membrane_integration',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self, integral_membrane_protein, product, complex=None, **kwargs
    ):
        """Generate species for membrane protein integration.

        Creates species for monomers, oligomeric complexes (if needed), and
        the integrated membrane protein product.

        Parameters
        ----------
        integral_membrane_protein : Species
            The membrane protein monomer that will be integrated. Must have
            a size attribute specifying oligomeric state.
        product : Species
            The integrated membrane protein product after insertion.
        complex : Species, optional
            Pre-specified oligomeric complex. If None and size > 1,
            automatically creates a Complex of size monomers. Ignored for
            monomeric proteins (size = 1).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [integral_membrane_protein, product, complex]
            where complex is None for monomers or a Complex species for
            oligomers.

        Notes
        -----
        For monomeric proteins (size = 1), no oligomeric complex is formed
        and the complex element in the return list is None.

        For oligomeric proteins (size > 1), a complex containing 'size'
        copies of the monomer is created or used if provided.

        """
        if complex is None:
            size = integral_membrane_protein.size
            if size > 1:
                complex1 = Complex(
                    [integral_membrane_protein] * size,
                    compartment=integral_membrane_protein.compartment,
                )
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
        **kwargs,
    ):
        """Generate reactions for membrane protein integration.

        Creates reactions for oligomerization (if needed) and membrane
        integration. For oligomeric proteins, generates both oligomerization
        and integration reactions. For monomers, generates only the
        integration reaction.

        Parameters
        ----------
        integral_membrane_protein : Species
            The membrane protein monomer. Must have a size attribute.
        product : Species
            The integrated membrane protein product.
        complex : Species, optional
            Pre-specified oligomeric complex. If None and size > 1,
            automatically created.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            For oligomers (size > 1): List of two reactions
            [binding_rxn1, integration_rxn1].
            For monomers (size = 1): List of one reaction [integration_rxn1].

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme depends on oligomeric state:

        For oligomers (size > 1):

        1. size * monomer <--> oligomer (mass-action rates: 'kb_oligomer',
            'ku_oligomer')
        2. oligomer --> product (ProportionalHillNegative with k=kex,
            d=complex, K=kcat, n=4, s1=product)

        For monomers (size = 1):

        1. monomer --> product (ProportionalHillNegative with k=kex,
           d=integral_membrane_protein, K=kcat, n=4, s1=product)

        The integration reaction uses `ProportionalHillNegative` kinetics with
        Hill coefficient n=4 and inhibitor species s1=product to model
        saturation and product inhibition.

        """
        # Get Parameters
        kb_oligomer = component.get_parameter(
            'kb_oligomer', part_id=part_id, mechanism=self
        )
        ku_oligomer = component.get_parameter(
            'ku_oligomer', part_id=part_id, mechanism=self
        )
        kex = component.get_parameter(
            'kex', part_id=part_id, mechanism=self
        )
        kcat = component.get_parameter(
            'kcat', part_id=part_id, mechanism=self
        )

        size = integral_membrane_protein.size

        if complex is None:
            if size > 1:
                complex1 = Complex(
                    [integral_membrane_protein] * size,
                    compartment=integral_membrane_protein.compartment,
                )
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


class Diffusion_Facilitated_Channel(Mechanism):
    """Diffusion mechanism facilitated by a membrane channel.

    A 'diffusion' mechanism that models passive, bidirectional diffusion of
    substrates through a membrane channel. Unlike simple diffusion,
    this mechanism requires a membrane channel but does not consume
    energy. The channel acts catalytically, binding substrate and product
    but not being consumed.

    The reaction follows the schema:
    $$
        'substrate' + 'membrane_channel' <--> 'product' + 'membrane_channel'
    $$

    Parameters
    ----------
    name : str, default='diffusion_facilitated_channel'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='diffusion'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('diffusion').

    See Also
    --------
    Simple_Diffusion : Passive diffusion without proteins.
    Diffusion_Facilitated_Carrier : Transport with MM kinetics.
    Transport_PrimaryActive_ABCexporter : Energy-dependent active transport.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models diffusion through membrane channels such as
    ion channels, aquaporins, and other pore-forming proteins. The channel
    facilitates movement down concentration gradients without conformational
    changes or energy expenditure.

    Common examples include:

    - Ion channels (K+, Na+, Ca2+ channels)
    - Aquaporins for water transport
    - Gap junctions between cells
    - Porins in bacterial outer membranes

    The transport is bidirectional with equal forward and reverse rate
    constants, reflecting passive equilibration across the membrane.

    Required parameters for this mechanism:

    - 'k_diff' : Transport rate constant (same for both directions)

    Examples
    --------
    Model potassium transport through an ion channel:

    >>> protein = bcp.IntegralMembraneProtein(
    ...     membrane_protein='Knck1',
    ...     product='K_channel',
    ...     compartment='cytoplasm',
    ...     membrane_compartment='membrane',
    ... )
    >>> channel = bcp.MembraneChannel(
    ...     membrane_channel=protein,
    ...     substrate='K',
    ...     internal_compartment='cytoplasm',
    ...     external_compartment='external'
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[protein, channel],
    ...     mechanisms={
    ...         'membrane_integration': bcp.Integration_MembraneProtein(),
    ...         'diffusion': bcp.Diffusion_Facilitated_Channel(),
    ...     },
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='diffusion_facilitated_channel',
        mechanism_type='diffusion',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(self, membrane_channel, substrate, product, **kwargs):
        """Generate species for simple transport.

        Returns the membrane channel, substrate, and product species
        involved in the transport reaction.

        Parameters
        ----------
        membrane_channel : Species
            The membrane channel through which transport occurs.
        substrate : Species
            The substrate species being transported (typically intracellular
            side).
        product : Species
            The product species after transport (typically extracellular
            side).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [membrane_channel, substrate, product].

        """
        return [membrane_channel, substrate, product]

    def update_reactions(
        self,
        membrane_channel,
        substrate,
        product,
        component=None,
        part_id=None,
        k_diff=None,
        **kwargs,
    ):
        r"""Generate reactions for simple membrane protein transport.

        Creates a single reversible mass-action reaction representing
        passive transport through a membrane channel with equal forward and
        reverse rate constants. The channel acts catalytically and is not
        consumed.

        Parameters
        ----------
        membrane_channel : Species
            The membrane channel facilitating transport.
        substrate : Species
            The substrate species being transported.
        product : Species
            The product species after transport.
        component : Component, optional
            Component containing parameter values. Required if k_diff is
            not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is
            provided, defaults to component.name.
        k_diff : Parameter or float, optional
            Transport rate constant. If None, retrieved from component
            parameters. Used as both forward and reverse rate constant.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single reversible mass-action reaction for
            transport.

        Raises
        ------
        ValueError
            If component is None and k_diff is not provided.

        Notes
        -----
        The reaction has equal forward and reverse rate constants:
        $$
            'substrate' + 'membrane_channel' <--> 'product' +
                'membrane_channel' \quad ('rates': 'k_diff', 'k_diff')
        $$
        The membrane channel appears on both sides of the reaction,
        indicating it acts catalytically and is recycled.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_diff is None):
            raise ValueError(
                "Must pass in a Component or values for k_diff."
            )
        if k_diff is None:
            k_diff = component.get_parameter(
                'k_diff', part_id=part_id, mechanism=self
            )
        else:
            k_diff = k_diff

        # Simple membrane protein transport
        # Sub (Internal) <--> Product (External)
        Simple_Transport_rxn = Reaction.from_massaction(
            inputs=[substrate, membrane_channel],
            outputs=[product, membrane_channel],
            k_forward=k_diff,
            k_reverse=k_diff,
        )
        return [Simple_Transport_rxn]


class Diffusion_Facilitated_Carrier(Mechanism):
    r"""Facilitated diffusion mechanism facilitated by a membrane carrier.

    A 'diffusion' mechanism that models facilitated diffusion of substrates
    through membrane carrier proteins. Unlike simple channels, carriers
    undergo conformational changes to transport substrates across membranes.

    The mechanism follows a multi-step Michaelis-Menten kinetic scheme with
    the following reaction pathway:

    1. Extracellular substrate binding and unbinding:
    $$
        'sub\_out' + 'MC\_out' <--> 'sub:MC'
    $$

    2. Carrier conformational change (substrate translocation):
    $$
        'sub:MC' <--> 'prod:MC'
    $$

    3. Intracellular substrate release and binding:
    $$
        'prod:MC' <--> 'sub\_in' + 'MC\_in'
    $$

    4. Empty carrier conformational reset:
    $$
        'MC\_in' <--> 'MC\_out'
    $$

    where `MC_out` and `MC_in` represent the carrier protein in its
    outward-facing and inward-facing conformations, respectively, while
    `sub:MC` and `prod:MC` represent the substrate-bound carrier
    intermediates.

    Parameters
    ----------
    name : str, default='diffusion_facilitated_carrier'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='diffusion'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('diffusion').

    See Also
    --------
    Diffusion_Facilitated_Channel : Passive transport via membrane channels.
    Transport_PrimaryActive_ABCexporter : Energy-dependent active transport.
    MichaelisMenten : Enzyme mechanism with similar kinetics.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models facilitated diffusion by carrier proteins that
    alternate between substrate-bound and product-bound conformations. The
    carrier binds substrate on one side of the membrane, undergoes a
    conformational change to transport it across, releases it as product,
    and returns to the original conformation.

    Key characteristics:

    - Does not require ATP or other energy sources
    - Transport is driven by concentration gradients
    - Carrier proteins alternate between conformational states
    - Follows Michaelis-Menten saturation kinetics

    Common examples include:

    - GLUT transporters for glucose
    - Amino acid carriers
    - Nucleoside transporters
    - Urea transporters

    The binding steps use GeneralPropensity objects with logistic sigmoid
    functions to ensure proper directionality based on species concentration
    gradients.

    Required parameters for this mechanism:

    - 'kb_subMC' : Forward binding rate for outer substrate to carrier
        ('MC_out')
    - 'ku_subMC' : Unbinding rate for outer substrate from carrier complex
        ('sub:MC')
    - 'kf_trnspMC' : Forward conformational change rate (transport step)
    - 'kr_trnspMC' : Reverse conformational change rate
    - 'kb_prodMC' : Forward binding rate for inner substrate to carrier
        ('MC_in')
    - 'ku_prodMC' : Unbinding rate for inner substrate from carrier complex
        ('prod:MC')
    - 'k_out' : Rate of carrier resetting from inner to outer conformation
    - 'k_in' : Rate of carrier changing from outer to inner conformation

    Examples
    --------
    Model glucose transport through a GLUT transporter:

    >>> glc_in = bcp.Species('glucose', compartment='cytoplasm')
    >>> glc_out = bcp.Species('glucose', compartment='external')
    >>> carrier = bcp.MembraneCarrier(
    ...     membrane_carrier='GlucoseTransporter',
    ...     substrate=glc_out,
    ...     external_compartment='external',
    ...     internal_compartment='cytoplasm',
    ... )
    >>> mechanism = bcp.Diffusion_Facilitated_Carrier()
    >>> mixture = bcp.Mixture(
    ...     components=[carrier],
    ...     mechanisms={'diffusion': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='diffusion_facilitated_carrier',
        mechanism_type='diffusion',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        complex_dict=None,
        **kwargs,
    ):
        """Generate species for facilitated transport.

        Creates species for the outer carrier conformation ('out' attribute),
        inner carrier conformation ('in' attribute), substrate in/out species,
        and two intermediate complexes (`sub:MC` and `prod:MC`).

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein that facilitates transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species with keys 'sub:MC' and
            'prod:MC'. If None, complexes are automatically created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [membrane_carrier (out), carrier_in, substrate_in,
            substrate_out, sub:MC complex, prod:MC complex].

        Notes
        -----
        The method creates two complex species representing intermediates in
        the transport cycle:

        1. sub:MC : substrate_out:membrane_carrier complex
        2. prod:MC : substrate_in:carrier_in complex

        """
        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )
        membrane_carrier.add_attribute('out')

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['sub:MC'] = Complex(
                [substrate_out, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )
            # Complex2
            complex_dict['prod:MC'] = Complex(
                [substrate_in, carrier_in],
                compartment=membrane_carrier.compartment,
            )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return ([membrane_carrier, carrier_in, substrate_in, substrate_out]
                + complex_array)

    def update_reactions(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for facilitated transport.

        Creates six reactions representing the complete transport cycle:
        outer substrate binding, outer substrate unbinding, translocation
        conformational change, inner substrate unbinding, inner substrate
        binding, and empty carrier conformational reset.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein facilitating transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes
            are automatically created using the same logic as in
            update_species.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of six reactions: [binding_rxn1, unbinding_rxn1,
            transport_rxn, unbinding_rxn2, binding_rxn2, config_rxn].

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme follows this pathway:

        1. sub_out + MC_out --> sub:MC (GeneralPropensity with logistic
            function, rate: 'kb_subMC')
        2. sub:MC --> sub_out + MC_out (mass-action, rate: 'ku_subMC')
        3. sub:MC <--> prod:MC (reversible mass-action, rates: 'kf_trnspMC',
            'kr_trnspMC')
        4. prod:MC --> sub_in + MC_in (mass-action, rate: 'ku_prodMC')
        5. sub_in + MC_in --> prod:MC (GeneralPropensity with logistic
            function, rate: 'kb_prodMC')
        6. MC_in <--> MC_out (reversible mass-action, rates: 'k_out', 'k_in')

        The binding steps use a GeneralPropensity with a logistic sigmoid
        function to enforce concentration gradient-driven directionality.
        The logistic function provides a continuous approximation of a
        Heaviside step, ensuring binding preferentially occurs when the
        substrate concentration on the origin side exceeds that on the
        destination side.

        """
        # Get Parameters
        kb_subMC = component.get_parameter(
            'kb_subMC', part_id=part_id, mechanism=self
        )
        ku_subMC = component.get_parameter(
            'ku_subMC', part_id=part_id, mechanism=self
        )
        kf_trnspMC = component.get_parameter(
            'kf_trnspMC', part_id=part_id, mechanism=self
        )
        kr_trnspMC = component.get_parameter(
            'kr_trnspMC', part_id=part_id, mechanism=self
        )
        kb_prodMC = component.get_parameter(
            'kb_prodMC', part_id=part_id, mechanism=self
        )
        ku_prodMC = component.get_parameter(
            'ku_prodMC', part_id=part_id, mechanism=self
        )
        k_out = component.get_parameter(
            'k_out', part_id=part_id, mechanism=self
                )
        k_in = component.get_parameter(
            'k_in', part_id=part_id, mechanism=self
        )

        # Carrier in a differernt configuration
        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['sub:MC'] = Complex(
                [substrate_out, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )
            # Complex2
            complex_dict['prod:MC'] = Complex(
                [substrate_in, carrier_in],
                compartment=membrane_carrier.compartment,
            )

        # Facilitated membrane protein transport
        # Sub + MC --> Sub:MC
        prop_subMC = GeneralPropensity(
            f"kb_subMC * {substrate_out} * {membrane_carrier} *"
            f"(1 / (1 + exp(-10 * ({substrate_out} - {substrate_in}))))",
            propensity_species=[
                substrate_out, substrate_in, membrane_carrier
            ],
            propensity_parameters=[kb_subMC],
        )
        binding_rxn1 = Reaction(
            [substrate_out, membrane_carrier],
            [complex_dict['sub:MC']],
            propensity_type=prop_subMC,
        )

        # Sub:MC --> Sub + MC
        unbinding_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['sub:MC']],
            outputs=[membrane_carrier, substrate_out],
            k_forward=ku_subMC,
        )

        # Sub:MC <--> Prod:MC
        transport_rxn = Reaction.from_massaction(
            inputs=[complex_dict['sub:MC']],
            outputs=[complex_dict['prod:MC']],
            k_forward=kf_trnspMC, k_reverse=kr_trnspMC
        )

        # MC:Prod --> MC + Prod
        unbinding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict['prod:MC']],
            outputs=[substrate_in, carrier_in],
            k_forward=ku_prodMC,
        )
        prop_probMC = GeneralPropensity(
            f"kb_prodMC * {substrate_in} * {carrier_in} *"
            f" (1 / (1 + exp(-10 * ({substrate_in} - {substrate_out}))))",
            propensity_species=[substrate_in, substrate_out, carrier_in],
            propensity_parameters=[kb_prodMC],
        )

        # MC + Prod --> MC:Prod
        binding_rxn2 = Reaction(
            [substrate_in, carrier_in],
            [complex_dict['prod:MC']],
            propensity_type=prop_probMC,
        )

        # MC_in <--> MC
        config_rxn = Reaction.from_massaction(
            inputs=[carrier_in],
            outputs=[membrane_carrier],
            k_forward=k_out, k_reverse=k_in
        )

        return [binding_rxn1, unbinding_rxn1, transport_rxn,
                unbinding_rxn2, binding_rxn2, config_rxn]


class Transport_SecondaryActive_Symporter(Mechanism):
    r"""Secondary active transport mechanism enabled by a membrane carrier.

    A secondary active transport mechanism that models transport of substrates
    co-transported with driving ions down their electrochemical gradient
    across membranes.

    The mechanism follows a multi-step kinetic scheme with the following
    reaction pathway:

    1. Extracellular ion binding and unbinding:
    $$
        'ion\_out' + 'MC\_out' <--> 'ion\_out:MC'
    $$

    2. Extracellular substrate binding and unbinding:
    $$
        'sub\_out' + 'ion\_out:MC' <--> 'sub:ion:MC'
    $$

    3. Carrier conformational change (translocation):
    $$
        'sub:ion:MC' <--> 'prod:ion:MC'
    $$

    4. Intracellular substrate release and binding:
    $$
        'prod:ion:MC' <--> 'ion\_in:MC' + 'sub\_in'
    $$

    5. Intracellular ion release and binding:
    $$
        'ion\_in:MC' <--> 'MC\_in' + 'ion\_in'
    $$

    6. Empty carrier conformational reset:
    $$
        'MC\_in' <--> 'MC\_out'
    $$

    where `MC_out` and `MC_in` represent the carrier protein in its
    outward-facing and inward-facing conformations, respectively.

    Parameters
    ----------
    driving_ion : dict, optional
        Dictionary mapping driving ion names to stoichiometric ratios
        (e.g., `{'Na+': '3:1'}`).
    name : str, default='transport_secondaryactive_symporter'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    driving_ion : dict
        Dictionary mapping driving ions to stoichiometric ratio strings.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    Diffusion_Facilitated_Carrier : Passive carrier-mediated transport.
    Transport_PrimaryActive_ABCexporter : Energy-dependent active transport.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models secondary active symport where substrate transport
    is coupled to driving ion movement across the membrane. The binding steps
    for driving ions use GeneralPropensity objects with logistic sigmoid
    functions to ensure directionality driven by concentration gradients.

    Required parameters for this mechanism:

    - 'kb_ionMC_out' : Forward binding rate for extracellular driving ion to
        'MC_out'
    - 'ku_ionMC_out' : Unbinding rate for extracellular driving ion from
        'ion_out:MC'
    - 'kb_subMC' : Forward binding rate for extracellular substrate to
        'ion_out:MC'
    - 'ku_subMC' : Unbinding rate for extracellular substrate from
        'sub:ion:MC'
    - 'kf_trnspMC' : Forward conformational change rate
    - 'kr_trnspMC' : Reverse conformational change rate
    - 'ku_prodMC' : Unbinding rate for intracellular substrate from
        'prod:ion:MC'
    - 'kb_prodMC' : Binding rate for intracellular substrate to 'ion_in:MC'
    - 'ku_ionMC_in' : Unbinding rate for intracellular driving ion from
        'ion_in:MC'
    - 'kb_ionMC_in' : Forward binding rate for intracellular driving ion to
        'MC_in'
    - 'k_out' : Rate of carrier resetting from inner to outer conformation
    - 'k_in' : Rate of carrier changing from outer to inner conformation

    Examples
    --------
    Model Na+/glucose symport across a membrane:

    >>> sglt1_carrier = bcp.MembraneCarrier(
    ...     membrane_carrier='sglt1',
    ...     substrate='glucose',
    ...     internal_compartment='cytoplasm',
    ...     external_compartment='extracellular'
    ... )
    >>> driving_ions = {'Na0': '2:1', 'Na1': '1:1'}
    >>> mechanism = bcp.Transport_SecondaryActive_Symporter(
    ...     driving_ion=driving_ions)
    >>> mixture = bcp.Mixture(
    ...    "transport_GLLUT1-glucose",
    ...    components = [sglt1_carrier],
    ...    mechanisms = {'transport': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        driving_ion=None,
        name='transport_secondaryactive_symporter',
        mechanism_type='transport',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )
        self.driving_ion = driving_ion

    def update_species(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        driving_ion=None,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for secondary active symport.

        Creates species for the outer carrier conformation ('out' attribute),
        inner carrier conformation ('in' attribute), substrate in/out species,
        driving ion in/out species, and intermediate complexes formed during
        transport.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein that facilitates transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        driving_ion : dict, optional
            Dictionary mapping driving ion names to stoichiometric ratio
            strings (e.g., `{'Na+': '3:1'}`). If None, defaults to
            `self.driving_ion`.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes are
            automatically created.
        component : Component, optional
            Component containing parameters (unused in species generation).
        part_id : str, optional
            Identifier for parameter lookup (unused in species generation).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [membrane_carrier (out), carrier_in, substrate_in,
            substrate_out, ions_in, ions_out, intermediate complexes].

        Raises
        ------
        ValueError
            If `driving_ion` is not provided and `self.driving_ion` is None.

        Notes
        -----
        For each driving ion specified, four complex species are created
        representing intermediates in the transport cycle:

        1. key_out:MC : ion_out:membrane_carrier complex
        2. sub:key:MC : substrate_out:key_out:MC complex
        3. key_in:MC : ion_in:carrier_in complex
        4. prod:key:MC : substrate_in:key_in:MC complex

        """
        if driving_ion is None:
            driving_ion = self.driving_ion
        if driving_ion is None:
            raise ValueError("A driving_ion must be provided to the " \
            "mechanism.")

        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )
        membrane_carrier.add_attribute('out')

        # Create empty lists
        ions_in = []
        ions_out = []

        for ion in driving_ion.keys():
            ion_in = Species(ion,
                            compartment=substrate_in.compartment,
            )
            ion_out = Species(ion,
                            compartment=substrate_out.compartment,
            )
            ions_in.append(ion_in)
            ions_out.append(ion_out)

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            for key, ionI, ionO in zip(
                driving_ion.keys(), ions_in, ions_out
            ):
                # Parse dictionary for ratios
                ratio_parts = driving_ion[key].split(':')
                # Convert the string pieces into integers
                ion_num = int(ratio_parts[0])
                sub_num = int(ratio_parts[1])

                # Complex 1
                complex_dict[f'{key}_out:MC'] = Complex(
                    [ion_num*[ionO], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )
                # Complex2
                complex_dict[f'sub:{key}:MC'] = Complex(
                    [sub_num*[substrate_out], complex_dict[f'{key}_out:MC']],
                    compartment=membrane_carrier.compartment,
                )
                # Complex4
                complex_dict[f'{key}_in:MC'] = Complex(
                    [ion_num*[ionI], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex3
                complex_dict[f'prod:{key}:MC'] = Complex(
                    [sub_num*[substrate_in], complex_dict[f'{key}_in:MC']],
                    compartment=membrane_carrier.compartment,
                )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [membrane_carrier, carrier_in, substrate_in,
                substrate_out] + ions_in + ions_out + complex_array

    def update_reactions(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        driving_ion=None,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for secondary active symport.

        Creates reactions representing the complete symport cycle per driving
        ion: extracellular ion binding, extracellular ion unbinding, substrate
        binding/unbinding, translocation conformational change, intracellular
        substrate release/binding, intracellular ion release, intracellular
        ion binding, and empty carrier reset.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein facilitating transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        driving_ion : dict, optional
            Dictionary mapping driving ion names to stoichiometric ratio
            strings. If None, defaults to `self.driving_ion`.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes are
            automatically created using the same logic as in update_species.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of reactions generated for each driving ion plus the empty
            carrier  reset reaction (`config_rxn`).

        Raises
        ------
        ValueError
            If `driving_ion` is not provided and `self.driving_ion` is None.
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        For each driving ion in `driving_ion`, the reaction scheme follows
        this pathway:

        1. ion_out + MC_out --> key_out:MC
            (GeneralPropensity with logistic function, rate: 'kb_ionMC_out')
        2. key_out:MC --> ion_out + MC_out (mass-action, rate: 'ku_ionMC_out')
        3. sub_out + key_out:MC <--> sub:key:MC
            (reversible mass-action, rates: 'kb_subMC', 'ku_subMC')
        4. sub:key:MC <--> prod:key:MC
            (reversible mass-action, rates: 'kf_trnspMC', 'kr_trnspMC')
        5. prod:key:MC <--> key_in:MC + sub_in
            (reversible mass-action, rates: 'ku_prodMC', 'kb_prodMC')
        6. key_in:MC --> ion_in + MC_in
            (mass-action, rate: 'ku_ionMC_in')
        7. ion_in + MC_in --> key_in:MC
            (GeneralPropensity with logistic function, rate: 'kb_ionMC_in')

        Additionally, a final empty carrier reset reaction is generated:
        8. MC_in <--> MC_out (reversible mass-action, rates: 'k_out', 'k_in')

        Ion binding steps use GeneralPropensity with logistic sigmoid
        functions to enforce concentration gradient-driven
        directionality.

        """
        if driving_ion is None:
            driving_ion = self.driving_ion
        if driving_ion is None:
            raise ValueError("A driving_ion must be provided to the" \
            "mechanism.")

        # Get Parameters
        kb_ionMC_out = component.get_parameter(
            'kb_ionMC_out', part_id=part_id, mechanism=self
        )
        ku_ionMC_out = component.get_parameter(
            'ku_ionMC_out', part_id=part_id, mechanism=self
        )
        kb_subMC = component.get_parameter(
            'kb_subMC', part_id=part_id, mechanism=self
        )
        ku_subMC = component.get_parameter(
            'ku_subMC', part_id=part_id, mechanism=self
        )
        kf_trnspMC = component.get_parameter(
            'kf_trnspMC', part_id=part_id, mechanism=self
        )
        kr_trnspMC = component.get_parameter(
            'kr_trnspMC', part_id=part_id, mechanism=self
        )
        kb_prodMC = component.get_parameter(
            'kb_prodMC', part_id=part_id, mechanism=self
        )
        ku_prodMC = component.get_parameter(
            'ku_prodMC', part_id=part_id, mechanism=self
        )
        kb_ionMC_in = component.get_parameter(
            'kb_ionMC_in', part_id=part_id, mechanism=self
        )
        ku_ionMC_in = component.get_parameter(
            'ku_ionMC_in', part_id=part_id, mechanism=self
        )
        k_out = component.get_parameter(
            'k_out', part_id=part_id, mechanism=self
                )
        k_in = component.get_parameter(
            'k_in', part_id=part_id, mechanism=self
        )

        # Carrier in a differernt configuration
        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )
        membrane_carrier.add_attribute('out')

        # Create empty list
        ions_in = []
        ions_out = []

        for ion in driving_ion.keys():
            ion_in = Species(ion,
                            compartment=substrate_in.compartment,
            )
            ion_out = Species(ion,
                            compartment=substrate_out.compartment,
            )
            ions_in.append(ion_in)
            ions_out.append(ion_out)

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            for key, ionI, ionO in zip(
                driving_ion.keys(), ions_in, ions_out
            ):
                # Parse dictionary for ratios
                ratio_parts = driving_ion[key].split(':')
                # Convert the string pieces into integers
                ion_num = int(ratio_parts[0])
                sub_num = int(ratio_parts[1])

                # Complex 1
                complex_dict[f'{key}_out:MC'] = Complex(
                    [ion_num*[ionO], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )
                # Complex2
                complex_dict[f'sub:{key}:MC'] = Complex(
                    [sub_num*[substrate_out], complex_dict[f'{key}_out:MC']],
                    compartment=membrane_carrier.compartment,
                )
                # Complex4
                complex_dict[f'{key}_in:MC'] = Complex(
                    [ion_num*[ionI], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex3
                complex_dict[f'prod:{key}:MC'] = Complex(
                    [sub_num*[substrate_in], complex_dict[f'{key}_in:MC']],
                    compartment=membrane_carrier.compartment,
                )

        # Secondary Active Transport
        secondaryactive_rxns = [] # Create an empty list

        for key, ionI, ionO in zip(
            driving_ion.keys(), ions_in, ions_out
        ):
            # Parse dictionary for ratios
            ratio_parts = driving_ion[key].split(':')
            # Convert the string pieces into integers
            ion_num = int(ratio_parts[0])
            sub_num = int(ratio_parts[1])

            # ion_out + MC --> ion:MC
            prop_ionMC_out = GeneralPropensity(
                f"kb_ionMC_out * {ionO} * {membrane_carrier} *"
                f"(1 / (1 + exp(-10 * ({ionO} - {ionI}))))",
                propensity_species=[
                    ionO, ionI, membrane_carrier
                ],
                propensity_parameters=[kb_ionMC_out],
            )
            rxn_binding_ionMC = Reaction(
                        [ion_num*[ionO], membrane_carrier],
                        [complex_dict[f'{key}_out:MC']],
                        propensity_type=prop_ionMC_out,
                    )
            secondaryactive_rxns.append(rxn_binding_ionMC)

            # ion_out + MC <-- ion:MC
            rxn_unbinding_ionMC = Reaction.from_massaction(
                [complex_dict[f'{key}_out:MC']],
                [ion_num*[ionO], membrane_carrier],
                k_forward=ku_ionMC_out,
            )
            secondaryactive_rxns.append(rxn_unbinding_ionMC)

            # sub_out + ion:MC <--> sub_out:ion:MC
            rxn_subMC = Reaction.from_massaction(
                [sub_num*[substrate_out], complex_dict[f'{key}_out:MC']],
                [complex_dict[f'sub:{key}:MC']],
                k_forward=kb_subMC,
                k_reverse=ku_subMC,
            )
            secondaryactive_rxns.append(rxn_subMC)

            # sub_out:ion:MC <--> sub_in:ion:MC
            rxn_transport = Reaction.from_massaction(
                [complex_dict[f'sub:{key}:MC']],
                [complex_dict[f'prod:{key}:MC']],
                k_forward=kf_trnspMC,
                k_reverse=kr_trnspMC,
            )
            secondaryactive_rxns.append(rxn_transport)

            # sub_in:ion:MC <--> sub_in + ion:MC
            rxn_subrelease = Reaction.from_massaction(
                [complex_dict[f'prod:{key}:MC']],
                [complex_dict[f'{key}_in:MC'], sub_num*[substrate_in]],
                k_forward=ku_prodMC,
                k_reverse=kb_prodMC,
            )
            secondaryactive_rxns.append(rxn_subrelease)

            # ion:MC --> ion_in + MC
            rxn_ionrelease = Reaction.from_massaction(
                [complex_dict[f'{key}_in:MC']],
                 [ion_num*[ionI], carrier_in],
                k_forward=ku_ionMC_in,
            )
            secondaryactive_rxns.append(rxn_ionrelease)

            # ion:MC <-- ion_in + MC
            prop_ionMC_in = GeneralPropensity(
                f"kb_ionMC_in * {ionI} * {carrier_in} *"
                f" (1 / (1 + exp(-10 * ({ionI} - {ionO}))))",
                propensity_species=[ionI, ionO, carrier_in],
                propensity_parameters=[kb_ionMC_in],
            )
            rxn_binding_ionMC2 = Reaction(
                [ion_num*[ionI], carrier_in],
                [complex_dict[f'{key}_in:MC']],
                propensity_type=prop_ionMC_in,
            )
            secondaryactive_rxns.append(rxn_binding_ionMC2)

        # MC_in <--> MC
        config_rxn = Reaction.from_massaction(
            inputs=[carrier_in],
            outputs=[membrane_carrier],
            k_forward=k_out, k_reverse=k_in
        )

        return secondaryactive_rxns + [config_rxn]


class Transport_SecondaryActive_Antiporter(Mechanism):
    r"""Secondary active transport mechanism enabled by a membrane carrier.

    A 'transport' mechanism that models secondary active antiport transport,
    where the translocation of a substrate across a membrane is coupled to the
    co-transport of a driving ion in the opposite direction.

    The mechanism follows an 8-step kinetic scheme per driving ion:

    1. Extracellular driving ion binding to outward carrier:
    $$
        n \text{ion\_out} + \text{MC\_out} \longrightarrow \text{key\_out:MC}
    $$

    2. Extracellular driving ion unbinding:
    $$
        \text{key\_out:MC} \longrightarrow n \text{ion\_out} + \text{MC\_out}
    $$

    3. Ion-bound carrier conformational change:
    $$
        \text{key\_out:MC} \longleftrightarrow \text{key\_in:MC}
    $$

    4. Intracellular driving ion unbinding from inward carrier:
    $$
        \text{key\_in:MC} \longrightarrow n \text{ion\_in} + \text{MC\_in}
    $$

    5. Intracellular driving ion binding:
    $$
        n \text{ion\_in} + \text{MC\_in} \longrightarrow
        \text{key\_in:MC}
    $$

    6. Intracellular substrate binding and unbinding:
    $$
        n \text{sub\_in} + \text{MC\_in} \longleftrightarrow \text{sub:MC}
    $$

    7. Substrate-bound carrier conformational change (translocation):
    $$
        \text{sub:MC} \longleftrightarrow \text{prod:MC}
    $$

    8. Extracellular substrate release and binding:
    $$
        \text{prod:MC} \longleftrightarrow n
        \text{sub\_out} + \text{MC\_out}
    $$

    where `MC_out` and `MC_in` represent the carrier protein in outward- and
    inward-facing states, respectively. `key_out:MC` and `key_in:MC` are
    ion-bound carrier intermediates, and `sub:MC` and `prod:MC` are
    substrate-bound carrier intermediates.

    Parameters
    ----------
    driving_ion : dict, optional
        Dictionary mapping driving ion species names to stoichiometric ratios
        in string format `'ion_ratio:sub_ratio'` (e.g., `{'Na': '3:1'}`).
    name : str, default='transport_secondaryactive_antiporter'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    driving_ion : dict
        Dictionary of driving ion stoichiometry configurations.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    Diffusion_Facilitated_Carrier : Passive carrier-mediated transport.
    Transport_PrimaryActive_ABCexporter : Energy-dependent active transport.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    Required parameters for this mechanism:

    - 'kb_ionMC_out' : Outer driving ion forward binding rate
    - 'ku_ionMC_out' : Outer driving ion unbinding rate
    - 'kf_ionX'      : Forward rate of ion-bound carrier conformational flip
    - 'kr_ionX'      : Reverse rate of ion-bound carrier conformational flip
    - 'kb_subMC'     : Inner substrate binding rate
    - 'ku_subMC'     : Inner substrate unbinding rate
    - 'kf_trnspMC'   : Forward rate of substrate-bound carrier translocation
    - 'kr_trnspMC'   : Reverse rate of substrate-bound carrier translocation
    - 'kb_prodMC'    : Outer substrate reverse binding rate
    - 'ku_prodMC'    : Outer substrate release rate
    - 'kb_ionMC_in'  : Inner driving ion reverse binding rate
    - 'ku_ionMC_in'  : Inner driving ion release rate

    Examples
    --------
    Model a Na+/Ca2+ antiporter mechanism:

    >>> ncx_carrier = bcp.MembraneCarrier(
    ...     membrane_carrier='ncx',
    ...     substrate='Ca2',
    ...     internal_compartment='cytoplasm',
    ...     external_compartment='extracellular'
    ... )
    >>> driving_ions = {'Na': '3:1'}
    >>> mechanism = bcp.Transport_SecondaryActive_Antiporter(
    ...     driving_ion=driving_ions)
    >>> mixture = bcp.Mixture(
    ...     "transport_GLLUT1-glucose",
    ...     components=[ncx_carrier],
    ...     mechanisms={'transport': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        driving_ion=None,
        name='transport_secondaryactive_antiporter',
        mechanism_type='transport',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )
        self.driving_ion = driving_ion

    def update_species(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        driving_ion=None,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for secondary active antiport transport.

        Creates species for outward membrane carrier (`MC_out`), inward
        membrane carrier (`MC_in`), substrate species in both compartments
        (`substrate_in`, `substrate_out`), driving ion species in both
        compartments (`ions_in`, `ions_out`), and four complex species per
        driving ion (`key_out:MC`, `key_in:MC`, `sub:MC`, `prod:MC`).

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein facilitating transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        driving_ion : dict, optional
            Dictionary mapping driving ion names to stoichiometry ratios
            `'ion:sub'` (e.g., `{'Na': '3:1'}`). If None, uses
            `self.driving_ion`.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes are
            automatically created.
        component : Component, optional
            Component containing parameters (unused in update_species).
        part_id : str, optional
            Identifier for parameter lookup (unused in update_species).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing `[membrane_carrier, carrier_in, substrate_in,
            substrate_out]` followed by `ions_in`, `ions_out`, and all
            generated complex species.

        Raises
        ------
        ValueError
            If no `driving_ion` dictionary is provided or set on the instance.

        Notes
        -----
        The method creates four intermediate complex species per driving ion
        entry:

        1. key_out:MC : Driving ion bound to outer-facing carrier
            (`ion_out:MC`)
        2. key_in:MC  : Driving ion bound to inner-facing carrier
            (`ion_in:MC`)
        3. sub:MC     : Intracellular substrate bound to inner-facing carrier
            (`sub_in:MC`)
        4. prod:MC    : Extracellular substrate bound to outer-facing carrier
            (`sub_out:MC`)

        """
        if driving_ion is None:
            driving_ion = self.driving_ion
        if driving_ion is None:
            raise ValueError("A driving_ion must be provided to the " \
            "mechanism.")

        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )
        membrane_carrier.add_attribute('out')

        # Create empty lists
        ions_in = []
        ions_out = []

        for ion in driving_ion.keys():
            ion_in = Species(ion,
                            compartment=substrate_in.compartment,
            )
            ion_out = Species(ion,
                            compartment=substrate_out.compartment,
            )
            ions_in.append(ion_in)
            ions_out.append(ion_out)

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            for key, ionI, ionO in zip(
                driving_ion.keys(), ions_in, ions_out
            ):
                # Parse dictionary for ratios
                ratio_parts = driving_ion[key].split(':')
                # Convert the string pieces into integers
                ion_num = int(ratio_parts[0])
                sub_num = int(ratio_parts[1])

                # Complex1
                complex_dict[f'{key}_out:MC'] = Complex(
                    [ion_num*[ionO], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )
                # Complex2
                complex_dict[f'{key}_in:MC'] = Complex(
                    [ion_num*[ionI], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex3
                complex_dict['sub:MC'] = Complex(
                    [sub_num*[substrate_in], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex4
                complex_dict['prod:MC'] = Complex(
                    [sub_num*[substrate_out], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [membrane_carrier, carrier_in, substrate_in,
                substrate_out] + ions_in + ions_out + complex_array

    def update_reactions(
        self,
        membrane_carrier,
        substrate_in,
        substrate_out,
        driving_ion=None,
        complex_dict=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for secondary active antiport transport.

        Creates eight reactions per driving ion representing the complete
        antiport transport cycle: outer ion binding/unbinding, ion-carrier
        translocation, inner ion unbinding/binding, inner substrate
        binding/unbinding, substrate-carrier translocation, and outer
        substrate release/binding.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein facilitating transport.
        substrate_in : Species
            The intracellular substrate species.
        substrate_out : Species
            The extracellular substrate species.
        driving_ion : dict, optional
            Dictionary mapping driving ion names to stoichiometry ratios
            `'ion:sub'`. If None, uses `self.driving_ion`.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes
            are automatically created using the same logic as in
            `update_species`.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of eight reactions per driving ion:
            `[rxn_binding_ionMC, rxn_unbinding_ionMC, rxn_ionx,
            rxn_ionMC_inF,` `rxn_ionMC_inR, rxn_subMC, rxn_transport,
            rxn_prodMC]`.

        Raises
        ------
        ValueError
            If `driving_ion` is None.
        AttributeError
            If `component` or `part_id` is None (required for parameter
            lookup).

        Notes
        -----
        The reaction scheme per driving ion follows this pathway:

        1. ion_out + MC_out --> key_out:MC
            (GeneralPropensity with logistic function, rate: 'kb_ionMC_out')
        2. key_out:MC --> ion_out + MC_out (mass-action, rate: 'ku_ionMC_out')
        3. key_out:MC <--> key_in:MC
            (reversible mass-action, rates: 'kf_ionX', 'kr_ionX')
        4. key_in:MC --> ion_in + MC_in (mass-action, rate: 'ku_ionMC_in')
        5. ion_in + MC_in --> key_in:MC
            (GeneralPropensity with logistic function, rate: 'kb_ionMC_in')
        6. sub_in + MC_in <--> sub:MC
            (reversible mass-action, rates: 'kb_subMC', 'ku_subMC')
        7. sub:MC <--> prod:MC
            (reversible mass-action, rates: 'kf_trnspMC', 'kr_trnspMC')
        8. prod:MC <--> sub_out + MC_out
            (reversible mass-action, rates: 'ku_prodMC', 'kb_prodMC')

        The ion binding steps (1 and 5) use `GeneralPropensity` objects with
        logistic sigmoid functions to enforce continuous concentration
        gradient-driven directionality.

        """
        if driving_ion is None:
            driving_ion = self.driving_ion
        if driving_ion is None:
            raise ValueError("A driving_ion must be provided to the " \
            "mechanism.")

        # Get Parameters
        kb_ionMC_out = component.get_parameter(
            'kb_ionMC_out', part_id=part_id, mechanism=self
        )
        ku_ionMC_out = component.get_parameter(
            'ku_ionMC_out', part_id=part_id, mechanism=self
        )
        kf_ionX = component.get_parameter(
            'kf_ionX', part_id=part_id, mechanism=self
        )
        kr_ionX = component.get_parameter(
            'kr_ionX', part_id=part_id, mechanism=self
        )
        kb_subMC = component.get_parameter(
            'kb_subMC', part_id=part_id, mechanism=self
        )
        ku_subMC = component.get_parameter(
            'ku_subMC', part_id=part_id, mechanism=self
        )
        kf_trnspMC = component.get_parameter(
            'kf_trnspMC', part_id=part_id, mechanism=self
        )
        kr_trnspMC = component.get_parameter(
            'kr_trnspMC', part_id=part_id, mechanism=self
        )
        kb_prodMC = component.get_parameter(
            'kb_prodMC', part_id=part_id, mechanism=self
        )
        ku_prodMC = component.get_parameter(
            'ku_prodMC', part_id=part_id, mechanism=self
        )
        kb_ionMC_in = component.get_parameter(
            'kb_ionMC_in', part_id=part_id, mechanism=self
        )
        ku_ionMC_in = component.get_parameter(
            'ku_ionMC_in', part_id=part_id, mechanism=self
        )

        # Carrier in a differernt configuration
        carrier_in = Species(
            membrane_carrier.name,
            material_type='protein',
            compartment=membrane_carrier.compartment,
            attributes=['in']
        )
        membrane_carrier.add_attribute('out')

        # Create empty list
        ions_in = []
        ions_out = []

        for ion in driving_ion.keys():
            ion_in = Species(ion,
                            compartment=substrate_in.compartment,
            )
            ion_out = Species(ion,
                            compartment=substrate_out.compartment,
            )
            ions_in.append(ion_in)
            ions_out.append(ion_out)

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            for key, ionI, ionO in zip(
                driving_ion.keys(), ions_in, ions_out
            ):
                # Parse dictionary for ratios
                ratio_parts = driving_ion[key].split(':')
                # Convert the string pieces into integers
                ion_num = int(ratio_parts[0])
                sub_num = int(ratio_parts[1])

                # Complex1
                complex_dict[f'{key}_out:MC'] = Complex(
                    [ion_num*[ionO], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )
                # Complex2
                complex_dict[f'{key}_in:MC'] = Complex(
                    [ion_num*[ionI], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex3
                complex_dict['sub:MC'] = Complex(
                    [sub_num*[substrate_in], carrier_in],
                    compartment=membrane_carrier.compartment,
                )
                # Complex4
                complex_dict['prod:MC'] = Complex(
                    [sub_num*[substrate_out], membrane_carrier],
                    compartment=membrane_carrier.compartment,
                )

        # Secondary Active Transport
        secondaryactive_rxns = [] # Create an empty list

        for key, ionI, ionO in zip(
            driving_ion.keys(), ions_in, ions_out
        ):
            # Parse dictionary for ratios
            ratio_parts = driving_ion[key].split(':')
            # Convert the string pieces into integers
            ion_num = int(ratio_parts[0])
            sub_num = int(ratio_parts[1])

            # ion_out + MC <--> ion_out:MC
            prop_ionMC_out = GeneralPropensity(
                f"kb_ionMC_out * {ionO} * {membrane_carrier} *"
                f"(1 / (1 + exp(-10 * ({ionO} - {ionI}))))",
                propensity_species=[
                    ionO, ionI, membrane_carrier
                ],
                propensity_parameters=[kb_ionMC_out],
            )

            rxn_binding_ionMC = Reaction(
                [ion_num*[ionO], membrane_carrier],
                [complex_dict[f'{key}_out:MC']],
                propensity_type=prop_ionMC_out,
            )
            secondaryactive_rxns.append(rxn_binding_ionMC)

            # ion_out + MC <-- ion_out:MC
            rxn_unbinding_ionMC = Reaction.from_massaction(
                [complex_dict[f'{key}_out:MC']],
                [ion_num*[ionO], membrane_carrier],
                k_forward=ku_ionMC_out,
            )
            secondaryactive_rxns.append(rxn_unbinding_ionMC)

            # ion_out:MC_out <--> ion_in:MC_in
            rxn_ionx = Reaction.from_massaction(
                [complex_dict[f'{key}_out:MC']],
                [complex_dict[f'{key}_in:MC']],
                k_forward=kf_ionX,
                k_reverse=kr_ionX,
            )
            secondaryactive_rxns.append(rxn_ionx)

            # ion_in:MC_in --> ion_in + MC_in
            rxn_ionMC_inF = Reaction.from_massaction(
                [complex_dict[f'{key}_in:MC']],
                [ion_num*[ionI], carrier_in],
                k_forward=ku_ionMC_in,
            )
            secondaryactive_rxns.append(rxn_ionMC_inF)
            # ion_in:MC_in <-- ion_in + MC_in
            prop_ionMC_in = GeneralPropensity(
                f"kb_ionMC_in * {ionI} * {carrier_in} *"
                f" (1 / (1 + exp(-10 * ({ionI} - {ionO}))))",
                propensity_species=[ionI, ionO, carrier_in],
                propensity_parameters=[kb_ionMC_in],
            )
            rxn_ionMC_inR = Reaction(
                [ion_num*[ionI], carrier_in],
                [complex_dict[f'{key}_in:MC']],
                propensity_type=prop_ionMC_in,
            )
            secondaryactive_rxns.append(rxn_ionMC_inR)

            # sub_in + MC_in <--> sub_in:MC_in
            rxn_subMC = Reaction.from_massaction(
                [sub_num*[substrate_in], carrier_in],
                [complex_dict['sub:MC']],
                k_forward=kb_subMC,
                k_reverse=ku_subMC,
            )
            secondaryactive_rxns.append(rxn_subMC)

            # sub_in:MC_in <--> sub_out:MC_out
            rxn_transport = Reaction.from_massaction(
                [complex_dict['sub:MC']],
                [complex_dict['prod:MC']],
                k_forward=kf_trnspMC,
                k_reverse=kr_trnspMC,
            )
            secondaryactive_rxns.append(rxn_transport)

            # sub_out:MC_out <--> sub_out + MC_out
            rxn_prodMC = Reaction.from_massaction(
                [complex_dict['prod:MC']],
                [sub_num*[substrate_out], membrane_carrier],
                k_forward=ku_prodMC,
                k_reverse=kb_prodMC,
            )
            secondaryactive_rxns.append(rxn_prodMC)

        return secondaryactive_rxns


class Transport_PrimaryActive_ABCexporter(Mechanism):
    r"""Primary active transport mechanism with ATP-dependent pumping.

    A 'transport' mechanism that models primary active transport where
    substrates are moved against their concentration gradients using energy
    from ATP hydrolysis. The mechanism follows Michaelis-Menten kinetics
    with explicit binding, ATP hydrolysis, conformational change, and
    product release steps.

    The reaction pathway follows this scheme:

    1. Substrate binding and unbinding:
    $$
        'Sub' + 'MP' <--> 'MP:sub'
    $$

    2. ATP binding and unbinding:
    $$
        'MP:sub' + n'E' <--> 'MP:sub:ATP'
    $$

    3. Pump conformational change (substrate translocation):
    $$
        'MP:sub:ATP' <--> 'MP:prod:ATP'
    $$

    4. Product release and binding:
    $$
        'MP:prod:ATP' <--> 'MP:ATP' + 'Prod'
    $$

    5. ATP hydrolysis step:
    $$
        'MP:ATP' --> 'MP:ADP'
    $$

    6. ADP release and binding (empty pump reset):
    $$
        'MP:ADP' <--> 'MP' + n'W'
    $$

    where `MP` represents the membrane pump, `E` represents ATP (energy),
    and `W` represents ADP (waste).

    Parameters
    ----------
    name : str, default='transport_primaryactive_abcexporter'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/transport_parameters.tsv'
        Path to file containing default parameter values for this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    Diffusion_Facilitated_Carrier : Passive facilitated diffusion.
    Diffusion_Facilitated_Channel : Passive transport via membrane channels.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models primary active transporters such as P-type ATPases
    (e.g., Na+/K+-ATPase, Ca2+-ATPase), ABC transporters, and other pumps
    that directly couple ATP hydrolysis to substrate transport. The pump
    undergoes conformational changes driven by ATP binding and hydrolysis to
    move substrates against concentration gradients.

    Key characteristics:

    - Requires ATP or other energy source
    - Can transport substrates against concentration gradients
    - Undergoes ATP-dependent conformational changes
    - Follows Michaelis-Menten saturation kinetics

    Common examples include:

    - Na+/K+-ATPase (maintains ion gradients in animal cells)
    - Ca2+-ATPase (SERCA pump in muscle cells)
    - H+-ATPases (proton pumps in various organisms)
    - ABC transporters (drug efflux pumps)

    The mechanism requires the membrane pump to have an ATP attribute
    (membrane_pump.ATP) that specifies the number of ATP molecules required
    per transport cycle.

    The binding steps use GeneralPropensity objects with Heaviside step
    functions to ensure proper directionality based on species availability.

    Required parameters for this mechanism:

    - 'kb_subMP' : Forward binding rate for substrate to membrane pump ('MP')
    - 'ku_subMP' : Unbinding rate for substrate from pump complex ('MP:sub')
    - 'kb_subMPnATP' : Forward binding rate for ATP to substrate:pump complex
    - 'ku_subMPnATP' : Unbinding rate for ATP from substrate:pump complex
    - 'kf_trnspMP' : Forward conformational change rate (transport step)
    - 'kr_trnspMP' : Reverse conformational change rate
    - 'ku_prodMP' : Unbinding rate for product from pump complex
        ('MP:prod:ATP')
    - 'kb_prodMP' : Binding rate for product to pump complex ('MP:ATP')
    - 'ku_MP' : Unbinding/hydrolysis rate for ATP/ADP steps ('MP:ATP',
        'MP:ADP')
    - 'kb_MP' : Binding rate for ADP to membrane pump

    Examples
    --------
    Model active sodium transport by Na+/K+-ATPase:

    >>> pump = bcp.MembranePump(
    ...     membrane_pump='NaK_ATPase',
    ...     substrate='Na',
    ...     direction='exporter',
    ...     ATP=1
    ... )
    >>> mechanism = bcp.Transport_PrimaryActive_ABCexporter()
    >>> mixture = bcp.Mixture(
    ...     components=[pump],
    ...     mechanisms={'transport': mechanism},
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='transport_primaryactive_abcexporter',
        mechanism_type='transport',
        parameter_file='mechanisms/transport_parameters.tsv',
        **kwargs,
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        membrane_pump,
        substrate,
        product,
        energy,
        waste,
        complex_dict=None,
        **kwargs,
    ):
        """Generate species for primary active transport.

        Creates species for the membrane pump, substrate, product, ATP/ADP
        energy species, and all intermediate complexes formed during the
        ATP-driven transport cycle.

        Parameters
        ----------
        membrane_pump : Species
            The membrane pump protein that transports substrates using ATP.
            Must have an ATP attribute specifying the number of ATP
            molecules required per transport cycle.
        substrate : Species
            The substrate species being transported (typically intracellular
            side).
        product : Species
            The product species after transport (typically extracellular
            side). Usually the same molecular species as substrate but in a
            different compartment.
        energy : Species
            ATP species used to drive active transport.
        waste : Species
            ADP species produced after ATP hydrolysis.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species with keys 'MP:sub',
            'MP:sub:ATP', 'MP:prod:ATP', 'MP:ATP', and 'MP:ADP'. If None,
            complexes are automatically created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [membrane_pump, substrate, product, energy,
            waste, complex_array] where complex_array is a list of five
            Complex species generated.

        Notes
        -----
        The method creates five complex species representing intermediates
        in the active transport cycle:

        1. MP:sub : membrane_pump:substrate complex
        2. MP:sub:ATP : membrane_pump:substrate:nATP complex
        3. MP:prod:ATP : membrane_pump:product:nATP complex
        4. MP:ATP : membrane_pump:nATP complex
        5. MP:ADP : membrane_pump:nADP complex

        The number of ATP/ADP molecules (nATP) is determined by the
        membrane_pump.ATP attribute.

        """
        nATP = membrane_pump.ATP

        if 'exporter' not in membrane_pump.attributes:
            warnings.warn(
                "This mechanism is defined as an exporter, but is" \
                "currently being used as an importer.",
            )

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['MP:sub'] = Complex(
                [substrate, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex2
            complex_dict['MP:sub:ATP'] = Complex(
                [nATP * [energy], complex_dict['MP:sub']],
                compartment=membrane_pump.compartment,
            )
            # Complex3
            complex_dict['MP:prod:ATP'] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex4
            complex_dict['MP:ATP'] = Complex(
                [nATP * [energy], membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex5
            complex_dict['MP:ADP'] = Complex(
                [nATP * [waste], membrane_pump],
                compartment=membrane_pump.compartment,
            )

        # Make dictionary into array
        complex_array = [value for value in complex_dict.values()]

        return [
            membrane_pump,
            substrate,
            product,
            energy,
            waste] + complex_array

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
        **kwargs,
    ):
        """Generate reactions for primary active transport.

        Creates eight reactions representing the complete ATP-driven
        transport cycle: substrate binding, substrate unbinding, ATP
        binding, ATP unbinding, transport step, product release,
        product binding, ATP hydrolysis, and ADP release/reset.

        Parameters
        ----------
        membrane_pump : Species
            The membrane pump protein. Must have an ATP attribute.
        substrate : Species
            The substrate species being transported.
        product : Species
            The product species after transport.
        energy : Species
            ATP species used for active transport.
        waste : Species
            ADP species produced after ATP hydrolysis.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species. If None, complexes
            are automatically created using the same logic as in
            update_species.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of eight reactions: [rxn1_SubBinding, rxn1_SubUnbinding,
            rxn2_ATPbinding, rxn2_ATPunbinding, rxn3_transport,
            rxn4_SubRelease, rxn5_atpHydro, rxn6_reset].

        Raises
        ------
        AttributeError
            If `component` or `part_id` is None (required for parameter
            lookup).

        Notes
        -----
        The reaction scheme follows this pathway:

        1. Sub + MP --> MP:sub (`GeneralPropensity` with Heaviside using
           rate 'kb_subMP')
        2. MP:sub --> Sub + MP (mass-action, rate: 'ku_subMP')
        3. MP:sub + nATP --> MP:sub:ATP (`GeneralPropensity` with Heaviside
           using rate 'kb_subMPnATP')
        4. MP:sub:ATP --> MP:sub + nATP (mass-action, rate: 'ku_subMPnATP')
        5. MP:sub:ATP <--> MP:prod:ATP (reversible mass-action, rates:
           'kf_trnspMP', 'kr_trnspMP')
        6. MP:prod:ATP <--> MP:ATP + Prod (reversible mass-action, rates:
           'ku_prodMP', 'kb_prodMP')
        7. MP:ATP --> MP:ADP (mass-action, rate: 'ku_MP')
        8. MP:ADP <--> MP + nADP (reversible mass-action, rates:
           'ku_MP', 'kb_MP')

        The binding steps use `GeneralPropensity` with Heaviside functions to
        enforce proper directionality. The Heaviside functions ensure that
        reactions only proceed when the required species are present.

        The number of ATP/ADP molecules (nATP) is determined by
        membrane_pump.ATP attribute.

        """
        # Get Parameters
        kb_subMP = component.get_parameter(
            'kb_subMP', part_id=part_id, mechanism=self
        )
        ku_subMP = component.get_parameter(
            'ku_subMP', part_id=part_id, mechanism=self
        )
        kb_subMPnATP = component.get_parameter(
            'kb_subMPnATP', part_id=part_id, mechanism=self
        )
        ku_subMPnATP = component.get_parameter(
            'ku_subMPnATP', part_id=part_id, mechanism=self
        )
        kf_trnspMP = component.get_parameter(
            'kf_trnspMP', part_id=part_id, mechanism=self
        )
        kr_trnspMP = component.get_parameter(
            'kr_trnspMP', part_id=part_id, mechanism=self
        )
        ku_prodMP = component.get_parameter(
            'ku_prodMP', part_id=part_id, mechanism=self
        )
        kb_prodMP = component.get_parameter(
            'kb_prodMP', part_id=part_id, mechanism=self
        )
        ku_MP = component.get_parameter(
            'ku_MP', part_id=part_id, mechanism=self
        )
        kb_MP = component.get_parameter(
            'kb_MP', part_id=part_id, mechanism=self
        )

        nATP = membrane_pump.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}

            # Complex1
            complex_dict['MP:sub'] = Complex(
                [substrate, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            complex1 = complex_dict['MP:sub']

            # Complex2
            complex_dict['MP:sub:ATP'] = Complex(
                [nATP * [energy], complex_dict['MP:sub']],
                compartment=membrane_pump.compartment,
            )
            # Complex3
            complex_dict['MP:prod:ATP'] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex4
            complex_dict['MP:ATP'] = Complex(
                [nATP * [energy], membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex5
            complex_dict['MP:ADP'] = Complex(
                [nATP * [waste], membrane_pump],
                compartment=membrane_pump.compartment,
            )

        # Active membrane protein transport
        # Sub + MP<--> Sub:MP
        prop_subMP = GeneralPropensity(
            f"kb_subMP * {substrate} * {membrane_pump} * "
            f"Heaviside({membrane_pump})",
            propensity_species=[substrate, membrane_pump],
            propensity_parameters=[kb_subMP],
        )
        rxn1_SubBinding = Reaction(
            [substrate, membrane_pump],
            [complex_dict['MP:sub']],
            propensity_type=prop_subMP,
        )

        rxn1_SubUnbinding = Reaction.from_massaction(
            inputs=[complex_dict['MP:sub']],
            outputs=[substrate, membrane_pump],
            k_forward=ku_subMP,
        )

        # Sub:MP + E <--> Sub:MP:E
        prop_subMPnATP = GeneralPropensity(
            f"kb_subMPnATP*{complex1}*{energy}*Heaviside({complex1})",
            propensity_species=[complex1, energy],
            propensity_parameters=[kb_subMPnATP],
        )
        rxn2_ATPbinding = Reaction(
            [complex1, nATP * [energy]],
            [complex_dict['MP:sub:ATP']],
            propensity_type=prop_subMPnATP,
        )

        rxn2_ATPunbinding = Reaction.from_massaction(
            inputs=[complex_dict['MP:sub:ATP']],
            outputs=[complex_dict['MP:sub'], nATP * [energy]],
            k_forward=ku_subMPnATP,
        )

        # Sub:MP:E <--> Prod:MP:E
        rxn3_transport = Reaction.from_massaction(
            inputs=[complex_dict['MP:sub:ATP']],
            outputs=[complex_dict['MP:prod:ATP']],
            k_forward=kf_trnspMP,
            k_reverse=kr_trnspMP
        )

        # Prod:MP:E <--> Prod + MP:E
        rxn4_SubRelease = Reaction.from_massaction(
            inputs=[complex_dict['MP:prod:ATP']],
            outputs=[complex_dict['MP:ATP'], product],
            k_forward=ku_prodMP,
            k_reverse=kb_prodMP,
        )

        # MP:E --> MP:W
        rxn5_atpHydro = Reaction.from_massaction(
            inputs=[complex_dict['MP:ATP']],
            outputs=[complex_dict['MP:ADP']],
            k_forward=ku_MP,
        )

        # MP:W <--> MP + W
        rxn6_reset = Reaction.from_massaction(
            inputs=[complex_dict['MP:ADP']],
            outputs=[nATP * [waste], membrane_pump],
            k_forward=ku_MP,
            k_reverse=kb_MP,
        )

        return [
            rxn1_SubBinding,
            rxn1_SubUnbinding,
            rxn2_ATPbinding,
            rxn2_ATPunbinding,
            rxn3_transport,
            rxn4_SubRelease,
            rxn5_atpHydro,
            rxn6_reset,
        ]
