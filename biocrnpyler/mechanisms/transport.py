# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.propensities import GeneralPropensity, ProportionalHillNegative
from ..core.reaction import Reaction
from ..core.species import Complex


class SimpleDiffusion(Mechanism):
    """Passive diffusion mechanism for substrate transport across membranes.

    A 'diffusion' mechanism that models simple passive diffusion of
    substrates through a membrane without requiring membrane proteins or
    energy. The transport is bidirectional and follows Fick's law of
    diffusion with equal forward and reverse rate constants.

    The reaction follows the schema:

    substrate <--> product

    where substrate and product represent the same species on opposite sides
    of the membrane.

    Parameters
    ----------
    name : str, default='simple_diffusion'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='diffusion'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('diffusion').

    See Also
    --------
    SimpleTransport : Passive transport through membrane channels.
    FacilitatedTransport_MM : Facilitated diffusion with carriers.
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
    >>> mechanism = bcp.SimpleDiffusion()
    >>> mixture = bcp.Mixture(
    ...     components=[O2],
    ...     mechanisms={'diffusion': mechanism},
    ...     parameters={'k_diff': 0.1}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name='simple_diffusion', mechanism_type='diffusion', **kwargs
    ):
        Mechanism.__init__(self, name, mechanism_type)

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
        """Generate reaction for simple diffusion.

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

        substrate <--> product (rates: 'k_diff', 'k_diff')

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


class MembraneProteinIntegration(Mechanism):
    """Membrane protein integration mechanism for protein insertion.

    A 'membrane_insertion' mechanism that models the integration of newly
    synthesized proteins into cellular membranes. Supports both monomeric
    and oligomeric membrane proteins, handling oligomerization before
    membrane insertion when required.

    The reaction schema depends on protein oligomeric state:

    For monomers (size = 1):
        monomer --> integral membrane protein

    For oligomers (size > 1):
        monomer * size <--> oligomer --> integral membrane protein

    Parameters
    ----------
    name : str, default='membrane_protein_integration'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='membrane_insertion'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('membrane_insertion').

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
    ...     direction='passive'
    ... )
    >>> mechanism = bcp.MembraneProteinIntegration()
    >>> mixture = bcp.Mixture(
    ...     components=[channel],
    ...     mechanisms={'membrane_insertion': mechanism},
    ...     parameters={
    ...         'kb_oligomer': 1.0, 'ku_oligomer': 0.1,
    ...         'kex': 0.5, 'kcat': 10.0
    ...     }
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='membrane_protein_integration',
        mechanism_type='membrane_insertion',
        **kwargs,
    ):
        Mechanism.__init__(self, name, mechanism_type)

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
            [oligomerization, integration].
            For monomers (size = 1): List of one reaction [integration].

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme depends on oligomeric state:

        For oligomers (size > 1):

        1. size * monomer <--> oligomer (rates: 'kb_oligomer',
           'ku_oligomer')
        2. oligomer --> product (ProportionalHillNegative with 'kex',
           'kcat')

        For monomers (size = 1):

        1. monomer --> product (ProportionalHillNegative with 'kex', 'kcat')

        The integration reaction uses `ProportionalHillNegative` kinetics with
        Hill coefficient n=4 to model saturation and product inhibition.

        """
        # Get Parameters
        kb_oligomer = component.get_parameter(
            'kb_oligomer', part_id=part_id, mechanism=self
        )
        ku_oligomer = component.get_parameter(
            'ku_oligomer', part_id=part_id, mechanism=self
        )
        kex = component.get_parameter('kex', part_id=part_id, mechanism=self)
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


class SimpleTransport(Mechanism):
    """Passive transport mechanism through membrane channel proteins.

    A 'transport' mechanism that models passive, bidirectional transport of
    substrates through membrane channel proteins. Unlike simple diffusion,
    this mechanism requires a membrane channel protein but does not consume
    energy. The channel acts catalytically, binding substrate and product
    but not being consumed.

    The reaction follows the schema:

    membrane_channel + substrate <--> membrane_channel + product

    Parameters
    ----------
    name : str, default='simple_membrane_protein_transport'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    SimpleDiffusion : Passive diffusion without proteins.
    FacilitatedTransport_MM : Transport with MM kinetics.
    PrimaryActiveTransport_MM : Energy-dependent active transport.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models passive transport through channel proteins such as
    ion channels, aquaporins, and other pore-forming proteins. The channel
    facilitates movement down concentration gradients without conformational
    changes or energy expenditure.

    The mechanism requires the membrane channel to have the 'passive'
    attribute, distinguishing it from active transporters and carriers that
    require different mechanisms.

    Common examples include:

    - Ion channels (K+, Na+, Ca2+ channels)
    - Aquaporins for water transport
    - Gap junctions between cells
    - Porins in bacterial outer membranes

    The transport is bidirectional with equal forward and reverse rate
    constants, reflecting passive equilibration across the membrane.

    Required parameters for this mechanism:

    - 'k_trnsp' : Transport rate constant (same for both directions)

    Examples
    --------
    Model potassium transport through an ion channel:

    >>> protein = bcp.IntegralMembraneProtein(
    ...     membrane_protein='Knck1',
    ...     product='K_channel',
    ...     direction='passive',
    ...     compartment='cytoplasm',
    ...     membrane_compartment='membrane',
    ...     attributes=['passive']
    ... )
    >>> channel = bcp.MembraneChannel(
    ...     integral_membrane_protein=protein.membrane_protein,
    ...     substrate='K',
    ...     direction='passive',
    ...     internal_compartment='cytoplasm',
    ...     external_compartment='external'
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[protein, channel],
    ...     mechanisms={
    ...         'membrane_insertion': bcp.MembraneProteinIntegration(),
    ...         'transport': bcp.SimpleTransport(),
    ...     },
    ...     parameters={'k_trnsp': 1.0},
    ...     parameter_file='mechanisms/transport_parameters.tsv',
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='simple_membrane_protein_transport',
        mechanism_type='transport',
        **kwargs,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, membrane_channel, substrate, product, **kwargs):
        """Generate species for simple transport.

        Returns the membrane channel, substrate, and product species
        involved in the transport reaction. Validates that the channel has
        the 'passive' attribute.

        Parameters
        ----------
        membrane_channel : Species
            The membrane channel protein through which transport occurs.
            Must have 'passive' as its first attribute.
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

        Raises
        ------
        ValueError
            If membrane_channel does not have 'passive' as its first
            attribute, indicating it should use FacilitatedTransport_MM
            instead.

        """
        if membrane_channel.attributes[0] != 'passive':
            raise ValueError(
                "Protein is not classified as a channel with passive "
                "transport of small molecules. Use mechanism "
                "FacilitatedTransport_MM instead"
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
        **kwargs,
    ):
        """Generate reaction for simple membrane protein transport.

        Creates a single reversible mass-action reaction representing
        passive transport through a membrane channel with equal forward and
        reverse rate constants. The channel acts catalytically and is not
        consumed.

        Parameters
        ----------
        membrane_channel : Species
            The membrane channel protein facilitating transport.
        substrate : Species
            The substrate species being transported.
        product : Species
            The product species after transport.
        component : Component, optional
            Component containing parameter values. Required if k_trnsp is
            not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is
            provided, defaults to component.name.
        k_trnsp : Parameter or float, optional
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
            If component is None and k_trnsp is not provided.

        Notes
        -----
        The reaction has equal forward and reverse rate constants:

        membrane_channel + substrate <--> membrane_channel + product
        (rates: 'k_trnsp', 'k_trnsp')

        The membrane channel appears on both sides of the reaction,
        indicating it acts catalytically and is recycled.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        if component is None and (k_trnsp is None):
            raise ValueError(
                "Must pass in a Component or values for k_trnsp."
            )
        if k_trnsp is None:
            k_trnsp = component.get_parameter(
                'k_trnsp', part_id=part_id, mechanism=self
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


class FacilitatedTransport_MM(Mechanism):
    """Facilitated diffusion mechanism with Michaelis-Menten kinetics.

    A 'transport' mechanism that models facilitated diffusion of substrates
    through membrane carrier proteins. Unlike simple channels, carriers
    undergo conformational changes to transport substrates across membranes.
    The mechanism follows Michaelis-Menten kinetics with explicit substrate
    and product binding steps.

    The reaction follows the schema:

    Sub + MC <--> Sub:MC --> Prod:MC --> Prod + MC

    where MC is the membrane carrier protein.

    Parameters
    ----------
    name : str, default='facilitated_membrane_protein_transport'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    SimpleTransport : Passive transport through channels.
    PrimaryActiveTransport_MM : Energy-dependent active transport.
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

    The mechanism uses a GeneralPropensity with a Heaviside function for
    the initial binding step to enforce directionality based on
    concentration gradients.

    Required parameters for this mechanism:

    - 'kb_subMC' : Forward binding rate for substrate to membrane carrier
    - 'ku_subMC' : Unbinding rate for substrate from carrier
    - 'k_trnspMC' : Conformational change rate (transport step)
    - 'ku_prodMC' : Unbinding rate for product from carrier

    Examples
    --------
    Model glucose transport through a GLUT transporter:

    >>> glc_in = bcp.Species('glucose', compartment='cytoplasm')
    >>> glc_out = bcp.Species('glucose', compartment='external')
    >>> carrier = bcp.MembraneChannel(
    ...     integral_membrane_protein='GlucoseTransporter',
    ...     substrate=glc_out,
    ...     external_compartment='external',
    ...     internal_compartment='cytoplasm',
    ...     direction='importer'
    ... )
    >>> mechanism = bcp.FacilitatedTransport_MM()
    >>> mixture = bcp.Mixture(
    ...     components=[carrier],
    ...     mechanisms={'transport': mechanism},
    ...     parameters={
    ...         'kb_subMC': 1.0, 'ku_subMC': 0.5,
    ...         'k_trnspMC': 0.8, 'ku_prodMC': 0.5
    ...     }
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='facilitated_membrane_protein_transport',
        mechanism_type='transport',
        **kwargs,
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self,
        membrane_carrier,
        substrate,
        product,
        complex_dict=None,
        **kwargs,
    ):
        """Generate species for facilitated transport.

        Creates species for the membrane carrier, substrate, product, and
        the two intermediate complexes formed during the transport cycle.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein that facilitates transport.
        substrate : Species
            The substrate species being transported (typically intracellular
            side).
        product : Species
            The product species after transport (typically extracellular
            side). Usually the same molecular species as substrate but in a
            different compartment.
        complex_dict : dict, optional
            Pre-defined dictionary of complex species with keys 'sub:MC' and
            'prod:MC'. If None, complexes are automatically created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [membrane_carrier, substrate, product,
            complex_array] where complex_array is a list of two Complex
            species: [substrate:carrier, product:carrier].

        Notes
        -----
        The method creates two complex species representing intermediates in
        the transport cycle:

        1. sub:MC : substrate:membrane_carrier complex
        2. prod:MC : product:membrane_carrier complex

        """
        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['sub:MC'] = Complex(
                [substrate, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )
            # Complex2
            complex_dict['prod:MC'] = Complex(
                [product, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )

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
        **kwargs,
    ):
        """Generate reactions for facilitated transport.

        Creates four reactions representing the complete transport cycle:
        substrate binding, substrate unbinding, conformational change
        (transport), and product release.

        Parameters
        ----------
        membrane_carrier : Species
            The membrane carrier protein facilitating transport.
        substrate : Species
            The substrate species being transported.
        product : Species
            The product species after transport.
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
            List of four reactions: [substrate_binding, substrate_unbinding,
            transport_step, product_release].

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme follows this pathway:

        1. MC + Sub <--> MC:Sub (GeneralPropensity with Heaviside function
           using 'kb_subMC')
        2. MC:Sub --> MC + Sub (irreversible, rate: 'ku_subMC')
        3. MC:Sub --> MC:Prod (conformational change, rate: 'k_trnspMC')
        4. MC:Prod --> MC + Prod (irreversible, rate: 'ku_prodMC')

        The initial binding step uses a GeneralPropensity with a Heaviside
        function to enforce concentration gradient-driven directionality.
        The Heaviside function ensures transport only occurs when substrate
        concentration exceeds product concentration.

        """
        # Get Parameters
        kb_subMC = component.get_parameter(
            'kb_subMC', part_id=part_id, mechanism=self
        )
        ku_subMC = component.get_parameter(
            'ku_subMC', part_id=part_id, mechanism=self
        )
        k_trnspMC = component.get_parameter(
            'k_trnspMC', part_id=part_id, mechanism=self
        )
        ku_prodMC = component.get_parameter(
            'ku_prodMC', part_id=part_id, mechanism=self
        )

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['sub:MC'] = Complex(
                [substrate, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )
            # Complex2
            complex_dict['prod:MC'] = Complex(
                [product, membrane_carrier],
                compartment=membrane_carrier.compartment,
            )

        # Facilitated membrane protein transport
        # Sub + MC --> Sub:MC
        prop_subMC = GeneralPropensity(
            f"kb_subMC * {substrate} * {membrane_carrier} * "
            f"Heaviside({substrate}-{product}) "
            f"- kb_subMC * {product} * {membrane_carrier} * "
            f"Heaviside({substrate}-{product})",
            propensity_species=[product, substrate, membrane_carrier],
            propensity_parameters=[kb_subMC],
        )
        binding_rxn1 = Reaction(
            [substrate, membrane_carrier],
            [complex_dict['sub:MC']],
            propensity_type=prop_subMC,
        )

        # Sub:MC --> Sub + MC
        unbinding_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['sub:MC']],
            outputs=[membrane_carrier, substrate],
            k_forward=ku_subMC,
        )

        # Sub:MC --> Prod:MC
        transport_rxn = Reaction.from_massaction(
            inputs=[complex_dict['sub:MC']],
            outputs=[complex_dict['prod:MC']],
            k_forward=k_trnspMC,
        )

        # MC:Prod --> MC + Prod
        unbinding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict['prod:MC']],
            outputs=[product, membrane_carrier],
            k_forward=ku_prodMC,
        )

        return [binding_rxn1, unbinding_rxn1, transport_rxn, unbinding_rxn2]


class PrimaryActiveTransport_MM(Mechanism):
    """Primary active transport mechanism with ATP-dependent pumping.

    A 'transport' mechanism that models primary active transport where
    substrates are moved against their concentration gradients using energy
    from ATP hydrolysis. The mechanism follows Michaelis-Menten kinetics
    with explicit binding, ATP hydrolysis, conformational change, and
    product release steps.

    The reaction follows the schema:

    Sub + MP <--> Sub:MP + E --> Sub:MP:E --> MP:Prod:E
    --> Prod + MP:W --> Prod + MP + W

    where MP is the membrane pump, E is ATP (energy), and W is ADP (waste).

    Parameters
    ----------
    name : str, default='active_membrane_protein_transport'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transport'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transport').

    See Also
    --------
    FacilitatedTransport_MM : Passive facilitated diffusion.
    SimpleTransport : Passive channel transport.
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

    The binding steps use GeneralPropensity with Heaviside functions to
    ensure proper directionality based on species concentrations.

    Required parameters for this mechanism:

    - 'kb_subMP' : Forward binding rate for substrate to membrane pump
    - 'ku_subMP' : Unbinding rate for substrate from pump
    - 'kb_subMPnATP' : Forward binding rate for ATP to substrate:pump
      complex
    - 'ku_subMPnATP' : Unbinding rate for ATP from substrate:pump complex
    - 'k_trnspMP' : Conformational change rate (transport step)
    - 'ku_prodMP' : Unbinding rate for product from pump
    - 'ku_MP' : Unbinding rate for ADP from pump

    Examples
    --------
    Model active sodium transport by Na+/K+-ATPase:

    >>> pump = bcp.MembranePump(
    ...     membrane_pump='NaK_ATPase',
    ...     substrate='Na',
    ...     direction='exporter',
    ...     ATP=1
    ... )
    >>> mechanism = bcp.PrimaryActiveTransport_MM()
    >>> mixture = bcp.Mixture(
    ...     components=[pump],
    ...     mechanisms={'transport': mechanism},
    ...     parameters={
    ...         'kb_subMP': 1.0, 'ku_subMP': 0.1,
    ...         'kb_subMPnATP': 1.0, 'ku_subMPnATP': 0.1,
    ...         'k_trnspMP': 0.5, 'ku_prodMP': 1.0,
    ...         'ku_MP': 1.0
    ...     }
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='active_membrane_protein_transport',
        mechanism_type='transport',
        **kwargs,
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
            Pre-defined dictionary of complex species with keys 'Pump:Sub',
            'Pump:Sub:ATP', 'Pump:Prod:ATP', and 'Pump:ADP'. If None,
            complexes are automatically created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list
            List containing [membrane_pump, substrate, product, energy,
            waste, complex_array] where complex_array is a list of four
            Complex species generated.

        Notes
        -----
        The method creates four complex species representing intermediates
        in the active transport cycle:

        1. Pump:Sub : membrane_pump:substrate complex
        2. Pump:Sub:ATP : membrane_pump:substrate:nATP complex
        3. Pump:Prod:ATP : membrane_pump:product:nATP complex
        4. Pump:ADP : membrane_pump:nADP complex

        The number of ATP/ADP molecules (nATP) is determined by the
        membrane_pump.ATP attribute.

        """
        nATP = membrane_pump.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}
            # Complex1
            complex_dict['Pump:Sub'] = Complex(
                [substrate, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex2
            complex_dict['Pump:Sub:ATP'] = Complex(
                [nATP * [energy], complex_dict['Pump:Sub']],
                compartment=membrane_pump.compartment,
            )
            # Complex3
            complex_dict['Pump:Prod:ATP'] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            # Complex4
            complex_dict['Pump:ADP'] = Complex(
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
            waste,
            complex_array,
        ]

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

        Creates seven reactions representing the complete ATP-driven
        transport cycle: substrate binding/unbinding, ATP
        binding/unbinding, conformational change (transport), product
        release, and ADP release.

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
            List of seven reactions: [substrate_binding,
            substrate_unbinding, ATP_binding, ATP_unbinding, transport_step,
            product_release, ADP_release].

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The reaction scheme follows this pathway:

        1. MP + Sub <--> MP:Sub (`GeneralPropensity with Heaviside using
           'kb_subMP', reverse rate: 'ku_subMP')
        2. MP:Sub + nATP <--> MP:Sub:nATP (GeneralPropensity with Heaviside
           using 'kb_subMPnATP', reverse rate: 'ku_subMPnATP')
        3. MP:Sub:nATP --> MP:Prod:nATP (conformational change, rate:
           'k_trnspMP')
        4. MP:Prod:nATP --> MP:nADP + Prod (product release, rate:
           'ku_prodMP')
        5. MP:nADP --> MP + nADP (ADP release, rate: 'ku_MP')

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
        k_trnspMP = component.get_parameter(
            'k_trnspMP', part_id=part_id, mechanism=self
        )
        ku_prodMP = component.get_parameter(
            'ku_prodMP', part_id=part_id, mechanism=self
        )
        ku_MP = component.get_parameter(
            'ku_MP', part_id=part_id, mechanism=self
        )

        nATP = membrane_pump.ATP

        if complex_dict is None:
            # Create empty dictionary for complexes
            complex_dict = {}

            # Complex1
            complex_dict['Pump:Sub'] = Complex(
                [substrate, membrane_pump],
                compartment=membrane_pump.compartment,
            )
            complex1 = complex_dict['Pump:Sub']

            # Complex2
            complex_dict['Pump:Sub:ATP'] = Complex(
                [nATP * [energy], complex_dict['Pump:Sub']],
                compartment=membrane_pump.compartment,
            )
            # Complex3
            complex_dict['Pump:Prod:ATP'] = Complex(
                [nATP * [energy], product, membrane_pump],
                compartment=membrane_pump.compartment,
            )

            # Complex4
            complex_dict['Pump:ADP'] = Complex(
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
        binding_rxn1 = Reaction(
            [substrate, membrane_pump],
            [complex_dict['Pump:Sub']],
            propensity_type=prop_subMP,
        )

        unbinding_rxn1 = Reaction.from_massaction(
            inputs=[complex_dict['Pump:Sub']],
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
            [complex_dict['Pump:Sub:ATP']],
            propensity_type=prop_subMPnATP,
        )

        unbinding_rxn2 = Reaction.from_massaction(
            inputs=[complex_dict['Pump:Sub:ATP']],
            outputs=[complex_dict['Pump:Sub'], nATP * [energy]],
            k_forward=ku_subMPnATP,
        )

        # Sub:MP:E --> Prod:MP:E
        transport_rxn = Reaction.from_massaction(
            inputs=[complex_dict['Pump:Sub:ATP']],
            outputs=[complex_dict['Pump:Prod:ATP']],
            k_forward=k_trnspMP,
        )
        # Prod:MP:E--> Prod+MP:W
        unbinding_rxn3 = Reaction.from_massaction(
            inputs=[complex_dict['Pump:Prod:ATP']],
            outputs=[complex_dict['Pump:ADP'], product],
            k_forward=ku_prodMP,
        )
        # MP:W --> MP+W
        unbinding_rxn4 = Reaction.from_massaction(
            inputs=[complex_dict['Pump:ADP']],
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


# Legacy class names
Simple_Diffusion = SimpleDiffusion
Membrane_Protein_Integration = MembraneProteinIntegration
Simple_Transport = SimpleTransport
Facilitated_Transport_MM = FacilitatedTransport_MM
Primary_Active_Transport_MM = PrimaryActiveTransport_MM
