# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from typing import List, Tuple, Union
from warnings import resetwarnings, warn

from ..mechanisms.global_mechanisms import GlobalMechanism
from .chemical_reaction_network import ChemicalReactionNetwork
from .compartment import Compartment
from .component import Component
from .mechanism import Mechanism
from .parameter import ParameterDatabase
from .reaction import Reaction
from .species import Species


class Mixture(object):
    """Container for components, mechanisms, and parameters in a CRN model.

    A Mixture holds together all components (DNA, RNA, Protein, etc.),
    mechanisms (transcription, translation, binding, etc.), and parameters
    needed to compile a chemical reaction network (CRN). Mixtures can include
    default mechanisms that apply to all components, as well as global
    mechanisms that affect all species (e.g., dilution, degradation).

    Parameters
    ----------
    name : str, default=''
        Name of the mixture for identification and parameter lookup.
    mechanisms : dict, list, or Mechanism, optional
        Default mechanisms for components in this mixture. Can be a dict with
        mechanism types (str) as keys and `Mechanism` objects as values, a
        list of `Mechanism` objects, or a single `Mechanism`.
    components : list of Component or Component, optional
        Components to include in the mixture. Components are deep-copied when
        added to prevent modification of original objects.
    parameters : dict, optional
        Dictionary of parameter values. Keys follow the format
        (mechanism, part_id, param_name).
    compartment : Compartment, optional
        Default compartment for all components and species in this mixture.
    parameter_file : str, optional
        Path to a CSV or TSV file containing parameters to load.
    overwrite_parameters : bool, default=False
        If True, parameters from file/dict overwrite existing parameters.
        If False, existing parameters are preserved.
    global_mechanisms : dict, list, or GlobalMechanism, optional
        Global mechanisms that apply to all species after component
        compilation (e.g., dilution, global degradation). Can be a dict,
        list, or single `GlobalMechanism`.
    species : list of Species or Species, optional
        Additional species to add directly to the CRN without going through
        component compilation.
    initial_condition_dictionary : dict, optional
        Dictionary mapping species to initial concentration values. Deprecated
        in favor of using parameters with mechanism='initial concentration'.
    global_component_enumerators : list, optional
        List of global component enumerators for advanced component generation
        patterns (e.g., creating all pairwise interactions).
    global_recursion_depth : int, default=4
        Maximum recursion depth for global component enumeration during
        compilation.
    local_recursion_depth : int, optional
        Maximum recursion depth for local component enumeration. If None,
        defaults to `global_recursion_depth + 2`.

    Attributes
    ----------
    name : str
        Name of the mixture.
    compartment : Compartment or None
        Default compartment for the mixture.
    components : list of Component
        List of components in the mixture (deep copies of added components).
    mechanisms : dict
        Dictionary of default mechanisms, keyed by mechanism type (str).
    global_mechanisms : dict
        Dictionary of global mechanisms, keyed by mechanism type (str).
    parameter_database : ParameterDatabase
        Database storing all parameters for this mixture.
    added_species : list of Species
        List of species added directly to the mixture.
    global_component_enumerators : list
        List of global component enumerators.
    global_recursion_depth : int
        Recursion depth for global component enumeration.
    local_recursion_depth : int
        Recursion depth for local component enumeration.
    crn : ChemicalReactionNetwork or None
        The compiled CRN, created by calling `compile_crn`.

    See Also
    --------
    Component : Base class for biomolecular components.
    Mechanism : Base class for reaction generation schemas.
    GlobalMechanism : Mechanisms that apply to all species.
    ChemicalReactionNetwork : Container for species and reactions.

    Notes
    -----
    Components added to a Mixture are deep-copied to prevent unintended
    modifications. Each component's `mixture` attribute is set to reference
    this Mixture, allowing components to access default mechanisms and
    parameters.

    The compilation process follows these steps:

    1. Global component enumeration (e.g., generating all interactions)
    2. Local component enumeration (e.g., DNA --> RNA --> Protein)
    3. Species generation from all enumerated components
    4. Reaction generation from all enumerated components
    5. Application of global mechanisms to all species

    Examples
    --------
    Create a basic mixture with components and mechanisms:

    >>> # Create a simple transcription-translation mixture
    >>> mixture = bcp.Mixture(
    ...     name="txtl_extract",
    ...     mechanisms={
    ...         'transcription': bcp.SimpleTranscription(),
    ...         'translation': bcp.SimpleTranslation()
    ...     },
    ...     parameters={'ktx': 0.05, 'ktl': 0.01}
    ... )

    Add components to the mixture:

    >>> promoter = bcp.Promoter("pconst")
    >>> gene = bcp.DNAconstruct([promoter, bcp.RBS('rbs'), bcp.CDS('gfp')])
    >>> mixture.add_component(gene)

    Compile the CRN:

    >>> crn = mixture.compile_crn()
    >>> print(
    ...     f"CRN has {len(crn.species)} species "
    ...     f"and {len(crn.reactions)} reactions")

    """

    def __init__(
        self,
        name='',
        mechanisms=None,
        components=None,
        parameters=None,
        compartment: Compartment = None,
        parameter_file=None,
        overwrite_parameters=False,
        global_mechanisms=None,
        species=None,
        initial_condition_dictionary=None,
        global_component_enumerators=None,
        global_recursion_depth=4,
        local_recursion_depth=None,
    ):
        # Initialize instance variables
        self.name = name  # Save the name of the mixture
        self.compartment = compartment

        # recursion depth for global component enumeration
        self.global_recursion_depth = global_recursion_depth

        if local_recursion_depth is None:
            self.local_recursion_depth = self.global_recursion_depth + 2

        # process the components
        if components is None and not hasattr(self, '_components'):
            self.components = []
        else:
            self.add_components(components)

        # process mechanisms:
        if mechanisms is None and not hasattr(self, '_mechanisms'):
            self.mechanisms = {}
        else:
            self.add_mechanisms(mechanisms, overwrite=True)

        # process global_mechanisms:
        # Global mechanisms are applied just once ALL species generated from
        # components inside a mixture
        # Global mechanisms should be used rarely, and with care. An example
        # usecase is degradation via dilution.
        if global_mechanisms is None and not hasattr(
            self, '_global_mechanisms'
        ):
            self.global_mechanisms = {}
        elif global_mechanisms is not None:
            self.add_mechanisms(global_mechanisms)

        # global component enumerators
        if global_component_enumerators is None:
            self.global_component_enumerators = []
        else:
            self.global_component_enumerators = global_component_enumerators

        # process the species
        self.add_species(species)

        # Create a paraemter database
        self.parameter_database = ParameterDatabase(
            parameter_file=parameter_file,
            parameter_dictionary=parameters,
            overwrite_parameters=overwrite_parameters,
        )

        # CRN is stored here during compilation
        self.crn = None

    def add_species(self, species: Union[List[Species], Species]):
        """Add species directly to the mixture without component compilation.

        Parameters
        ----------
        species : Species or list of Species
            Species object(s) to add directly to the mixture. These species
            will be included in the CRN during compilation.

        Raises
        ------
        AssertionError
            If any element in the list is not a Species object.

        Notes
        -----
        Species added this way bypass component enumeration and are added
        directly to the CRN during `compile_crn`.

        """
        if not hasattr(self, 'added_species'):
            self.added_species = []

        if species is not None:
            if not isinstance(species, list):
                species_list = [species]
            else:
                species_list = species

            assert all(
                isinstance(x, Species) for x in species_list
            ), 'only Species type is accepted!'

            self.added_species += species_list

    def set_species(
        self,
        species: Union[Species, str],
        material_type=None,
        attributes=None,
    ):
        """Convert various inputs into Species objects.

        Parameters
        ----------
        species : Species, str, or Component
            The species to convert. Can be a `Species` object (returned
            as-is), a string (creates new Species), or a `Component` (extracts
            its species).
        material_type : str, optional
            Material type for the species (e.g., 'dna', 'rna', 'protein').
            Only used when creating new Species from strings.
        attributes : list of str, optional
            Attributes to assign to the species. Only used when creating
            new Species from strings.

        Returns
        -------
        Species
            The converted Species object.

        Raises
        ------
        ValueError
            If the input cannot be converted to a valid Species.

        """
        if isinstance(species, Species):
            return species
        elif isinstance(species, str):
            return Species(
                name=species,
                material_type=material_type,
                attributes=attributes,
            )
        elif (
            isinstance(species, Component)
            and species.get_species() is not None
        ):
            return species.get_species()
        else:
            raise ValueError(
                "Invalid Species: string, chemical_reaction_network.Species "
                "or Component with implemented .get_species() required as "
                "input."
            )

    @property
    def components(self):
        return self._components

    @components.setter
    def components(self, components):
        self._components = []
        self.add_components(components)

    def add_component(self, component):
        """Add a single component to the mixture.

        Parameters
        ----------
        component : Component or list of Component
            Component object to add to the mixture. If a list is provided,
            calls `add_components` instead. The component is deep-copied
            before being added.

        Raises
        ------
        AssertionError
            If the component is not a Component object.
        ValueError
            If a component with the same type and name already exists in
            the mixture.

        Notes
        -----
        Components are deep-copied when added to prevent modification of the
        original component. The copied component's `mixture` attribute is set
        to this Mixture, and its `compartment` is set to the mixture's
        compartment.

        """
        if not hasattr(self, '_components'):
            self.components = []

        if isinstance(component, list):
            self.add_components(component)
        else:
            assert isinstance(component, Component), (
                'the object: %s passed into mixture as component must '
                + 'be of the class Component' % str(component)
            )

            # Check if component is already in self._components
            for comp in self._components:
                if (
                    type(comp) is type(component)
                    and comp.name == component.name
                ):
                    raise ValueError(
                        f"{comp} of same type and name already in Mixture!"
                    )
            else:
                # Components are copied before being added to Mixtures
                component_copy = copy.deepcopy(component)
                component_copy.set_mixture(self)
                component_copy.compartment = self.compartment
                self.components.append(component_copy)

    def get_mechanism(self, mechanism_type):
        """Retrieve a mechanism by type from the mixture.

        Parameters
        ----------
        mechanism_type : str
            The type identifier of the mechanism to retrieve (e.g.,
            'transcription', 'translation', 'binding').

        Returns
        -------
        Mechanism or None
            The requested mechanism object, or None if not found.

        Raises
        ------
        TypeError
            If `mechanism_type` is not a string.

        """
        if not isinstance(mechanism_type, str):
            raise TypeError(
                f"mechanism_type must be a string. Received {mechanism_type}."
            )

        if mechanism_type in self.mechanisms:
            return self.mechanisms[mechanism_type]
        else:
            return None

    def add_components(self, components: Union[List[Component], Component]):
        """Add multiple components to the mixture.

        Parameters
        ----------
        components : Component or list of Component
            Component object(s) to add to the mixture. Each component is
            deep-copied before being added.

        Raises
        ------
        ValueError
            If `components` is not a Component, list of Components, or if
            any duplicate components are detected.

        See Also
        --------
        add_component : Add a single component to the mixture.

        """
        if isinstance(components, Component):
            self.add_component(components)
        elif isinstance(components, List):
            for component in components:
                self.add_component(component)
        else:
            raise ValueError(
                f"add_components expected a list of Components. "
                f"Received {components}"
            )

    def get_component(self, component=None, name=None, index=None):
        """Retrieve components from the mixture by various criteria.

        Exactly one of the three parameters must be provided.

        Parameters
        ----------
        component : Component, optional
            A component instance to search for. Returns components with
            matching type and name.
        name : str, optional
            Name of the component to search for. Returns all components
            with this name.
        index : int, optional
            Index of the component in the mixture's component list.

        Returns
        -------
        Component, list of Component, or None
            - Single Component if exactly one match is found or index is used
            - List of Components if multiple matches are found
            - None if no matches are found

        Raises
        ------
        ValueError
            If zero or more than one parameter is provided, or if parameters
            are of incorrect types.

        """
        if [component, name, index].count(None) != 2:
            raise ValueError(
                f"get_component requires a single keyword. "
                f"Received component={component}, name={name}, index={index}."
            )
        if not (isinstance(component, Component) or component is None):
            raise ValueError(
                f"component must be of type Component. Received {component}."
            )
        if not (isinstance(name, str) or name is None):
            raise ValueError(f"name must be of type str. Received {name}.")
        if not (isinstance(index, int) or index is None):
            raise ValueError(f"index must be of type int. Received {index}.")

        matches = []
        if index is not None:
            matches.append(self.components[index])
        else:
            for comp in self.components:
                if component is not None:
                    if (
                        type(comp) is type(component)
                        and comp.name == component.name
                    ):
                        matches.append(comp)
                elif name is not None:
                    if comp.name == name:
                        matches.append(comp)
        if len(matches) == 0:
            return None
        elif len(matches) == 1:
            return matches[0]
        else:
            warn(
                "get_component found multiple matching components. "
                "A list has been returned."
            )
            return matches

    @property
    def mechanisms(self):
        """Mechanism: Stores mixture mechanisms."""
        return self._mechanisms

    @mechanisms.setter
    def mechanisms(self, mechanisms):
        self._mechanisms = {}
        self.add_mechanisms(mechanisms, overwrite=True)

    def add_mechanism(self, mechanism, mech_type=None, overwrite=False):
        """Add a mechanism to the mixture's mechanism dictionary.

        Parameters
        ----------
        mechanism : Mechanism or GlobalMechanism
            The mechanism object to add. If a `GlobalMechanism` is provided,
            it is automatically added to the global mechanisms dictionary.
        mech_type : str, optional
            The type key under which to store the mechanism. If None, uses
            the mechanism's `mechanism_type` attribute.
        overwrite : bool, default=False
            If True, replaces any existing mechanisms with the same keys. If
            False, raises ValueError when keys already exist. If None, ignores
            mechanisms that already exist.

        Raises
        ------
        TypeError
            If `mechanism` is not a Mechanism object, or if `mech_type` is
            not a string.
        ValueError
            If mechanism key already exists and `overwrite` is None.

        See Also
        --------
        add_global_mechanism : Add a global mechanism specifically.

        """
        if not hasattr(self, '_mechanisms'):
            self._mechanisms = {}

        if not isinstance(mechanism, Mechanism):
            raise TypeError(
                f"mechanism must be a Mechanism. Received {mechanism}."
            )

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(
                f"mechanism keys must be strings. Received {mech_type}"
            )

        if isinstance(mechanism, GlobalMechanism):
            self.add_global_mechanism(mechanism, mech_type, overwrite)
        elif isinstance(mechanism, Mechanism):
            if mech_type in self._mechanisms and not overwrite:
                if overwrite is False:
                    raise ValueError(
                        f"mech_type {mech_type} already in Mixture {self}. "
                        f"To overwrite, use keyword overwrite=True."
                    )
            else:
                self._mechanisms[mech_type] = copy.deepcopy(mechanism)

    def add_mechanisms(self, mechanisms, overwrite=False):
        """Add multiple mechanisms to the mixture.

        Accepts mechanisms as a single object, list, or dictionary and adds
        them to the mixture's mechanism dictionary. Can handle both regular
        `Mechanism` and `GlobalMechanism` objects.

        Parameters
        ----------
        mechanisms : Mechanism, GlobalMechanism, dict, or list
            The mechanism(s) to add. Can be a single mechanism, a dict with
            mechanism types as keys and mechanisms as values, or a list of
            mechanisms.
        overwrite : bool, default=False
            If True, replaces any existing mechanisms with the same keys. If
            False, raises ValueError when keys already exist. If None, ignores
            mechanisms that already exist.

        Raises
        ------
        ValueError
            If `mechanisms` is not a valid type, or if mechanism key conflicts
            occur with `overwrite=False`.

        See Also
        --------
        add_mechanism : Add a single mechanism to the mixture.

        """
        if isinstance(mechanisms, Mechanism):
            self.add_mechanism(mechanisms, overwrite=overwrite)
        elif isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_mechanism(
                    mechanisms[mech_type], mech_type, overwrite=overwrite
                )
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_mechanism(mech, overwrite=overwrite)
        else:
            raise ValueError(
                f"add_mechanisms expected a list of Mechanisms. "
                f"Recieved {mechanisms}"
            )

    @property
    def global_mechanisms(self):
        """Mechanism: Stores global mechanisms in the mixture."""
        return self._global_mechanisms

    @global_mechanisms.setter
    def global_mechanisms(self, mechanisms):
        self._global_mechanisms = {}
        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_global_mechanism(
                    mechanisms[mech_type], mech_type, overwrite=True
                )
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_global_mechanism(mech, overwrite=True)

    def add_global_mechanism(
        self, mechanism, mech_type=None, overwrite=False
    ):
        """Add a global mechanism to the mixture.

        Global mechanisms are applied to all species after component
        compilation, enabling operations like dilution or global degradation.

        Parameters
        ----------
        mechanism : GlobalMechanism
            The global mechanism object to add. Must be a `GlobalMechanism`
            instance.
        mech_type : str, optional
            The type key under which to store the mechanism. If None, uses
            the mechanism's `mechanism_type` attribute.
        overwrite : bool, default=False
            If True, replaces any existing global mechanism with the same key.
            If False, raises ValueError when key already exists.

        Raises
        ------
        TypeError
            If `mechanism` is not a GlobalMechanism object, or if `mech_type`
            is not a string.
        ValueError
            If mechanism key already exists and `overwrite` is False.

        Notes
        -----
        Global mechanisms are applied during `compile_crn` after all
        component reactions have been generated.

        """
        if not hasattr(self, '_global_mechanisms'):
            self._global_mechanisms = {}

        if not isinstance(mechanism, GlobalMechanism):
            raise TypeError(
                f"mechanism must be a GlobalMechanism. Recieved {mechanism}."
            )

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(
                f"mechanism keys must be strings. Recieved {mech_type}"
            )

        if mech_type in self._mechanisms and not overwrite:
            raise ValueError(
                f"mech_type {mech_type} already in Mixture {self}. "
                f"To overwrite, use keyword overwrite = True."
            )
        else:
            self._global_mechanisms[mech_type] = copy.deepcopy(mechanism)

    def update_parameters(
        self, parameter_file=None, parameters=None, overwrite_parameters=True
    ):
        """Update the parameter database with new parameters.

        Parameters
        ----------
        parameter_file : str, optional
            Path to a CSV or TSV file containing parameters to load.
        parameters : dict, optional
            Dictionary of parameters to add. Keys follow the format
            (mechanism, part_id, param_name).
        overwrite_parameters : bool, default=True
            If True, new parameter values overwrite existing ones. If False,
            existing parameters are preserved.

        """
        if parameter_file is not None:
            self.parameter_database.load_parameters_from_file(
                parameter_file, overwrite_parameters=overwrite_parameters
            )

        if parameters is not None:
            self.parameter_database.load_parameters_from_dictionary(
                parameters, overwrite_parameters=overwrite_parameters
            )

    def get_parameter(self, mechanism, part_id, param_name):
        """Retrieve a parameter from the mixture's parameter database.

        Parameters
        ----------
        mechanism : str
            Mechanism identifier for the parameter lookup key.
        part_id : str
            Part identifier for the parameter lookup key.
        param_name : str
            Name of the parameter to retrieve.

        Returns
        -------
        Parameter or None
            The parameter object, or None if not found.

        """
        param = self.parameter_database.find_parameter(
            mechanism, part_id, param_name
        )
        return param

    def get_initial_concentration(
        self, S: Union[List, Species], component=None
    ):
        """Determine initial concentrations using parameter hierarchy.

        Searches for initial concentration parameters for species following a
        hierarchical lookup strategy, defaulting to 0 if not found.

        Parameters
        ----------
        S : Species or list of Species
            Species object(s) for which to find initial concentrations. Lists
            are automatically flattened.
        component : Component, optional
            The component that generated the species. Used for
            component-specific parameter lookup.

        Returns
        -------
        dict
            Dictionary mapping each Species object to its initial
            concentration value (float).

        Raises
        ------
        ValueError
            If any element in `S` is not a Species object.

        Notes
        -----
        The parameter lookup hierarchy is:

        1. Component's `ParameterDatabase` with `mechanism='initial
           concentration'`, `part_id=mixture.name`, and parameter names:
           `str(s)`, `s.name`, or `component.name` (where `s` is the
           component's primary species)

        2. Mixture's `ParameterDatabase` with the same keys

        3. Defaults to 0 if not found

        """
        if isinstance(S, Species):
            S = [S]

        # flatten the species list
        S = Species.flatten_list(S)

        init_conc_dict = {}
        for s in S:
            if not isinstance(s, Species):
                raise ValueError(
                    f"{s} is not a Species! Can only find initial "
                    f"concentration of a Species."
                )
            init_conc = None
            # 1 Check the component
            if component is not None:
                init_conc = component.get_parameter(
                    param_name=str(s),
                    part_id=self.name,
                    mechanism='initial concentration',
                    check_mixture=False,
                    return_none=True,
                )
                if init_conc is None:
                    init_conc = component.get_parameter(
                        param_name=s.name,
                        part_id=self.name,
                        mechanism='initial concentration',
                        check_mixture=False,
                        return_none=True,
                    )
                if init_conc is None and component.get_species() == s:
                    init_conc = component.get_parameter(
                        param_name=component.name,
                        part_id=self.name,
                        mechanism='initial concentration',
                        check_mixture=False,
                        return_none=True,
                    )

            # 2 Check self
            if init_conc is None:
                init_conc = self.get_parameter(
                    param_name=str(s),
                    part_id=self.name,
                    mechanism='initial concentration',
                )
                if init_conc is None:
                    init_conc = self.get_parameter(
                        param_name=s.name,
                        part_id=self.name,
                        mechanism='initial concentration',
                    )
                if (
                    init_conc is None
                    and component is not None
                    and component.get_species() == s
                ):
                    init_conc = self.get_parameter(
                        param_name=component.name,
                        part_id=self.name,
                        mechanism='initial concentration',
                    )

            if init_conc is None:
                init_conc = 0

            init_conc_dict[s] = init_conc

        return init_conc_dict

    def add_species_to_crn(
        self,
        new_species,
        component=None,
        no_initial_concentrations=False,
        copy_species=True,
        compartment=None,
    ):
        """Add species to the CRN with initial concentrations.

        Helper method that adds species to the CRN and automatically looks up
        and assigns their initial concentrations.

        Parameters
        ----------
        new_species : Species or list of Species
            Species to add to the CRN.
        component : Component, optional
            The component that generated these species. Used for
            component-specific initial concentration lookup.
        no_initial_concentrations : bool, default=False
            If True, skips initial concentration lookup and assignment.
        copy_species : bool, default=True
            If True, deep-copies species before adding them to the CRN.
        compartment : Compartment, optional
            Compartment to assign to the species. Overrides species'
            existing compartments.

        Returns
        -------
        list of Species
            All species in the CRN after addition (may include pre-existing
            species).

        Notes
        -----
        This method tracks which species are newly added and only assigns
        initial concentrations to those new species, preventing overwriting
        of previously set initial concentrations.

        """
        if self.crn is None:
            self.crn = ChemicalReactionNetwork(species=[], reactions=[])

        if isinstance(new_species, Species):
            new_species = [new_species]
        # 1) snapshot what’s already in the CRN
        before = set(self.crn._species_set)

        # 2) add (deep-copied) new_species into the CRN, with
        # compartment override
        all_species = self.crn.add_species(
            new_species, copy_species=copy_species, compartment=compartment
        )

        # 3) compute which Species objects we just inserted
        after = set(self.crn._species_set)
        just_added = list(after - before)

        # 4) assign initial concentrations only on those new copies
        if not no_initial_concentrations and just_added:
            ic_dict = self.get_initial_concentration(just_added, component)
            # merge into any existing initial_concentration_dict
            existing = (
                getattr(self.crn, 'initial_concentration_dict', {}) or {}
            )
            existing.update(ic_dict)
            self.crn.initial_concentration_dict = existing

        # 5) return the updated species list as before
        return all_species
        # self.crn.add_species(new_species, copy_species = copy_species,
        #                      compartment = compartment)

        # if not no_initial_concentrations:
        #     init_conc_dict = self.get_initial_concentration(
        #         remove_bindloc(new_species), component)
        #     self.crn.initial_concentration_dict = init_conc_dict

    def apply_global_mechanisms(
        self, species, compartment=None
    ) -> Tuple[List[Species], List[Reaction]]:
        """Apply all global mechanisms to a set of species.

        Calls each global mechanism's `update_species_global` and
        `update_reactions_global` methods, then adds the resulting species
        and reactions to the CRN.

        Parameters
        ----------
        species : list of Species
            Species to which global mechanisms should be applied.
        compartment : Compartment, optional
            Compartment for newly generated species and reactions.

        Returns
        -------
        tuple of (list of Species, list of Reaction)
            New species and reactions generated by global mechanisms.

        Notes
        -----
        Global mechanisms are typically used for operations that affect all
        species uniformly, such as dilution, global degradation, or
        compartment transport.

        """
        global_mech_species = []
        global_mech_reactions = []
        if self.global_mechanisms:
            for mech in self.global_mechanisms:
                # Update Global Mechanisms
                global_mech_species += self.global_mechanisms[
                    mech
                ].update_species_global(species, self, compartment)
                global_mech_reactions += self.global_mechanisms[
                    mech
                ].update_reactions_global(species, self, compartment)

        self.add_species_to_crn(
            global_mech_species, component=None, compartment=compartment
        )
        self.crn.add_reactions(global_mech_reactions, compartment=compartment)
        return global_mech_species, global_mech_reactions

    def component_enumeration(
        self, comps_to_enumerate=None, recursion_depth=10
    ) -> List[Component]:
        """Recursively enumerate components to generate derived components.

        Calls each component's `enumerate_components` method repeatedly to
        expand high-level components into their constituent parts (e.g.,
        DNAconstruct --> RNAconstruct --> Protein).

        Parameters
        ----------
        comps_to_enumerate : list of Component, optional
            Initial components to enumerate. If None, uses all components in
            the mixture.
        recursion_depth : int, default=10
            Maximum number of enumeration iterations. Prevents infinite
            recursion.

        Returns
        -------
        list of Component
            All components including the original components and all derived
            components generated through enumeration.

        Warns
        -----
        UserWarning
            Warns if unenumerated components remain after reaching the
            recursion depth limit.

        """
        all_components = []
        new_components = []
        if comps_to_enumerate is None:
            comps_to_enumerate = list(self.components)

        if recursion_depth == 0:
            all_components = comps_to_enumerate
        else:
            # Recursion depth
            for a in range(recursion_depth):
                for component in comps_to_enumerate:
                    component.set_mixture(self)
                    component.compartment = self.compartment
                    enumerated = component.enumerate_components(
                        previously_enumerated=all_components + new_components
                    )
                    new_components += enumerated

                all_components += comps_to_enumerate
                comps_to_enumerate = list(new_components)
                new_components = []

        if len(comps_to_enumerate) > 0:
            warn(
                "Mixture was left with unenumerated components "
                + str(', '.join([str(c) for c in comps_to_enumerate]))
            )
        return all_components

    def global_component_enumeration(
        self, comps_to_enumerate=None, recursion_depth=None
    ) -> List[Component]:
        """Apply global component enumerators to generate new components.

        Global component enumerators create new components based on patterns
        across all components (e.g., generating all pairwise binding
        interactions between proteins).

        Parameters
        ----------
        comps_to_enumerate : list of Component, optional
            Initial components to pass to enumerators. If None, uses all
            components in the mixture.
        recursion_depth : int, optional
            Maximum number of enumeration iterations. If None, uses
            `self.global_recursion_depth`.

        Returns
        -------
        list of Component
            All components including original and newly generated components
            from global enumeration.

        Notes
        -----
        This method is called during `compile_crn` before local component
        enumeration. Global enumerators are useful for creating complex
        interaction networks without manually specifying every interaction.

        """
        if recursion_depth is None:
            recursion_depth = self.global_recursion_depth

        if comps_to_enumerate is None:
            # these go into the ComponentEnuemrators
            comps_to_enumerate = list(self.components)

        # Recursion depth
        enumerated_components = list(
            comps_to_enumerate
        )  # These will be returned
        for global_enumerator in self.global_component_enumerators:
            # reset the enumeration if there's any stored info
            global_enumerator.reset(enumerated_components)
        for a in range(recursion_depth):
            # these will be added to comps_to_enumerate at the end of
            # the iteration
            new_comps_to_enumerate = []

            for global_enumerator in self.global_component_enumerators:
                enumerated = global_enumerator.enumerate_components(
                    comps_to_enumerate,
                    # this should be only NEWLY CREATED components
                    previously_enumerated=enumerated_components,
                )
                for c in enumerated:
                    # These components are passed into the enumerator
                    # next recursion
                    if (
                        c not in new_comps_to_enumerate
                        and c not in comps_to_enumerate
                    ):
                        new_comps_to_enumerate.append(c)

                    # These components are returend
                    if c not in enumerated_components:
                        enumerated_components.append(c)

            # Update comps_to_enumerate
            comps_to_enumerate += new_comps_to_enumerate

        return enumerated_components

    def compile_crn(
        self,
        recursion_depth: int = None,
        initial_concentration_dict: dict = None,
        return_enumerated_components: bool = False,
        initial_concentrations_at_end: bool = False,
        copy_objects: bool = True,
        add_reaction_species: bool = True,
        compartment: Compartment = None,
    ) -> ChemicalReactionNetwork:
        """Compile a chemical reaction network from the mixture.

        Enumerates components, generates species and reactions from each
        component, applies global mechanisms, and returns a complete CRN.

        Parameters
        ----------
        recursion_depth : int, optional
            Maximum recursion depth for both local and global component
            enumeration. If None, uses `self.global_recursion_depth`.
        initial_concentration_dict : dict, optional
            Dictionary mapping species to initial concentrations. This
            overrides all other initial concentration settings and is applied
            at the very end of compilation.
        return_enumerated_components : bool, default=False
            If True, returns a tuple of (CRN, enumerated_components) instead
            of just the CRN.
        initial_concentrations_at_end : bool, default=False
            If True, initial concentrations are only set at the end using the
            mixture's parameter database, ignoring component-specific initial
            concentrations during compilation.
        copy_objects : bool, default=True
            If True, species and reactions are deep-copied when added to the
            CRN. Protects CRN validity at the expense of compilation speed.
        add_reaction_species : bool, default=True
            If True, species appearing in reactions are automatically added to
            the CRN. Ensures no missing species at the expense of compilation
            speed.
        compartment : Compartment, optional
            Compartment to assign to all species and reactions that are not
            already assigned to a compartment. If None, uses
            `self.compartment`.

        Returns
        -------
        ChemicalReactionNetwork or tuple
            If `return_enumerated_components` is False, returns the compiled
            `ChemicalReactionNetwork`. If True, returns a tuple of
            (ChemicalReactionNetwork, list of enumerated Components).

        Notes
        -----
        The compilation process follows these steps:

        1. Add any directly-added species to the CRN
        2. Global component enumeration (generates component interactions)
        3. Local component enumeration (e.g., DNA --> RNA --> Protein)
        4. Generate species from all enumerated components
        5. Generate reactions from all enumerated components
        6. Apply global mechanisms to all species
        7. Set initial concentrations

        Examples
        --------
        Basic compilation:

        >>> gene = bcp.DNAassembly(
        ...     'GFP', promoter='pconst', rbs='RBS', protein='GFP')
        >>> mixture = bcp.Mixture(
        ...     name="txtl_extract",
        ...     components=[gene],
        ...     mechanisms={
        ...         'transcription': bcp.SimpleTranscription(),
        ...         'translation': bcp.SimpleTranslation()
        ...     },
        ...     parameters={'ktx': 0.05, 'ktl': 0.01}
        ... )
        >>> crn = mixture.compile_crn()

        Compilation with custom initial concentrations:

        >>> crn = mixture.compile_crn(
        ...     initial_concentration_dict={gene.dna: 1, gene.transcript: 50}
        ... )

        Get both CRN and enumerated components:

        >>> crn, components = mixture.compile_crn(
        ...     return_enumerated_components=True
        ... )

        """
        resetwarnings()

        if compartment is None and hasattr(self, 'compartment'):
            compartment = self.compartment

        self.crn = ChemicalReactionNetwork([], [])

        # helper to flatten & drop duplicates via ==, retagging
        # default -> compartment
        def _filter_new_by_eq(cands, existing_list):
            # flatten nested lists
            flat = []
            for x in cands:
                if isinstance(x, (list, tuple)):
                    flat.extend(x)
                else:
                    flat.append(x)

            new = []
            for s in flat:
                # retag any default compartment on the candidate itself
                if (
                    compartment
                    and getattr(s, 'compartment', None)
                    and s.compartment.name == 'default'
                ):
                    s.compartment = compartment

                # use __eq__ to check if already in existing_list
                if not any(s == e for e in existing_list):
                    new.append(s)
                    existing_list.append(s)
            return new

        # add the extra species
        if self.added_species:
            self.add_species_to_crn(
                self.added_species,
                component=None,
                no_initial_concentrations=initial_concentrations_at_end,
                copy_species=copy_objects,
                compartment=compartment,
            )

        # enumerate components
        if recursion_depth is None:
            recursion_depth = self.global_recursion_depth
        global_comps = self.global_component_enumeration(
            recursion_depth=recursion_depth
        )
        enumerated = self.component_enumeration(
            global_comps, recursion_depth=self.local_recursion_depth
        )

        for comp in enumerated:
            comp.set_mixture(self)
            if compartment is not None:
                comp.compartment = compartment

        # component species
        for comp in enumerated:
            sp_list = comp.update_species()
            # Use _filter_new_by_eq to prevent duplicates, just like
            # for reaction species
            existing = list(self.crn._species_set)
            fresh = _filter_new_by_eq(sp_list, existing)
            if fresh:
                self.add_species_to_crn(
                    fresh,
                    component=comp,
                    no_initial_concentrations=initial_concentrations_at_end,
                    copy_species=copy_objects,
                    compartment=compartment,
                )

        # component reactions & their species
        for comp in enumerated:
            rxns = comp.update_reactions()
            self.crn.add_reactions(
                rxns,
                copy_reactions=copy_objects,
                add_species=False,
                compartment=compartment,
            )

            # retag species with compartments inside these new reactions
            new_rxns = self.crn._reactions[-len(rxns) :]
            for r in new_rxns:
                for w in r.inputs + r.outputs:
                    if (
                        compartment
                        and w.species.compartment.name == 'default'
                    ):
                        w.species.compartment = compartment

            # collect & add only new species
            if add_reaction_species:
                participants = [
                    w.species
                    for r in new_rxns
                    for w in (r.inputs + r.outputs)
                ]
                existing = list(self.crn._species_set)
                fresh = _filter_new_by_eq(participants, existing)
                if fresh:
                    self.add_species_to_crn(
                        fresh,
                        component=comp,
                        no_initial_concentrations=initial_concentrations_at_end,
                        copy_species=copy_objects,
                        compartment=compartment,
                    )

        # global mechanisms
        global_sp, global_rxns = [], []
        for mech in self.global_mechanisms.values():
            global_sp.extend(
                mech.update_species_global(list(self.crn._species_set), self)
            )
            global_rxns.extend(
                mech.update_reactions_global(
                    list(self.crn._species_set), self
                )
            )

        if global_sp:
            existing = list(self.crn._species_set)
            fresh = _filter_new_by_eq(global_sp, existing)
            if fresh:
                self.add_species_to_crn(
                    fresh,
                    component=None,
                    no_initial_concentrations=initial_concentrations_at_end,
                    copy_species=copy_objects,
                    compartment=compartment,
                )

        if global_rxns:
            self.crn.add_reactions(
                global_rxns,
                copy_reactions=copy_objects,
                add_species=False,
                compartment=compartment,
            )
            new_grx = self.crn._reactions[-len(global_rxns) :]
            for r in new_grx:
                for w in r.inputs + r.outputs:
                    if (
                        compartment
                        and w.species.compartment.name == 'default'
                    ):
                        w.species.compartment = compartment

            if add_reaction_species:
                participants = [
                    w.species for r in new_grx for w in (r.inputs + r.outputs)
                ]
                existing = list(self.crn._species_set)
                fresh = _filter_new_by_eq(participants, existing)
                if fresh:
                    self.add_species_to_crn(
                        fresh,
                        component=None,
                        no_initial_concentrations=initial_concentrations_at_end,
                        copy_species=copy_objects,
                        compartment=compartment,
                    )

        # initial‐conc handling
        if initial_concentrations_at_end:
            self.crn.initial_concentration_dict = (
                self.get_initial_concentration(
                    list(self.crn._species_set), component=None
                )
            )
        if initial_concentration_dict is not None:
            self.crn.initial_concentration_dict = initial_concentration_dict

        # return
        if return_enumerated_components:
            return self.crn, enumerated
        return self.crn

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self) + '\n'
        if self.components:
            txt += 'Components = ['
            for comp in self.components:
                txt += '\n\t' + str(comp)
        if self.mechanisms:
            txt += ' ]\nMechanisms = {'
            for mech in self.mechanisms:
                txt += '\n\t' + mech + ':' + self.mechanisms[mech].name
        if self.global_mechanisms:
            txt += ' }\nGlobal Mechanisms = {'
            for mech in self.global_mechanisms:
                txt += '\n\t' + mech + ':' + self.global_mechanisms[mech].name
        txt += ' }'
        return txt
