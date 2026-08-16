# module.py - grouping of components, mechanisms and parameters
# RMM, 16 Aug 2026

import copy
from typing import List, Union
from warnings import warn

from .component import Component
from .mixture import Mixture
from .parameter import ParameterKey
from .species import Species


class Module(Component):
    """A reusable group of components, mechanisms, and parameters.

    A Module bundles together the components that make up a subsystem
    (e.g. a signaling pathway or a genetic circuit) along with the
    mechanisms and parameters that subsystem needs. Modules can be
    combined with each other and added to a `Mixture` using ``+``::

        system = mixture + module1 + module2

    A Module is itself a `Component`, so it can be placed anywhere a
    component can go, including inside another Module. During
    compilation a Module contributes no reactions of its own: it
    enumerates into its member components, each of which receives the
    Module's mechanisms and parameters for any mechanism type or
    parameter key it does not already define itself.

    This gives a three-level lookup hierarchy for the components inside
    a Module:

    1. the component's own mechanisms and parameters
    2. the Module's mechanisms and parameters
    3. the Mixture's mechanisms and parameters

    so two Modules in the same Mixture may use different mechanisms for
    the same mechanism type.

    Parameters
    ----------
    name : str
        Name of the module, used for identification and parameter lookup.
    components : list of Component or Component, optional
        Components belonging to this module. Components are deep-copied
        when added to prevent modification of the original objects.
    mechanisms : dict, list, or Mechanism, optional
        Mechanisms used by the components in this module. Can be a dict
        with mechanism types (str) as keys and mechanism objects as
        values, or a list of mechanism objects. These take precedence
        over the Mixture's mechanisms but not over mechanisms defined on
        an individual component.
    parameters : dict, optional
        Dictionary of parameter values. Keys follow the format
        (mechanism, part_id, param_name).
    parameter_file : str, optional
        Path to a CSV or TSV file containing parameters to load.
    species : list of Species or Species, optional
        Species added directly to the CRN by this module, without going
        through component compilation.
    **keywords
        Additional keyword arguments passed to `Component`.

    Attributes
    ----------
    name : str
        Name of the module.
    components : list of Component
        Member components of the module (deep copies of those added).
    mechanisms : dict
        Dictionary of the module's mechanisms, keyed by mechanism type.
    parameter_database : ParameterDatabase
        Database storing the module's parameters.
    species : list of Species
        Species added directly by this module.

    See Also
    --------
    Component : Base class for biomolecular components.
    Mixture : The context in which modules and components are compiled.

    Notes
    -----
    Unlike a `Mixture`, a Module does not hold global mechanisms and
    cannot be compiled on its own; it must be placed in a Mixture::

        crn = Mixture('test', components=[my_module]).compile_crn()

    Components in different modules that share a type and name compile
    to the same species, and are compiled only once. This is how
    modules are wired together: a protein produced in one module and
    consumed in another is written with the same name in both.

    A module can be used more than once in the same mixture by making
    copies of it with `instance`, which renames the components that
    should differ between the copies and leaves the rest shared.

    Precedence goes by level rather than by how specifically a
    parameter is keyed: a parameter written on a component beats one
    from its module, which beats one from the mixture, even if the
    outer one names the parameter more precisely. Mechanisms follow the
    same order. Note that this means an outer level cannot override an
    inner one; to change a parameter, set it where it was defined, with
    `update_parameters` and ``overwrite_parameters=True``.

    Examples
    --------
    Group components into a module and add it to a mixture:

    >>> reporter = bcp.Module(
    ...     name='reporter',
    ...     components=[bcp.DNAassembly('GFP', promoter='P', rbs='B')],
    ...     mechanisms={'transcription': bcp.SimpleTranscription(),
    ...                 'translation': bcp.SimpleTranslation()},
    ...     parameters={'ktx': 0.05, 'ktl': 0.01},
    ... )
    >>> system = bcp.Mixture('extract') + reporter
    >>> crn = system.compile_crn()

    Combine two modules before adding them to a mixture:

    >>> system = bcp.Mixture('extract') + (reporter + sensor)

    Override a mechanism on one module without affecting the others:

    >>> reporter.add_mechanism(bcp.SimpleTranscription(), overwrite=True)

    Use the same module twice, with different inputs and outputs:

    >>> system = (bcp.Mixture('extract')
    ...           + signaling.instance('s1', {'ligand': 'IPTG'})
    ...           + signaling.instance('s2', {'ligand': 'aTc'}))

    """

    def __init__(
        self,
        name: str,
        components=None,
        mechanisms=None,
        parameters=None,
        parameter_file=None,
        species=None,
        **keywords,
    ):
        super().__init__(
            name=name,
            mechanisms=mechanisms,
            parameters=parameters,
            parameter_file=parameter_file,
            **keywords,
        )

        self.components = []
        if components is not None:
            self.add_components(components)

        self.species = []
        self.add_species(species)

    def add_component(self, component: Component):
        """Add a single component to the module.

        Parameters
        ----------
        component : Component or list of Component
            Component to add to the module. If a list is provided, calls
            `add_components` instead. The component is deep-copied before
            being added.

        Raises
        ------
        ValueError
            If `component` is not a Component, or if a component of the
            same type and name is already in the module.

        """
        if isinstance(component, list):
            self.add_components(component)
            return

        if not isinstance(component, Component):
            raise ValueError(
                f"add_component expected a Component. Received {component}."
            )

        if self._matching_component(component, self.components) is not None:
            raise ValueError(
                f"{component} of same type and name already in Module!"
            )

        # Components are copied before being added to Modules
        self.components.append(copy.deepcopy(component))

    def add_components(self, components: Union[List[Component], Component]):
        """Add multiple components to the module.

        Parameters
        ----------
        components : Component or list of Component
            Component(s) to add to the module. Each component is
            deep-copied before being added.

        Raises
        ------
        ValueError
            If `components` is not a Component or a list of Components.

        See Also
        --------
        add_component : Add a single component to the module.

        """
        if isinstance(components, Component):
            self.add_component(components)
        elif isinstance(components, list):
            for component in components:
                self.add_component(component)
        else:
            raise ValueError(
                f"add_components expected a list of Components. "
                f"Received {components}."
            )

    def add_species(self, species: Union[List[Species], Species]):
        """Add species directly to the module without component compilation.

        Parameters
        ----------
        species : Species or list of Species
            Species object(s) added to the CRN when this module is
            compiled.

        Raises
        ------
        ValueError
            If any element is not a Species object.

        """
        if species is None:
            return

        if isinstance(species, Species):
            species = [species]

        if not all(isinstance(s, Species) for s in species):
            raise ValueError(
                f"add_species expected a list of Species. Received {species}."
            )

        self.species += species

    @staticmethod
    def _duplicate_component(component, components):
        """Return the component in `components` that would compile the same.

        Matches on type and name like `_matching_component`, but also
        requires the same species, so that two components which share a
        name but produce different species are both compiled.

        Notes
        -----
        Some components derive their name from a species and the
        compartment it is in, and keep the name they were built with
        (`IntegralMembraneProtein` is one). A renamed copy of such a
        component therefore still carries the original's name while
        producing different species, and matching on name alone would
        drop it from the model.

        """
        for comp in Module._all_matching_components(component, components):
            if comp.get_species() == component.get_species():
                return comp
        return None

    @staticmethod
    def _matching_component(component, components):
        """Return the component in `components` matching type and name."""
        for comp in Module._all_matching_components(component, components):
            return comp
        return None

    @staticmethod
    def _all_matching_components(component, components):
        """Yield the components matching `component` by type and name."""
        for comp in components:
            if type(comp) is type(component) and comp.name == component.name:
                yield comp

    def apply_context(self, component: Component) -> Component:
        """Return a copy of a component carrying this module's context.

        The module's mechanisms and parameters are copied onto the
        component for every mechanism type and parameter key the
        component does not already define itself. Because a Component
        searches its own mechanisms and parameters before those of its
        Mixture, this gives the module's settings precedence over the
        Mixture's while leaving the component's own settings untouched.

        Parameters
        ----------
        component : Component
            The component to place in this module's context.

        Returns
        -------
        Component
            A deep copy of `component` with the module's mechanisms and
            parameters filled in.

        See Also
        --------
        enumerate_components : Applies this to every component in the module.

        """
        new_component = copy.deepcopy(component)

        # Module mechanisms fill in the ones the component does not define.
        # add_mechanism is overridden by components that hold sub-components
        # (e.g. DNAassembly), so this reaches those sub-components too.
        for mech_type, mech in self.mechanisms.items():
            new_component.add_mechanism(
                mech, mech_type, overwrite=False, optional_mechanism=True
            )

        # Module parameters fill in the ones the component does not define
        self._merge_parameters(
            new_component, list(self.parameter_database), set()
        )

        return new_component

    @staticmethod
    def _merge_parameters(component, entries, visited):
        """Add parameters to a component and its sub-components.

        A parameter is added only if the component cannot already answer
        that lookup for itself, so a component's own parameters take
        precedence over its module's.

        Parameters
        ----------
        component : Component
            The component to add parameters to.
        entries : list of ParameterEntry
            The parameter entries to add.
        visited : set of int
            Object ids already visited, used to stop recursion when
            sub-components refer back to their parent.

        Notes
        -----
        Whether a component can answer a lookup is decided with
        `ParameterDatabase.find_parameter`, which falls back through
        progressively more general keys, rather than by comparing keys
        exactly. A component's loosely keyed parameter therefore takes
        precedence over a specifically keyed one from its module: once
        the module's parameters are copied in they sit in the same
        database as the component's own, where the more specific key
        would otherwise win. Deciding by lookup keeps the precedence
        between a component and its module the same whatever keys the
        two happen to use, and matches how mechanisms behave.

        Components that build sub-components (a `DNAassembly` builds a
        `Promoter` and an `RBS`, for example) copy their parameters to
        those sub-components when they are constructed, which happens
        before a Module applies its context. Since the sub-components
        are what resolve parameters during compilation, the module's
        parameters have to be added to them directly.

        The recursion stops at a nested Module, which passes parameters
        to its own components when it is enumerated. Stopping here keeps
        the nested Module's parameters ahead of the enclosing one's.

        The entries themselves are shared rather than copied, matching
        `ParameterDatabase.load_parameters_from_database`, so units and
        parameter origins are preserved.

        """
        if id(component) in visited:
            return
        visited.add(id(component))

        database = component.parameter_database

        # Decide everything against the component's own parameters,
        # before any of the module's have been added: an entry added
        # early would otherwise answer the lookup for a later one and
        # keep it out.
        to_add = [
            entry
            for entry in entries
            if database.find_parameter(
                entry.parameter_key.mechanism,
                entry.parameter_key.part_id,
                entry.parameter_key.name,
            )
            is None
        ]
        for entry in to_add:
            database[entry.parameter_key] = entry

        if isinstance(component, Module):
            return

        for value in vars(component).values():
            if isinstance(value, Component):
                Module._merge_parameters(value, entries, visited)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, Component):
                        Module._merge_parameters(item, entries, visited)

    def enumerate_components(self, previously_enumerated=None) -> List:
        """Enumerate the module into its member components.

        A component is skipped if one of the same type and name has
        already been enumerated, or is already a component of the
        mixture. Modules share components by name, so without this the
        shared component would be compiled once per module and its
        reactions would appear in the CRN more than once.

        Parameters
        ----------
        previously_enumerated : set or list, optional
            Components already enumerated, used to compile a component
            shared by several modules only once.

        Returns
        -------
        list of Component
            Copies of the module's components, each carrying the
            module's mechanisms and parameters.

        Notes
        -----
        A component declared by more than one module is compiled under
        the context of whichever module is enumerated first. If two
        modules need different mechanisms or parameters for it, rename
        it in one of them with `instance` so that they are separate
        components.

        See Also
        --------
        apply_context : Place a single component in the module's context.
        instance : Create a copy of a module with components renamed.

        """
        already_enumerated = list(previously_enumerated or [])
        if self.mixture is not None:
            already_enumerated += [
                comp for comp in self.mixture.components if comp is not self
            ]

        enumerated = []
        for component in self.components:
            if (
                self._duplicate_component(component, already_enumerated)
                is not None
            ):
                continue
            contextualized = self.apply_context(component)
            enumerated.append(contextualized)
            already_enumerated.append(contextualized)

        return enumerated

    def instance(self, name: str, rename=None):
        """Create a named copy of this module with components renamed.

        Used to place more than one copy of the same module in a
        mixture. Only the names given in `rename` change; anything left
        alone keeps its name and is therefore shared between the
        instances, which is how the copies are wired to common parts of
        the system.

        Renaming reaches components, the species they hold, and
        parameter keys that refer to them by name, wherever these
        appear in the module, including inside nested modules.

        Parameters
        ----------
        name : str
            Name of the new module. Instances need distinct names, since
            a Mixture holds only one component of each type and name.
        rename : dict, optional
            Dictionary mapping current names to new names.

        Returns
        -------
        Module
            A copy of this module under the new name, with the renaming
            applied. This module is left unchanged.

        Raises
        ------
        ValueError
            If `rename` is not a dictionary of strings, or if any name
            in it is not found in the module. An unused entry usually
            means a name was mistyped, which would silently leave the
            instances sharing a component they were meant to separate.

        See Also
        --------
        enumerate_components : Where shared components are compiled once.

        Examples
        --------
        Put two copies of a signaling module in one mixture:

        >>> signaling = bcp.Module(
        ...     'signaling', components=[receptor, kinase])
        >>> s1 = signaling.instance(
        ...     's1', rename={'ligand': 'IPTG', 'output': 'GFP'})
        >>> s2 = signaling.instance(
        ...     's2', rename={'ligand': 'aTc', 'output': 'RFP'})
        >>> system = bcp.BasicPURE() + s1 + s2

        The kinase is not renamed, so both instances share it.

        """
        if rename is None:
            rename = {}
        if not isinstance(rename, dict):
            raise ValueError(
                f"rename must be a dictionary of names. Received {rename}."
            )
        for old_name, new_name in rename.items():
            if not isinstance(old_name, str) or not isinstance(new_name, str):
                raise ValueError(
                    f"rename must map strings to strings. Received "
                    f"{old_name}: {new_name}."
                )

        new_module = copy.deepcopy(self)
        new_module.name = name

        renamed = set()
        visited = set()
        for component in new_module.components:
            self._rename(component, rename, visited, renamed)
        for species in new_module.species:
            self._rename(species, rename, visited, renamed)
        self._rename_parameter_keys(
            new_module.parameter_database, rename, renamed
        )

        unused = sorted(set(rename) - renamed)
        if unused:
            raise ValueError(
                f"{unused} not found in Module {self.name}, so nothing was "
                f"renamed for them. Check the names in rename."
            )

        return new_module

    @staticmethod
    def _rename(obj, rename, visited, renamed):
        """Rename a component or species, and everything it holds.

        Parameters
        ----------
        obj : Component or Species
            The object to rename.
        rename : dict
            Dictionary mapping current names to new names.
        visited : set of int
            Object ids already visited. Stops the recursion when
            sub-components refer back to their parent, and makes sure
            each object is renamed only once, so that a rename such as
            ``{'a': 'b', 'b': 'c'}`` does not rename anything twice.
        renamed : set of str
            Collects the names that were found, so that unused entries
            in `rename` can be reported.

        Notes
        -----
        Unlike `_merge_parameters` this recurses through nested Modules,
        since renaming is not a question of precedence: an instance
        renames everything it contains.

        """
        if id(obj) in visited:
            return
        visited.add(id(obj))

        if obj.name in rename:
            renamed.add(obj.name)
            obj.name = rename[obj.name]

        if isinstance(obj, Component):
            Module._rename_parameter_keys(
                obj.parameter_database, rename, renamed
            )

        # Species hold the species they are built from, so this reaches
        # into complexes as well as into sub-components
        for value in vars(obj).values():
            if isinstance(value, (Component, Species)):
                Module._rename(value, rename, visited, renamed)
            elif isinstance(value, list):
                for item in value:
                    if isinstance(item, (Component, Species)):
                        Module._rename(item, rename, visited, renamed)

    @staticmethod
    def _rename_parameter_keys(database, rename, renamed):
        """Rename the parameters that refer to a renamed component.

        Parameters
        ----------
        database : ParameterDatabase
            The database to rewrite the keys of.
        rename : dict
            Dictionary mapping current names to new names.
        renamed : set of str
            Collects the names that were found.

        Notes
        -----
        Parameters are looked up by `part_id` and by name, and a
        component's initial concentration is stored under the
        component's name, so both fields have to follow a rename.

        """
        new_parameters = {}
        for key, entry in database.parameters.items():
            new_key = ParameterKey(
                mechanism=key.mechanism,
                part_id=rename.get(key.part_id, key.part_id),
                name=rename.get(key.name, key.name),
            )
            if new_key != key:
                renamed.update({key.part_id, key.name} & set(rename))
                entry.parameter_key = new_key
            new_parameters[new_key] = entry

        database.parameters = new_parameters

    def update_species(self) -> List[Species]:
        """Return the species added directly by this module.

        Returns
        -------
        list of Species
            The module's directly added species. Species belonging to the
            module's components are generated by those components during
            enumeration.

        """
        return list(self.species)

    def update_reactions(self) -> List:
        """Return the reactions produced by the module itself.

        Returns
        -------
        list
            Always empty. Reactions come from the module's components,
            which are enumerated separately.

        """
        return []

    def __add__(self, other):
        """Combine two modules into a new module.

        The mechanisms and parameters of each module are applied to that
        module's own components before merging, so components keep the
        settings of the module they came from.

        Parameters
        ----------
        other : Module
            The module to combine with this one.

        Returns
        -------
        Module
            A new module containing the components and species of both.

        Warns
        -----
        UserWarning
            If both modules contain a component of the same type and
            name. The component from the left module is kept.

        """
        if not isinstance(other, Module):
            return NotImplemented

        new_module = Module(name=f"{self.name}_{other.name}")
        for module in [self, other]:
            for component in module.components:
                match = self._matching_component(
                    component, new_module.components
                )
                if match is not None:
                    warn(
                        f"{component} is in both Module {self.name} and "
                        f"Module {other.name}. The copy from Module "
                        f"{self.name} has been kept."
                    )
                    continue
                new_module.add_component(module.apply_context(component))
            new_module.add_species(module.species)

        return new_module

    def __radd__(self, other):
        """Add this module to a mixture, returning a new mixture.

        Supports the ``mixture + module`` syntax. `Mixture` does not
        define ``+``, so Python falls back to this method.

        Parameters
        ----------
        other : Mixture
            The mixture to add this module to. It is not modified.

        Returns
        -------
        Mixture
            A copy of `other` with this module added as a component.

        """
        if not isinstance(other, Mixture):
            return NotImplemented

        new_mixture = copy.deepcopy(other)
        # Any CRN compiled from the old mixture does not describe the new one
        new_mixture.crn = None
        new_mixture.add_component(self)
        return new_mixture

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self)
        if self.components:
            txt += '\nComponents = ['
            for comp in self.components:
                txt += '\n\t' + str(comp)
            txt += ' ]'
        if self.mechanisms:
            txt += '\nMechanisms = {'
            for mech in self.mechanisms:
                txt += '\n\t' + mech + ':' + self.mechanisms[mech].name
            txt += ' }'
        if self.species:
            txt += '\nSpecies = ['
            for spec in self.species:
                txt += '\n\t' + str(spec)
            txt += ' ]'
        return txt
