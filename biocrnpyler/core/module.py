# module.py - grouping of components, mechanisms and parameters
# RMM, 16 Aug 2026

import copy
from typing import List, Union
from warnings import warn

from .component import Component
from .mixture import Mixture
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
    to the same species. This is how modules are wired together: a
    protein produced in one module and consumed in another is written
    with the same name in both.

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
    def _matching_component(component, components):
        """Return the component in `components` matching type and name."""
        for comp in components:
            if type(comp) is type(component) and comp.name == component.name:
                return comp
        return None

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

        Parameters are added only for keys the component does not
        already define, so a component's own parameters take precedence
        over the module's.

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
        for entry in entries:
            if entry.parameter_key not in database.parameters:
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

        Parameters
        ----------
        previously_enumerated : set or list, optional
            Collection of components already enumerated. Unused by
            Module, which always enumerates into its members.

        Returns
        -------
        list of Component
            Copies of the module's components, each carrying the
            module's mechanisms and parameters.

        See Also
        --------
        apply_context : Place a single component in the module's context.

        """
        return [self.apply_context(comp) for comp in self.components]

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
