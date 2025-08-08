
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from typing import List, Union, Tuple
from warnings import resetwarnings, warn

from .chemical_reaction_network import ChemicalReactionNetwork
from .component import Component
from ..mechanisms.global_mechanisms import GlobalMechanism
from .mechanism import Mechanism
from .parameter import ParameterDatabase
from .reaction import Reaction
from .species import Species
from .compartment import Compartment
from ..utils.general import remove_bindloc


class Mixture(object):
    def __init__(self, name="", mechanisms=None, components=None, parameters=None, compartment: Compartment = None,
                 parameter_file=None, global_mechanisms=None, species=None, initial_condition_dictionary=None,
                 global_component_enumerators=None, global_recursion_depth=4, local_recursion_depth=None, **kwargs):
        """A Mixture object holds together all the components (DNA,Protein, etc), mechanisms (Transcription, Translation),
        and parameters related to the mixture itself (e.g. Transcription rate). Default components and mechanisms can be
        added as well as global mechanisms that impacts all species (e.g. cell growth).

        :param name: Name of the mixture
        :param mechanisms: Dictionary of mechanisms
        :param components: List of components in the mixture (list of Components)
        :param parameters: Dictionary of parameters (check parameters documentation for the keys)
        :param parameter_file: Parameters can be loaded from a parameter file
        :param default_mechanisms:
        :param global_mechanisms: dict of global mechanisms that impacts all species (e.g. cell growth)
        """
        # Initialize instance variables
        self.name = name  # Save the name of the mixture
        self.compartment = compartment

        # recursion depth for global component enumeration
        self.global_recursion_depth = global_recursion_depth

        if (local_recursion_depth is None):
            self.local_recursion_depth = self.global_recursion_depth+2

        #

        # process the components
        if components is None and not hasattr(self, "_components"):
            self.components = []
        else:
            self.add_components(components)

        # process mechanisms:
        if mechanisms is None and not hasattr(self, "_mechanisms"):
            self.mechanisms = {}
        else:
            self.add_mechanisms(mechanisms)

        # process global_mechanisms:
        # Global mechanisms are applied just once ALL species generated from
        # components inside a mixture
        # Global mechanisms should be used rarely, and with care. An example
        # usecase is degradation via dilution.
        if global_mechanisms is None and not hasattr(self, "_global_mechanisms"):
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
            parameter_file=parameter_file, parameter_dictionary=parameters, **kwargs)

        # CRN is stored here during compilation
        self.crn = None

    def add_species(self, species: Union[List[Species], Species]):
        if not hasattr(self, "added_species"):
            self.added_species = []

        if species is not None:
            if not isinstance(species, list):
                species_list = [species]
            else:
                species_list = species

            assert all(isinstance(x, Species)
                       for x in species_list), 'only Species type is accepted!'

            self.added_species += species_list

    def set_species(self, species: Union[Species, str], material_type=None, attributes=None):
        """Used to set internal species from strings, Species or Components

        :param species: name of a species or a species instance
        :param material_type: material type of a species as a string
        :param attributes: Species attribute
        :return: Species in the mixture
        """
        if isinstance(species, Species):
            return species
        elif isinstance(species, str):
            return Species(name=species, material_type=material_type, attributes=attributes)
        elif isinstance(species, Component) and species.get_species() is not None:
            return species.get_species()
        else:
            raise ValueError(
                "Invalid Species: string, chemical_reaction_network.Species or Component with implemented .get_species() required as input.")

    @property
    def components(self):
        return self._components

    @components.setter
    def components(self, components):
        self._components = []
        self.add_components(components)

    def add_component(self, component):
        """this function adds a single component to the mixture."""
        if not hasattr(self, "_components"):
            self.components = []

        if isinstance(component, list):
            self.add_components(component)
        else:
            assert isinstance(
                component, Component), "the object: %s passed into mixture as component must be of the class Component" % str(component)

            # Check if component is already in self._components
            for comp in self._components:
                if type(comp) == type(component) and comp.name == component.name:
                    raise ValueError(
                        f"{comp} of the same type and name already in Mixture!")
            else:
                # Components are copied before being added to Mixtures
                component_copy = copy.deepcopy(component)
                component_copy.set_mixture(self)
                component_copy.compartment = self.compartment
                self.components.append(component_copy)

    def get_mechanism(self, mechanism_type):
        """Searches the Mixture for a Mechanism of the correct type.

        If no Mechanism is found, None is returned.
        """
        if not isinstance(mechanism_type, str):
            raise TypeError(
                f"mechanism_type must be a string. Recievied {mechanism_type}.")

        if mechanism_type in self.mechanisms:
            return self.mechanisms[mechanism_type]
        else:
            return None

    def add_components(self, components: Union[List[Component], Component]):
        """This function adds a list of components to the mixture.
        """
        if isinstance(components, Component):
            self.add_component(components)
        elif isinstance(components, List):
            for component in components:
                self.add_component(component)
        else:
            raise ValueError(
                f"add_components expected a list of Components. Received {components}")

    def get_component(self, component=None, name=None, index=None):
        """Function to get components from Mixture._components.

        One of the 3 keywords must not be None.

        :param component: an instance of a component. Searches Mixture._components for a Component with the same type and name.
        :param name: str. Searches Mixture._components for a Component with the same name
        :param index: int. returns Mixture._components[index]
        :return: if nothing is found, returns None.
        """
        if [component, name, index].count(None) != 2:
            raise ValueError(
                f"get_component requires a single keyword. Received component={component}, name={name}, index={index}.")
        if not (isinstance(component, Component) or component is None):
            raise ValueError(
                f"component must be of type Component. Received {component}.")
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
                    if type(comp) == type(component) and comp.name == component.name:
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
                "get_component found multiple matching components. A list has been returned.")
            return matches

    @property
    def mechanisms(self):
        """mechanisms stores Mixture Mechanisms."""
        return self._mechanisms

    @mechanisms.setter
    def mechanisms(self, mechanisms):
        self._mechanisms = {}
        self.add_mechanisms(mechanisms, overwrite=True)

    def add_mechanism(self, mechanism, mech_type=None, overwrite=False):
        """adds a mechanism of type mech_type to the Mixture mechanism_dictionary.

        :param mechanism: a Mechanism instance
        :param mech_type: the type of mechanism. defaults to mechanism.mech_type if None
        :param overwrite: whether to overwrite existing mechanisms of the same type (default False)
        :return:
        """
        if not hasattr(self, "_mechanisms"):
            self._mechanisms = {}

        if not isinstance(mechanism, Mechanism):
            raise TypeError(
                f"mechanism must be a Mechanism. Received {mechanism}.")

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(
                f"mechanism keys must be strings. Received {mech_type}")

        if isinstance(mechanism, GlobalMechanism):
            self.add_global_mechanism(mechanism, mech_type, overwrite)
        elif isinstance(mechanism, Mechanism):
            if mech_type in self._mechanisms and not overwrite:
                raise ValueError(
                    f"mech_type {mech_type} already in Mixture {self}. To overwrite, use keyword overwrite = True.")
            else:
                self._mechanisms[mech_type] = copy.deepcopy(mechanism)

    def add_mechanisms(self, mechanisms, overwrite=False):
        """This function adds a list or dictionary of mechanisms to the mixture.

        Can take both GlobalMechanisms and Mechanisms

        :param mechanisms: a Mechanism instance
        :param overwrite: whether to overwrite existing mechanisms of the same type (default False)
        :return:
        """
        if isinstance(mechanisms, Mechanism):
            self.add_mechanism(mechanisms, overwrite=overwrite)
        elif isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_mechanism(
                    mechanisms[mech_type], mech_type, overwrite=overwrite)
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_mechanism(mech, overwrite=overwrite)
        else:
            raise ValueError(
                f"add_mechanisms expected a list of Mechanisms. Recieved {mechanisms}")

    @property
    def global_mechanisms(self):
        """global_mechanisms stores global Mechanisms in the Mixture."""
        return self._global_mechanisms

    @global_mechanisms.setter
    def global_mechanisms(self, mechanisms):
        self._global_mechanisms = {}
        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_global_mechanism(
                    mechanisms[mech_type], mech_type, overwrite=True)
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_global_mechanism(mech, overwrite=True)

    def add_global_mechanism(self, mechanism, mech_type=None, overwrite=False):
        """adds a mechanism of type mech_type to the Mixture global_mechanism dictonary.

        keywordS:
          mechanism: a Mechanism instance
          mech_type: the type of mechanism. defaults to mechanism.mech_type if None
          overwrite: whether to overwrite existing mechanisms of the same type (default False)
        """
        if not hasattr(self, "_global_mechanisms"):
            self._global_mechanisms = {}

        if not isinstance(mechanism, GlobalMechanism):
            raise TypeError(
                f"mechanism must be a GlobalMechanism. Recieved {mechanism}.")

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(
                f"mechanism keys must be strings. Recieved {mech_type}")

        if mech_type in self._mechanisms and not overwrite:
            raise ValueError(
                f"mech_type {mech_type} already in Mixture {self}. To overwrite, use keyword overwrite = True.")
        else:
            self._global_mechanisms[mech_type] = copy.deepcopy(mechanism)

    def update_parameters(self, parameter_file=None, parameters=None, overwrite_parameters=True):
        if parameter_file is not None:
            self.parameter_database.load_parameters_from_file(
                parameter_file, overwrite_parameters=overwrite_parameters)

        if parameters is not None:
            self.parameter_database.load_parameters_from_dictionary(
                parameters, overwrite_parameters=overwrite_parameters)

    def get_parameter(self, mechanism, part_id, param_name):
        param = self.parameter_database.find_parameter(
            mechanism, part_id, param_name)
        return param

    def get_initial_concentration(self, S: Union[List, Species], component=None):
        """
        Tries to find an initial condition of species s using the parameter hierarchy using the key:

        1. Searches Component's ParameterDatabase using the key:
            mechanisms = "initial concentration"
            part_id = mixture.name
            parameter_name = str(s)

            if s == component.get_species, also checks with parameter_name=component.name

        2. Searches the Mixture's ParameterDatabase using the key:
            mechanisms = "initial concentration"
            part_id = mixture.name
            parameter_name = str(s)

            if s == component.get_species, also checks with parameter_name=component.name

        3. Defaults to 0
        """
        if isinstance(S, Species):
            S = [S]

        # flatten the species list
        S = Species.flatten_list(S)

        init_conc_dict = {}
        for s in S:
            if not isinstance(s, Species):
                raise ValueError(
                    f"{s} is not a Species! Can only find initial concentration of a Species.")
            init_conc = None
            # 1 Check the component
            if component is not None:
                init_conc = component.get_parameter(param_name=str(s),
                                                    part_id=self.name,
                                                    mechanism="initial concentration",
                                                    check_mixture=False,
                                                    return_none=True)
                if init_conc is None:
                    init_conc = component.get_parameter(param_name=s.name,
                                                        part_id=self.name,
                                                        mechanism="initial concentration",
                                                        check_mixture=False,
                                                        return_none=True)
                if init_conc is None and component.get_species() == s:
                    init_conc = component.get_parameter(param_name=component.name,
                                                        part_id=self.name,
                                                        mechanism="initial concentration",
                                                        check_mixture=False,
                                                        return_none=True)

            # 2 Check self
            if init_conc is None:
                init_conc = self.get_parameter(param_name=str(s),
                                               part_id=self.name,
                                               mechanism="initial concentration")
                if init_conc is None:
                    init_conc = self.get_parameter(param_name=s.name,
                                                   part_id=self.name,
                                                   mechanism="initial concentration")
                if init_conc is None and component is not None and component.get_species() == s:
                    init_conc = self.get_parameter(param_name=component.name,
                                                   part_id=self.name,
                                                   mechanism="initial concentration")

            if init_conc is None:
                init_conc = 0

            init_conc_dict[s] = init_conc

        return init_conc_dict

    def add_species_to_crn(self, new_species, component=None,
                           no_initial_concentrations=False, copy_species=True,
                           compartment=None):

        if self.crn is None:
            self.crn = ChemicalReactionNetwork(species=[], reactions=[])

        if isinstance(new_species, Species):
            new_species = [new_species]
        # 1) snapshot what’s already in the CRN
        before = set(self.crn._species_set)

        # 2) add (deep-copied) new_species into the CRN, with compartment override
        all_species = self.crn.add_species(new_species,
                                           copy_species=copy_species,
                                           compartment=compartment)

        # 3) compute which Species objects we just inserted
        after = set(self.crn._species_set)
        just_added = list(after - before)

        # 4) assign initial concentrations only on those new copies
        if not no_initial_concentrations and just_added:
            ic_dict = self.get_initial_concentration(just_added, component)
            # merge into any existing initial_concentration_dict
            existing = getattr(
                self.crn, "initial_concentration_dict", {}) or {}
            existing.update(ic_dict)
            self.crn.initial_concentration_dict = existing

        # 5) return the updated species list as before
        return all_species
        # self.crn.add_species(new_species, copy_species = copy_species,
        #                      compartment = compartment)

        # if not no_initial_concentrations:
        #     init_conc_dict = self.get_initial_concentration(remove_bindloc(new_species), component)
        #     self.crn.initial_concentration_dict = init_conc_dict

    def apply_global_mechanisms(self, species, compartment=None) -> Tuple[List[Species], List[Reaction]]:
        # update with global mechanisms

        global_mech_species = []
        global_mech_reactions = []
        if self.global_mechanisms:
            for mech in self.global_mechanisms:
                # Update Global Mechanisms
                global_mech_species += self.global_mechanisms[mech].update_species_global(
                    species, self, compartment)
                global_mech_reactions += self.global_mechanisms[mech].update_reactions_global(
                    species, self, compartment)

        self.add_species_to_crn(global_mech_species,
                                component=None, compartment=compartment)
        self.crn.add_reactions(global_mech_reactions, compartment=compartment)
        return global_mech_species, global_mech_reactions

    def component_enumeration(self, comps_to_enumerate=None, recursion_depth=10) -> List[Component]:
        # Components that produce components through Component Enumeration

        all_components = []
        new_components = []
        if (comps_to_enumerate is None):
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
                        previously_enumerated=all_components+new_components)
                    new_components += enumerated

                all_components += comps_to_enumerate
                comps_to_enumerate = list(new_components)
                new_components = []

        if (len(comps_to_enumerate) > 0):
            warn("Mixture was left with unenumerated components " +
                 str(', '.join([str(c) for c in comps_to_enumerate])))
        return all_components

    def global_component_enumeration(self, comps_to_enumerate=None, recursion_depth=None) -> List[Component]:
        """components that produce other components infinitely"""
        if (recursion_depth is None):
            recursion_depth = self.global_recursion_depth

        if (comps_to_enumerate is None):
            # these go into the ComponentEnuemrators
            comps_to_enumerate = list(self.components)

        # Recursion depth
        enumerated_components = list(
            comps_to_enumerate)  # These will be returned
        for global_enumerator in self.global_component_enumerators:
            # reset the enumeration if there's any stored info
            global_enumerator.reset(enumerated_components)
        for a in range(recursion_depth):
            # these will be added to comps_to_enumerate at the end of the iteration
            new_comps_to_enumerate = []

            for global_enumerator in self.global_component_enumerators:
                enumerated = global_enumerator.enumerate_components(comps_to_enumerate,
                                                                    # this should be only NEWLY CREATED components
                                                                    previously_enumerated=enumerated_components)
                for c in enumerated:
                    # These components are passed into the enumerator next recursion
                    if c not in new_comps_to_enumerate and c not in comps_to_enumerate:
                        new_comps_to_enumerate.append(c)

                    # These components are returend
                    if c not in enumerated_components:
                        enumerated_components.append(c)

            # Update comps_to_enumerate
            comps_to_enumerate += new_comps_to_enumerate

        return enumerated_components

    def compile_crn(self,
                    recursion_depth: int = None,
                    initial_concentration_dict: dict = None,
                    return_enumerated_components: bool = False,
                    initial_concentrations_at_end: bool = False,
                    copy_objects: bool = True,
                    add_reaction_species: bool = True,
                    compartment: Compartment = None
                    ) -> ChemicalReactionNetwork:
        """
        Creates a chemical reaction network from the species
        and reactions associated with a mixture object (with compartment).
        :param initial_concentration_dict: a dictionary to overwrite initial 
               concentrations at the end of compile time
        :param recursion_depth: how deep to run the Local and 
                                Global Component Enumeration
        :param return_enumerated_components: returns a list of all enumerated
                                             components along with the CRN
        :param initial_concentrations_at_end: if True does not look in Components
                                              for Species' initial concentrations
                                              and only checks the Mixture 
                                              database at the end.
        :param copy_objects: Species and Reactions will be copied when placed
                             into the CRN. Protects CRN validity at the expense
                             of compilation speed.
        :param add_reaction_species: Species inside reactions will be 
                                     added to the CRN. Ensures no missing species
                                     at the expense of compilation speed.
        :return: ChemicalReactionNetwork
        """
        resetwarnings()

        if compartment is None and hasattr(self, "compartment"):
            compartment = self.compartment

        self.crn = ChemicalReactionNetwork([], [])

        # helper to flatten & drop duplicates via ==, retagging default→compartment
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
                if compartment and getattr(s, "compartment", None) and s.compartment.name == "default":
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
                compartment=compartment
            )

        # enumerate components
        if recursion_depth is None:
            recursion_depth = self.global_recursion_depth
        global_comps = self.global_component_enumeration(
            recursion_depth=recursion_depth)
        enumerated = self.component_enumeration(
            global_comps,
            recursion_depth=self.local_recursion_depth
        )

        for comp in enumerated:
            comp.set_mixture(self)
            if compartment is not None:
                comp.compartment = compartment

        # component species
        for comp in enumerated:
            sp_list = comp.update_species()
            # Use _filter_new_by_eq to prevent duplicates, just like for reaction species
            existing = list(self.crn._species_set)
            fresh = _filter_new_by_eq(sp_list, existing)
            if fresh:
                self.add_species_to_crn(
                    fresh,
                    component=comp,
                    no_initial_concentrations=initial_concentrations_at_end,
                    copy_species=copy_objects,
                    compartment=compartment
                )

        # component reactions & their species
        for comp in enumerated:
            rxns = comp.update_reactions()
            self.crn.add_reactions(
                rxns,
                copy_reactions=copy_objects,
                add_species=False,
                compartment=compartment
            )

            # retag species with compartments inside these new reactions
            new_rxns = self.crn._reactions[-len(rxns):]
            for r in new_rxns:
                for w in (r.inputs + r.outputs):
                    if compartment and w.species.compartment.name == "default":
                        w.species.compartment = compartment

            # collect & add only new species
            if add_reaction_species:
                participants = [
                    w.species for r in new_rxns for w in (r.inputs + r.outputs)]
                existing = list(self.crn._species_set)
                fresh = _filter_new_by_eq(participants, existing)
                if fresh:
                    self.add_species_to_crn(
                        fresh,
                        component=comp,
                        no_initial_concentrations=initial_concentrations_at_end,
                        copy_species=copy_objects,
                        compartment=compartment
                    )

        # global mechanisms
        global_sp, global_rxns = [], []
        for mech in self.global_mechanisms.values():
            global_sp.extend(mech.update_species_global(list(self.crn._species_set),
                                                        self))
            global_rxns.extend(mech.update_reactions_global(list(self.crn._species_set),
                                                            self))

        if global_sp:
            existing = list(self.crn._species_set)
            fresh = _filter_new_by_eq(global_sp, existing)
            if fresh:
                self.add_species_to_crn(
                    fresh,
                    component=None,
                    no_initial_concentrations=initial_concentrations_at_end,
                    copy_species=copy_objects,
                    compartment=compartment
                )

        if global_rxns:
            self.crn.add_reactions(
                global_rxns,
                copy_reactions=copy_objects,
                add_species=False,
                compartment=compartment
            )
            new_grx = self.crn._reactions[-len(global_rxns):]
            for r in new_grx:
                for w in (r.inputs + r.outputs):
                    if compartment and w.species.compartment.name == "default":
                        w.species.compartment = compartment

            if add_reaction_species:
                participants = [
                    w.species for r in new_grx for w in (r.inputs + r.outputs)]
                existing = list(self.crn._species_set)
                fresh = _filter_new_by_eq(participants, existing)
                if fresh:
                    self.add_species_to_crn(
                        fresh,
                        component=None,
                        no_initial_concentrations=initial_concentrations_at_end,
                        copy_species=copy_objects,
                        compartment=compartment
                    )

        # initial‐conc handling
        if initial_concentrations_at_end:
            self.crn.initial_concentration_dict = self.get_initial_concentration(
                list(self.crn._species_set), component=None
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
        txt = str(self)+"\n"
        if self.components:
            txt += "Components = ["
            for comp in self.components:
                txt += "\n\t"+str(comp)
        if self.mechanisms:
            txt += " ]\nMechanisms = {"
            for mech in self.mechanisms:
                txt += "\n\t"+mech+":"+self.mechanisms[mech].name
        if self.global_mechanisms:
            txt += " }\nGlobal Mechanisms = {"
            for mech in self.global_mechanisms:
                txt += "\n\t"+mech+":"+self.global_mechanisms[mech].name
        txt += " }"
        return txt
