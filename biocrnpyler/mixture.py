
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from typing import List, Union
from warnings import resetwarnings, warn

from .chemical_reaction_network import ChemicalReactionNetwork
from .component import Component
from .global_mechanism import GlobalMechanism
from .mechanism import Mechanism
from .parameter import ParameterDatabase
from .reaction import Reaction
from .species import Species


class Mixture(object):
    def __init__(self, name="", mechanisms=None, components=None, parameters=None, parameter_file=None,
                 global_mechanisms=None, species=None, **kwargs):
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
        else:
            self.add_mechanisms(global_mechanisms)

        # process the species
        self.add_species(species)

        # Create a paraemter database
        self.parameter_database = ParameterDatabase(parameter_file = parameter_file, parameter_dictionary = parameters, **kwargs)
        
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

            assert all(isinstance(x, Species) for x in species_list), 'only Species type is accepted!'

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
            raise ValueError("Invalid Species: string, chemical_reaction_network.Species or Component with implemented .get_species() required as input.")

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
            assert isinstance(component, Component), "the object: %s passed into mixture as component must be of the class Component" % str(component)

            # Check if component is already in self._components
            for comp in self._components:
                if type(comp) == type(component) and comp.name == component.name:
                    raise ValueError(f"{comp} of the same type and name already in Mixture!")
            else:
                # Components are copied before being added to Mixtures
                component_copy = copy.deepcopy(component)
                component_copy.set_mixture(self)
                self.components.append(component_copy)

    def get_mechanism(self, mechanism_type):
        """Searches the Mixture for a Mechanism of the correct type.

        If no Mechanism is found, None is returned.
        """
        if not isinstance(mechanism_type, str):
            raise TypeError(f"mechanism_type must be a string. Recievied {mechanism_type}.")

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
            raise ValueError(f"add_components expected a list of Components. Received {components}")

    def get_component(self, component=None, name=None, index=None):
        """Function to get components from Mixture._components.

        One of the 3 keywords must not be None.

        :param component: an instance of a component. Searches Mixture._components for a Component with the same type and name.
        :param name: str. Searches Mixture._components for a Component with the same name
        :param index: int. returns Mixture._components[index]
        :return: if nothing is found, returns None.
        """
        if [component, name, index].count(None) != 2:
            raise ValueError(f"get_component requires a single keyword. Received component={component}, name={name}, index={index}.")
        if not (isinstance(component, Component) or component is None):
            raise ValueError(f"component must be of type Component. Received {component}.")
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
            warn("get_component found multiple matching components. A list has been returned.")
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
            raise TypeError(f"mechanism must be a Mechanism. Received {mechanism}.")

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(f"mechanism keys must be strings. Received {mech_type}")

        if isinstance(mechanism, GlobalMechanism):
            self.add_global_mechanism(mechanism, mech_type, overwrite)
        elif isinstance(mechanism, Mechanism):
            if mech_type in self._mechanisms and not overwrite:
                raise ValueError(f"mech_type {mech_type} already in Mixture {self}. To overwrite, use keyword overwrite = True.")
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
            self.add_mechanism(mechanisms, overwrite = overwrite)
        elif isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_mechanism(mechanisms[mech_type], mech_type, overwrite = overwrite)
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_mechanism(mech, overwrite = overwrite)
        else:
            raise ValueError(f"add_mechanisms expected a list of Mechanisms. Recieved {mechanisms}")

                
    @property
    def global_mechanisms(self):
        """global_mechanisms stores global Mechanisms in the Mixture."""
        return self._global_mechanisms

    @global_mechanisms.setter
    def global_mechanisms(self, mechanisms):
        self._global_mechanisms = {}
        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_global_mechanism(mechanisms[mech_type], mech_type, overwrite = True)
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_global_mechanism(mech, overwrite = True)

    def add_global_mechanism(self, mechanism, mech_type = None, overwrite = False):
        """adds a mechanism of type mech_type to the Mixture global_mechanism dictonary.

        keywordS:
          mechanism: a Mechanism instance
          mech_type: the type of mechanism. defaults to mechanism.mech_type if None
          overwrite: whether to overwrite existing mechanisms of the same type (default False)
        """
        if not hasattr(self, "_global_mechanisms"):
            self._global_mechanisms = {}

        if not isinstance(mechanism, GlobalMechanism):
            raise TypeError(f"mechanism must be a GlobalMechanism. Recieved {mechanism}.")

        if mech_type is None:
            mech_type = mechanism.mechanism_type
        if not isinstance(mech_type, str):
            raise TypeError(f"mechanism keys must be strings. Recieved {mech_type}")

        if mech_type in self._mechanisms and not overwrite:
            raise ValueError(f"mech_type {mech_type} already in Mixture {self}. To overwrite, use keyword overwrite = True.")
        else:
            self._global_mechanisms[mech_type] = copy.deepcopy(mechanism)

    def update_parameters(self, parameter_file = None, parameters = None, overwrite_parameters = True):
        if parameter_file is not None:
            self.parameter_database.load_parameters_from_file(parameter_file, overwrite_parameters = overwrite_parameters)

        if parameters is not None:
            self.parameter_database.load_parameters_from_dictionary(parameters, overwrite_parameters = overwrite_parameters)
    
    def get_parameter(self, mechanism, part_id, param_name):
        param = self.parameter_database.find_parameter(mechanism, part_id, param_name)

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

        init_conc_dict = {}
        for s in S:
            if not isinstance(s, Species):
                raise ValueError(f"{s} is not a Species! Can only find initial concentration of a Species.")

            init_conc = None
            #1 Check the component
            if component is not None:
                init_conc = component.get_parameter(param_name = str(s), part_id = self.name, mechanism = "initial concentration", check_mixture = False, return_none = True)

                if init_conc is None and component.get_species() == s:
                    init_conc = component.get_parameter(param_name = component.name, part_id = self.name, mechanism = "initial concentration", check_mixture = False, return_none = True)
                
            #2 Check self
            if init_conc is None:
                init_conc = self.get_parameter(param_name = str(s), part_id = self.name, mechanism = "initial concentration")

                if init_conc is None and component is not None and component.get_species() == s:
                    init_conc = self.get_parameter(param_name = component.name, part_id = self.name, mechanism = "initial concentration")

            if init_conc is None:
                init_conc = 0

            init_conc_dict[s] = init_conc

        return init_conc_dict

    def add_species_to_crn(self, new_species, component):

        if self.crn is None:
            self.crn = ChemicalReactionNetwork(species = [], reactions = [])

        if isinstance(new_species, Species):
            new_species = [new_species]

        for s in new_species:
            if isinstance(s, Species) or(isinstance(s, list) and(all(isinstance(ss, Species) for ss in s) or len(s) == 0)):
                init_conc_dict = self.get_initial_concentration(s, component)
                self.crn.add_species(s)
                self.crn.initial_concentration_dict = init_conc_dict
            elif s is not None:
                raise ValueError(f"Invalid Species Returned in {component}.update_species(): {s}.")

    def apply_global_mechanisms(self, species) -> (List[Species], List[Reaction]):
        # update with global mechanisms

        global_mech_species = []
        global_mech_reactions = []
        if self.global_mechanisms:
            for mech in self.global_mechanisms:
                # Update Global Mechanisms
                global_mech_species += self.global_mechanisms[mech].update_species_global(species, self)
                global_mech_reactions += self.global_mechanisms[mech].update_reactions_global(species, self)

        self.add_species_to_crn(global_mech_species, component = None)
        self.crn.add_reactions(global_mech_reactions)

    def compile_crn(self, initial_concentration_dict = None) -> ChemicalReactionNetwork:
        """Creates a chemical reaction network from the species and reactions associated with a mixture object.
        :param initial_concentration_dict: a dictionary to overwride initial concentrations at the end of compile time
        :return: ChemicalReactionNetwork
        """
        resetwarnings()#Reset warnings - better to toggle them off manually.

        #reset the Components' mixture to self - in case they have been added to other Mixtures
        for c in self.components:
            c.set_mixture(self)

        #Create a CRN to filter out duplicate species
        self.crn = ChemicalReactionNetwork([], [])

        #add the extra species to the CRN
        self.add_species_to_crn(self.added_species, component = None)

        #Append Species from each Component
        for component in self.components:
            self.add_species_to_crn(component.update_species(), component)

        #Append Reactions from each Component
        for component in self.components:
            self.crn.add_reactions(component.update_reactions())

        #global mechanisms are applied last and only to all the species
        #the reactions and species are added to the CRN
        self.apply_global_mechanisms(self.crn.species)

        #Manually change/override initial conditions at compile time
        if initial_concentration_dict is not None:
            self.crn.initial_concentration_dict = initial_concentration_dict

        return self.crn

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self)+"\n"
        if self.components:
            txt += "Components = ["
            for comp in self.components:
                txt+="\n\t"+str(comp)
        if self.mechanisms:
            txt+=" ]\nMechanisms = {"
            for mech in self.mechanisms:
                txt+="\n\t"+mech+":"+self.mechanisms[mech].name
        if self.global_mechanisms:
            txt+=" }\nGlobal Mechanisms = {"
            for mech in self.global_mechanisms:
                txt+="\n\t"+mech+":"+self.global_mechanisms[mech].name
        txt+=" }"
        return txt
