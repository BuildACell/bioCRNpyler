#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from warnings import warn as pywarn
from .chemical_reaction_network import Species, ComplexSpecies, Reaction
from .parameter import Parameter
from typing import List, Union


def warn(txt):
    pywarn(txt)


# Component class for core components
class Component(object):

    def __init__(self, name: Union[str, Species],
                 mechanisms={},  # custom mechanisms
                 parameters=None,  # parameter configuration
                 parameter_file=None, #custom parameter file
                 mixture=None,
                 attributes=[],
                 initial_conc=0,
                 parameter_warnings = True,
                 **keywords  # parameter keywords
                 ):

        if isinstance(name, Species):
            self.name = name.name
        elif isinstance(name, str):
            self.name = name
        else:
            raise ValueError("name must be a Species or string")

        # Toggles whether warnings will be sent when parameters aren't found by
        # the default name.
        self.set_parameter_warnings(parameter_warnings)
        self._initial_conc = initial_conc

        # Check to see if a subclass constructor has overwritten default
        # mechanisms.
        # Attributes can be used to store key words like protein deg-tags for
        # components that mimic CRN species.
        self.attributes = []
        self.set_attributes(attributes)

        # Check to see if a subclass constructor has overwritten default
        # mechanisms.
        if not hasattr(self, 'default_mechanisms'):
            self.default_mechanisms = {}

        self.custom_mechanisms = {}
        self.mechanisms = {}
        if mixture is not None:
            mixture_mechanisms = mixture.mechanisms
        else:
            mixture_mechanisms = {}
        self.update_mechanisms(mechanisms=mechanisms,
                               mixture_mechanisms=mixture_mechanisms)

        self.custom_parameters = {}
        self.parameters = {}
        if mixture is not None:
            mixture_parameters = mixture.parameters
        else:
            mixture_parameters = {}

        parameters = Parameter.create_parameter_dictionary(parameters, parameter_file)
        self.update_parameters(mixture_parameters=mixture_parameters, parameters=parameters)

    @property
    def initial_concentration(self) -> float:
        return self._initial_conc

    @initial_concentration.setter
    def initial_concentration(self, initial_conc: float):
        if initial_conc is None:
            self._initial_conc = initial_conc
        elif initial_conc < 0.0:
            raise ValueError("Initial concentration must be non-negative, "f"this was given: {initial_conc}")
        else:
            self._initial_conc = initial_conc

    # TODO implement as an abstractmethod
    def get_species(self) -> None:
        """
        the child class should implement this method
        :return: None
        """
        warn(f"get_species is not defined for component {self.name}, None returned.")
        return None

    #If allows species to be set from strings, species, or Components
    def set_species(self, species: Union[Species, str], material_type = None, attributes = None):
        if isinstance(species, Species):
                return species
        elif isinstance(species, str):
            return Species(name = species, material_type = material_type, attributes = attributes)
        elif isinstance(species, Component) and species.get_species() is not None:
            return species.get_species()
        else:
            raise ValueError("Invalid Species: string, chemical_reaction_network.Species or Component with implemented .get_species() required as input.")

    def __hash__(self):
        return str.__hash__(repr(self.get_species()))

    def set_attributes(self, attributes: List[str]):
        if attributes is not None:
            for attribute in attributes:
                self.add_attribute(attribute)

    def add_attribute(self, attribute: str):
        assert isinstance(attribute, str) and attribute is not None, "Attribute: %s must be a str" % attribute

        self.attributes.append(attribute)
        if hasattr(self, 'species') and self.species is not None:
            self.species.add_attribute(attribute)
        else:
            raise Warning(f"Component {self.name} has no internal species and therefore no attributes")

    #This function gives Components their own parameter dictionary. By default, components get all the parameters from Mixture
    #Parameters passed in from mixture are superceded by parameters passed in the parameters keyword
    #parameters already saved as custom parameters also supersede parameters from the Mixture
    #In other words, the component remembers custom changes to parameters and mixture parameters will never overwrite those changes
    #Mixture parameters are only ever used by default when no other component-level parameter has ever been given with the same key.
    #the overwrite_custom_parameters keyword lets this function overwrite the existing custom Component-level parameters,
    def update_parameters(self, mixture_parameters = {}, parameters = {},
                          overwrite_custom_parameters = True):


        for p in parameters:
            if overwrite_custom_parameters or p not in self.custom_parameters:
                self.parameters[p] = parameters[p]
                if p in self.custom_parameters:
                    self.custom_parameters[p] = parameters[p]


        for p in mixture_parameters:
            if p not in self.parameters:
                self.parameters[p] = mixture_parameters[p]

    def update_mechanisms(self, mixture_mechanisms={}, mechanisms={},
                          overwrite_custom_mechanisms=True):

        for mech_type in mixture_mechanisms:
            self.mechanisms[mech_type] = mixture_mechanisms[mech_type]

        # The mechanisms used during compilation are stored as their own
        # dictionary
        for mech_type in self.default_mechanisms:
            self.mechanisms[mech_type] = self.default_mechanisms[mech_type]

        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                if overwrite_custom_mechanisms \
                   or mech_type not in self.custom_mechanisms:
                    self.mechanisms[mech_type] = mechanisms[mech_type]
                    self.custom_mechanisms[mech_type] = mechanisms[mech_type]
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                if overwrite_custom_mechanisms \
                   or mech not in self.custom_mechanisms:
                    self.mechanisms[mech.mechanism_type] = mech
                    self.custom_mechanisms[mech.mechanism_type] = mech
        else:
            raise ValueError("Mechanisms must be passed as a list of "
                             "instantiated objects or a dictionary "
                             "{type:mechanism}")

    #Set get_parameter property
    def set_parameter_warnings(self, parameter_warnings):
        self.parameter_warnings = parameter_warnings

    # Get Parameter Hierarchy:
    def get_parameter(self, param_name: str, part_id=None, mechanism=None):
        return_val = None
        warning_txt = None

        # Ideally parameters can be found
        # (mechanism.name/mechanism_type, part_id, param_name) --> val
        if part_id is not None and mechanism is not None:
            if mechanism is not None \
               and (mechanism.name, part_id, param_name) in self.parameters:
                return_val = self.parameters[(mechanism.name,
                                              part_id, param_name)]
            elif mechanism is not None \
                 and (mechanism.mechanism_type, part_id, param_name) \
                     in self.parameters:
                return_val = self.parameters[(mechanism.mechanism_type, part_id,
                                              param_name)]

        # Next try (part_id, param_name) --> val
        if part_id is not None and return_val is None:
            if (part_id, param_name) in self.parameters:
                return_val = self.parameters[(part_id, param_name)]

                if mechanism is not None:
                    warning_txt = ("No Parameter found with "
                        f"param_name={param_name} and part_id={part_id} and "
                        f"mechanism={repr(mechanism)}. Parameter found under "
                        f"the key (part_id, param_name)=({part_id}, "
                        f"{param_name}).")
        # Next try (Mechanism.name/mechanism_type, param_name) --> val
        if mechanism is not None and return_val is None:
            if (mechanism.name, param_name) in self.parameters:
                return_val = self.parameters[((mechanism.name, param_name))]
                if part_id is not None:
                    warning_txt = ("No Parameter found with "
                        f"param_name={param_name} and part_id={part_id} and "
                        f"mechanism={repr(mechanism)}. Parameter found under "
                        f"the key (mechanism.name, "
                        f"param_name)=({mechanism.name}, {param_name})")
            elif (mechanism.mechanism_type, param_name) in self.parameters:
                return_val = self.parameters[(mechanism.mechanism_type,
                                              param_name)]
                if part_id is not None:
                    warning_txt = ("No Parameter found with "
                        f"param_name={param_name} and part_id={part_id} and "
                        f"mechanism={repr(mechanism)}. Parameter found under "
                        f"the key (mechanism.name, "
                        f"param_name)=({mechanism.name}, {param_name})")
        # Finally try (param_name) --> return val
        if param_name in self.parameters and return_val is None:
            return_val = self.parameters[param_name]
            if mechanism is not None or part_id is not None:
                warning_txt = (f"No parameter found with "
                               f"param_name={param_name} and part_id={part_id} "
                               f"and mechanism={repr(mechanism)}. Parameter "
                               f"found under the key param_name={param_name}")
        if return_val is None:
            raise ValueError("No parameters can be found that match the "
                             "(mechanism, part_id, "
                            f"param_name)=({repr(mechanism)}, {part_id}, "
                            f"{param_name})")

        else:
            if (warning_txt is not None) and self.parameter_warnings:
                warn(warning_txt)
            return return_val

    # TODO implement abstractmethod
    def update_species(self) -> List:
        """
        the child class should implement this method
        :return: empty list
        """
        species = []
        warn("Unsubclassed update_species called for " + repr(self))
        return species

    # TODO implement abstractmethod
    def update_reactions(self) -> List:
        """
        the child class should implement this method
        :return: empty list
        """
        reactions = []
        warn("Unsubclassed update_reactions called for " + repr(self))
        return reactions

    def __repr__(self):
        return type(self).__name__ + ": " + self.name

