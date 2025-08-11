#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
from numbers import Real
from typing import List, Union
from warnings import warn

from .mechanism import Mechanism
from .compartment import Compartment
from ..mechanisms.global_mechanisms import GlobalMechanism
from .parameter import Parameter, ParameterDatabase, ParameterKey
from .species import Species


class Component:
    """Component class for core components.

    These subclasses of Component represent different kinds of biomolecules.

    This class must be Subclassed to provide functionality with the functions get_species and get_reactions overwritten.
    """

    def __init__(
        self,
        name: Union[str, Species],
        mechanisms=None,  # custom mechanisms
        parameters=None,  # parameter configuration
        parameter_file=None,  # custom parameter file
        mixture=None,
        compartment=None,
        attributes=None,
        # This is added as a parameter ("initial concentration", None, self.name):initial_concentration
        initial_concentration=None,
        initial_condition_dictionary=None,
        **keywords,  # parameter keywords
    ):
        """Initializes a Component object.

        :param name:
        :param mechanisms:
        :param parameters:
        :param parameter_file:
        :param mixture:
        :param compartment:
        :param attributes:
        :param initial_concentration:
        :param initial_condition_dictionary:
        :param keywords:
        """
        if mechanisms is None:
            self.mechanisms = {}
        else:
            self.mechanisms = mechanisms
        if isinstance(name, Species):
            self.name = name.name
        elif isinstance(name, str):
            self.name = name
        else:
            raise ValueError("name must be a Species or string")

        # Attributes can be used to store key words like protein deg-tags for
        # components that mimic CRN species.
        self.attributes = []
        self.set_attributes(attributes)

        if mixture is not None:
            self.set_mixture(mixture)
        else:
            self.set_mixture(None)

        if compartment is not None:
            if not isinstance(compartment, Compartment):
                raise TypeError("The compartment must be a Compartment object")
            self.compartment = compartment
        else:
            # Don't create a default compartment here to avoid circular imports
            self._compartment = None

        self.parameter_database = ParameterDatabase(
            parameter_file=parameter_file, parameter_dictionary=parameters
        )
        self.initial_concentration = initial_concentration

        # Components can also store initial conditions, just like Mixtures
        if initial_condition_dictionary is None:
            self.initial_condition_dictionary = {}
        else:
            self.initial_condition_dictionary = dict(initial_condition_dictionary)

    @property
    def initial_concentration(self):
        return self._initial_concentration

    @initial_concentration.setter
    def initial_concentration(self, initial_concentration):
        if initial_concentration is not None:
            if initial_concentration < 0:
                raise ValueError(
                    f"Initial concentration must be non-negative, this was given: {initial_concentration}"
                )
            param_key = ParameterKey(
                mechanism="initial concentration", part_id=None, name=self.name
            )
            self.parameter_database.add_parameter(
                parameter_name="initial concentration",
                parameter_value=initial_concentration,
                parameter_key=param_key,
                overwrite_parameters=True,
            )

        self._initial_concentration = initial_concentration

    @property
    def compartment(self):
        """The compartment of the Component.

        :return: Compartment
        """
        return self._compartment

    @compartment.setter
    def compartment(self, compartment: Compartment) -> None:
        """Set the compartment of the Component.

        :param compartment:
        :return: None
        """
        if compartment is not None and not isinstance(compartment, Compartment):
            raise TypeError("The compartment must be a Compartment object")
        if compartment is not None:
            self._compartment = compartment
        if hasattr(self, "_compartment") and self._compartment:
            if hasattr(self, "species"):
                if isinstance(self.species, list):
                    for s in self.species:
                        if s.compartment.name == "default":
                            s.compartment = compartment
                else:
                    if self.species.compartment.name == "default":
                        self.species.compartment = compartment
        if compartment is not None and not isinstance(compartment, Compartment):
            raise TypeError("The compartment must be a Compartment object")

    def set_mixture(self, mixture) -> None:
        """Set the mixture the Component is in.

        :param mixture:
        :return: None
        """
        self.mixture = mixture

    # TODO implement as an abstractmethod
    def get_species(self) -> None:
        """The subclasses should implement this method!

        :return: None
        """
        return None

    @classmethod
    def set_species(
        self,
        species: Union[Species, str],
        material_type=None,
        compartment=None,
        attributes=None,
    ) -> Species:
        """Helper function that allows species to be set from strings, species, or Components

        :param species: Species, str, Component
        :param material_type:
        :param compartment:
        :param attributes:
        :return: Species
        """
        if isinstance(species, Species):
            return species
        elif isinstance(species, str):
            return Species(
                name=species,
                material_type=material_type,
                compartment=compartment,
                attributes=attributes,
            )
        elif isinstance(species, Component) and species.get_species() is not None:
            return species.get_species()
        elif isinstance(species, list):
            return [
                self.set_species(
                    s,
                    material_type=material_type,
                    compartment=compartment,
                    attributes=attributes,
                )
                for s in species
            ]
        else:
            raise ValueError(
                f"Invalid Species: string, chemical_reaction_network.Species or Component with implemented .get_species() required as input. Recieved {species}"
            )

    def __hash__(self):
        return str.__hash__(repr(self.get_species()))

    def set_attributes(self, attributes: List[str]):
        if attributes is not None:
            for attribute in attributes:
                self.add_attribute(attribute)

    def add_attribute(self, attribute: str):
        assert isinstance(attribute, str) and attribute is not None, (
            "Attribute: %s must be a str" % attribute
        )

        self.attributes.append(attribute)
        if hasattr(self, "species") and self.species is not None:
            self.species.add_attribute(attribute)
        else:
            raise Warning(
                f"Component {self.name} has no internal species and therefore no attributes"
            )

    def update_parameters(
        self,
        parameter_file=None,
        parameters=None,
        parameter_database=None,
        overwrite_parameters=True,
    ):
        """Updates the ParameterDatabase inside a Component

        Possible inputs:
            parameter_file (string)
            parameters (dict)
            parameter_database (ParameterDatabase)
        """

        if parameter_file is not None:
            self.parameter_database.load_parameters_from_file(
                parameter_file, overwrite_parameters=overwrite_parameters
            )

        if parameters is not None:
            self.parameter_database.load_parameters_from_dictionary(
                parameters, overwrite_parameters=overwrite_parameters
            )

        if parameter_database is not None:
            self.parameter_database.load_parameters_from_database(
                parameter_database, overwrite_parameters=overwrite_parameters
            )

    def get_mechanism(self, mechanism_type, optional_mechanism=False):
        """Searches the Component for a Mechanism of the correct type.

        If the Component does not have the mechanism, searches the Components' Mixture for the Mechanism.

        :param mechanism_type:
        :param optional_mechanism: toggles whether an error is thrown if no mechanism is found
        :return:
        """
        if not isinstance(mechanism_type, str):
            raise TypeError(
                f"mechanism_type must be a string. Received {mechanism_type}."
            )

        mech = None
        if mechanism_type in self.mechanisms:
            return self.mechanisms[mechanism_type]
        elif self.mixture is not None:
            mech = self.mixture.get_mechanism(mechanism_type)
        if mech is None and not optional_mechanism:
            raise KeyError(
                f"Unable to find mechanism of type {mechanism_type} in Component {self} or Mixture {self.mixture}."
            )
        else:
            return mech

    @property
    def mechanisms(self):
        return self._mechanisms

    @mechanisms.setter
    def mechanisms(self, mechanisms):
        self._mechanisms = {}
        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_mechanism(mechanisms[mech_type], mech_type, overwrite=True)
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_mechanism(mech, overwrite=True)

    def add_mechanism(
        self,
        mechanism: Mechanism,
        mech_type=None,
        overwrite=False,
        optional_mechanism=False,
    ):
        """adds a mechanism of type mech_type to the Component Mechanism dictionary.

        :param mechanism:
        :param mech_type:
        :param overwrite: toggles whether the mechanism is added overwriting any mechanism with the same key.
        :param optional_mechanism: toggles whether an error is thrown if a Mechanism is added that conflicts with an exising Mechanism
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

        if mech_type not in self._mechanisms or overwrite:
            self._mechanisms[mech_type] = copy.deepcopy(mechanism)
        elif not optional_mechanism:
            raise ValueError(
                f"mech_type {mech_type} already in component {self}. To overwrite, use keyword overwrite = True."
            )

    def add_mechanisms(
        self,
        mechanisms: Union[Mechanism, GlobalMechanism],
        overwrite=False,
        optional_mechanism=False,
    ):
        """This function adds a list or dictionary of mechanisms to the mixture.

        :param mechanisms: Can take both GlobalMechanisms and Mechanisms
        :param overwrite: toggles whether the mechanism is added overwriting any mechanism with the same key.
        :param optional_mechanism: toggles whether an error is thrown if a Mechanism is added that conflicts with an exising Mechanism
        :return:
        """

        if isinstance(mechanisms, Mechanism):
            self.add_mechanism(
                mechanisms, overwrite=overwrite, optional_mechanism=optional_mechanism
            )
        elif isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                self.add_mechanism(
                    mechanisms[mech_type],
                    mech_type,
                    overwrite=overwrite,
                    optional_mechanism=optional_mechanism,
                )
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                self.add_mechanism(
                    mech, overwrite=overwrite, optional_mechanism=optional_mechanism
                )
        else:
            raise ValueError(
                f"add_mechanisms expected a list of Mechanisms. Received {mechanisms}"
            )

    def get_parameter(
        self,
        param_name: str,
        part_id=None,
        mechanism=None,
        return_numerical=False,
        return_none=False,
        check_mixture=True,
    ) -> Union[Parameter, Real]:
        """Get a parameter from different objects that hold parameters.

        Hierarchy:
         1. tries to find the Parameter in Component.parameter_database
         2. tries to find the parameter in Component.mixture.parameter_database

        :param param_name:
        :param part_id:
        :param mechanism:
        :param return_numerical: numerical value or the parameter object is returned
        :param return_none: returns None instead of throwing an error if a parameter isn't found
        :param check_mixture: toggle whether or not to check the Component's Mixture as well
        :return: Parameter object or a Real number
        """
        # Try the Component ParameterDatabase
        param = self.parameter_database.find_parameter(mechanism, part_id, param_name)

        # Next try the Mixture ParameterDatabase
        if param is None and self.mixture is not None and check_mixture:
            param = self.mixture.get_parameter(mechanism, part_id, param_name)

        if param is None and not return_none:
            raise ValueError(
                "No parameters can be found that match the "
                "(mechanism, part_id, "
                f"param_name)=({repr(mechanism)}, {part_id}, "
                f"{param_name})"
            )
        elif return_numerical and isinstance(param, Parameter):
            return param.value
        else:
            return param

    # TODO implement abstractmethod
    def update_species(self) -> List[Species]:
        """The subclasses should implement this method!

        :return: empty list
        """
        species = []
        warn("Unsubclassed update_species called for " + repr(self))
        return species

    # TODO implement abstractmethod
    def update_reactions(self) -> List[Species]:
        """The subclasses should implement this method!

        :return: empty list
        """
        reactions = []
        warn("Unsubclassed update_reactions called for " + repr(self))
        return reactions

    def enumerate_components(self, **keywords) -> List:
        """this is for component enumeration. Usually you will return a list of components that are
        copies of existing ones (first list) and new components (second list). For example,
        A DNA_construct makes a list of copies of its parts as the first output, and a list of RNA_constructs
        as the second output.
        An RNA_construct will make a list of copies of its parts as the first output, and a list of Protein
        components as its second output (if it makes any proteins)"""
        return []

    def __repr__(self):
        return type(self).__name__ + ": " + self.name
