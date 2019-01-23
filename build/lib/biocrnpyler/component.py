# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn as pywarn


def warn(txt):
    pywarn(txt)


# import chemical_reaction_network as crn
from .chemical_reaction_network import Specie


# Component class for core components
class Component(object):

    def __init__(self, name,
                 mechanisms={},  # custom mechanisms
                 parameters={},  # parameter configuration
                 mixture=None,
                 attributes=[],
                 initial_conc=0,
                 parameter_warnings = True,
                 **keywords  # parameter keywords
                 ):

        self.name = name

        self.parameter_warnings = parameter_warnings#Toggles whether warnings will be sent when parameters aren't found by the default name
        self._initial_conc = initial_conc

        #Check to see if a subclass constructor has overwritten default mechanisms
        # Attributes can be used to store key words like protein deg-tags for components that mimic CRN species
        self.set_attributes(attributes)

        # Check to see if a subclass constructor has overwritten default mechanisms
        if not hasattr(self, 'default_mechanisms'):
            self.default_mechanisms = {}

        self.custom_mechanisms = {}
        self.mechanisms = {}

        if mixture is not None:
            mixture_mechanisms = mixture.mechanisms
        else:
            mixture_mechanisms = {}
        self.update_mechanisms(mechanisms=mechanisms, mixture_mechanisms=mixture_mechanisms)

        self.custom_parameters = {}
        self.parameters = {}
        if mixture is not None:
            mixture_parameters = mixture.parameters
        else:
            mixture_parameters = {}
        self.update_parameters(mixture_parameters=mixture_parameters, parameters=parameters)

    #@property
    def initial_concentration(self):
        return self._initial_conc

    #@initial_concentration.setter
    def initial_concentration(self, initial_conc):
        if initial_conc is None:
            self._initial_conc = initial_conc
        elif initial_conc < 0:
            raise ValueError('Initial concentration must be non-negative, this was given: %d' % initial_conc)
        else:
            self._initial_conc = initial_conc

    def get_specie(self):
        warn("get_specie() not defined for component " + self.name + " None returned.")
        return None

    def __hash__(self):
        return str.__hash__(repr(self.get_specie()))

    def set_attributes(self, attributes):
        self.attributes = []
        for attribute in attributes:
            if isinstance(attribute, str):
                self.attributes.append(attribute)
            elif attribute is not None:
                raise RuntimeError("Invalid Attribute: "+repr(attribute)+" attributes must be strings")

    def add_attribute(self, attribute):
        if isinstance(attribute, str):
            self.attributes.append(attribute)
            self.specie.add_attribute(attribute)
        else:
            raise ValueError("Attribute must be a str")

    def update_parameters(self, mixture_parameters = {}, parameters = {}, overwrite_custom_parameters = True):
        for p in parameters:
            if overwrite_custom_parameters or p not in self.custom_parameters:
                self.parameters[p] = parameters[p]
                if p in self.custom_parameters:
                    self.custom_parameters[p] = parameters[p]

        for p in mixture_parameters:
            if p not in self.parameters:
                self.parameters[p] = mixture_parameters[p]

    def update_mechanisms(self, mixture_mechanisms={}, mechanisms={}, overwrite_custom_mechanisms=True):

        for mech_type in mixture_mechanisms:
            self.mechanisms[mech_type] = mixture_mechanisms[mech_type]

        # The mechanisms used during compilation are stored as their own dictionary
        for mech_type in self.default_mechanisms:
            self.mechanisms[mech_type] = self.default_mechanisms[mech_type]

        if isinstance(mechanisms, dict):
            for mech_type in mechanisms:
                if overwrite_custom_mechanisms or mech_type not in self.custom_mechanisms:
                    self.mechanisms[mech_type] = mechanisms[mech_type]
                    self.custom_mechanisms[mech_type] = mechanisms[mech_type]
        elif isinstance(mechanisms, list):
            for mech in mechanisms:
                if overwrite_custom_mechanisms or mech not in self.custom_mechanisms:
                    self.mechanisms[mech.type] = mech
                    self.custom_mechanisms[mech.type] = mech
        else:
            raise ValueError(
                "Mechanisms must be passed as a list of instantiated objects or a dictionary {type:mechanism}")

    # Get Parameter Hierarchy:
    def get_parameter(self, param_name, part_id=None, mechanism=None):
        return_val = None
        warning_txt = None

        # Ideally parameters can be found (mechanism.name/type, part_id, param_name) --> val
        if part_id is not None and mechanism is not None:
            if mechanism is not None and (mechanism.name, part_id, param_name) in self.parameters:
                return_val = self.parameters[(mechanism.name, part_id, param_name)]
            elif mechanism is not None and (mechanism.type, part_id, param_name) in self.parameters:
                return_val = self.parameters[(mechanism.type, part_id, param_name)]

        # Next try (part_id, param_name) --> val
        if part_id is not None and return_val is None:
            if (part_id, param_name) in self.parameters:
                return_val = self.parameters[(part_id, param_name)]

                if mechanism is not None:
                    warning_txt = "No Parameter found with param_name=" + param_name + " and part_id=" + part_id + " and mechanism=" + repr(
                        mechanism) + ". Parameter found under the key (part_id, param_name)=(" + part_id + ", " + param_name + ")"
        # Next try (Mechanism.name/type, param_name) --> val
        if mechanism is not None and return_val is None:
            if (mechanism.name, param_name) in self.parameters:
                return_val = self.parameters[((mechanism.name, param_name))]
                if part_id is not None:
                    warning_txt = "No Parameter found with param_name=" + param_name + " and part_id=" + part_id + " and mechanism=" + repr(
                        mechanism) + ". Parameter found under the key (mechanism.name, Component.name, param_name)=(" \
                                  + mechanism.name + ", " + self.name + ", " + param_name + ")"
            elif (mechanism.type, param_name) in self.parameters:
                return_val = self.parameters[(mechanism.type, param_name)]
                if part_id is not None:
                    warning_txt = "No Parameter found with param_name=" + param_name + " and part_id=" + part_id + " and mechanism=" + repr(
                        mechanism) + ". Parameter found under the key (mechanism.name, param_name)=(" \
                                  + mechanism.name + ", " + param_name + ")"
        # Finally try (param_name) --> return val
        if param_name in self.parameters and return_val is None:
            return_val = self.parameters[param_name]
            if mechanism is not None or part_id is not None:
                warning_txt = "No Parameter found with param_name=" + param_name + " and part_id=" + str(
                    part_id) + " and mechanism=" + repr(
                    mechanism) + ". Parameter found under the key param_name=" + param_name
        if return_val is None:
            raise ValueError("No Parameters can be found that match the (mechanism, param_id, param_name)=( " + repr(
                mechanism) + ', ' + part_id + ", " + param_name + ")")
        else:
            if warning_txt is not None and self.parameter_warnings:
                warn(warning_txt)
            return return_val

    def update_species(self):
        species = []
        warn("Unsubclassed update_species called for " + repr(self))
        return species

    def update_reactions(self):
        reactions = []
        warn("Unsubclassed update_reactions called for " + repr(self))
        return reactions

    def __repr__(self):
        return type(self).__name__ + ": " + self.name


# These subclasses of Component represent different kinds of biomolecules.


class DNA(Component):
    """DNA class

    The DNA class is used to represent a DNA sequence that has a given
    length.  Its main purpose is as the parent object for DNA
    fragments and DNA assemblies.

    Note: for initialization of members of this class, the arguments
    should be as follows:

      DNA(name, length, [mechanisms], [config_file], [prefix])

        DNAtype(name, required_arguments, [length], [mechanisms],
                [config_file], [prefix], [optional_arguments])

          DNAelement(name, required_arguments, [length], [mechanisms],
                     [config_file], [optional_arguments])


    Data attributes
    ---------------
    name        Name of the sequence (str)
    length      Length of the sequence (int)
    assy        DNA assembly that we are part of
    mechanisms  Local mechanisms for this component (overrides defaults)
    parameters  Parameter dictionary for the DNA element

    """

    def __init__(
            self, name, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            initial_conc=None,
            parameter_warnings = True,
            **keywords
    ):
        self.specie = Specie(self.name, type="dna", attributes=list(attributes))
        self._length = length
        Component.__init__(self=self, name=name, mechanisms=mechanisms, parameters=parameters,
                           attributes=attributes, initial_conc=initial_conc, parameter_warnings = parameter_warnings, **keywords)

    def get_specie(self):
        return self.specie

    def update_species(self):
        species = [self.get_specie()]
        return species

    def update_reactions(self):
        return []

    @property
    def length(self):
        return  self._length

    @length.setter
    def length(self, dna_length):
        if dna_length >= 0:
            self._length = dna_length
        else:
            raise ValueError("Length cannot be negative!")


class RNA(Component):
    def __init__(
            self, name, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            initial_conc=None,
            **keywords
    ):
        self.length = length
        self.specie = Specie(self.name, type="rna", attributes=list(attributes))
        Component.__init__(self=self, name=name, mechanisms=mechanisms, parameters=parameters,
                           attributes=attributes, initial_conc=initial_conc, **keywords)

    def get_specie(self):
        return self.specie

    def update_species(self):
        species = [self.get_specie()]
        return species

    def update_reactions(self):
        return []


class Protein(Component):
    def __init__(
            self, name, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            degredation_tag=None,
            initial_conc=None,
            **keywords
    ):
        self.length = length
        self.degredation_tag = degredation_tag
        if degredation_tag not in attributes:
            attributes.append(degredation_tag)
        self.specie = Specie(name, type="protein", attributes=list(attributes))

        Component.__init__(self=self, name=name, mechanisms=mechanisms, parameters=parameters,
                           attributes=attributes, initial_conc=initial_conc, **keywords)

    def get_specie(self):
        return self.specie

    def update_species(self):
        species = [self.get_specie()]
        return species

    def update_reactions(self):
        return []

#A complex forms when two or more species bind together
#Complexes inherit the attributes of their species
class Complex(Component):
    def __init__(
            self, species,  # positional arguments
            name = None, #Override the default naming convention for a complex
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters,
            attributes=[],
            initial_conc=None,
            **keywords
    ):
        self.specie = Specie(self.name, type="complex", attributes=list(attributes))
        Component.__init__(self=self, name=name, mechanisms=mechanisms, parameters=parameters, attributes=attributes,
                           initial_conc=initial_conc, **keywords)

    def get_specie(self):
        return self.specie

    def update_species(self):
        species = [self.get_specie()]
        return species

    def update_reactions(self):
        return []
