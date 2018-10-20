from warnings import warn as pywarn
def warn(txt):
    pywarn(txt)
#from .sbmlutil import add_species, add_reaction
#from .parameter import eval_parameter
#from .parameter import get_parameters
import chemical_reaction_network as crn

# Component class for core components
class Component(object):

    def __init__(self, name,
                 mechanisms={},  # custom mechanisms
                 parameters={}, # parameter configuration
                 mixture = None,
                 attributes = [],
                 **keywords  # parameter keywords
                 ):

        self.name = name

        #Attributes can be used to store key words like protein deg-tags for components that mimic CRN species
        self.attributes = attributes

        #Check to see if a subclass constructor has overwritten default mechanisms
        if not hasattr(self, 'default_mechanisms'):
            self.default_mechanisms = {}

        self.custom_mechanisms = {}
        self.mechanisms = {}

        if mixture != None:
            mixture_mechanisms = mixture.mechanisms
        else:
            mixture_mechanisms = {}
        self.update_mechanisms(mechanisms = mechanisms, mixture_mechanisms = mixture_mechanisms)


        self.custom_parameters = {}
        self.parameters = {}
        if mixture != None:
            mixture_parameters = mixture.parameters
        else:
            mixture_parameters = {}
        self.update_parameters(mixture_parameters = mixture_parameters, parameters = parameters)

    def update_parameters(self, mixture_parameters = {}, parameters = {}, overwrite_custom_parameters = True):
        for p in parameters:
            if overwrite_custom_parameters or p not in self.custom_parameters:
                self.parameters[p] = parameters[p]
                if p in self.custom_parameters:
                    self.custom_parameters[p] = parameters[p]

        for p in mixture_parameters:
            if p not in self.parameters:
                self.parameters[p] = mixture_parameters[p]

    def update_mechanisms(self, mixture_mechanisms = {}, mechanisms = {}, overwrite_custom_mechanisms = True):


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
                if overwrite_custom_mechanisms or mech_type not in self.custom_mechanisms:
                    self.mechanisms[mech.type] = mech
                    self.custom_mechanisms[mech.type] = mech
        else:
            raise ValueError(
                "Mechanisms must be passed as a list of instantiated objects or a dictionary {type:mechanism}")


    #Get Parameter Hierarchy:
    def get_parameter(self, param_name, part_id = None, mechanism = None):
        if  part_id != None and mechanism != None:
            if mechanism != None and (mechanism.name, part_id, param_name) in self.parameters:
                return self.parameters[(mechanism.name, part_id, param_name)]
            elif mechanism != None and (mechanism.type, part_id, param_name) in self.parameters:
                return self.parameters[(mechanism.type, part_id, param_name)]
            else:
                warn("No Parameter found with param_name="+param_name+", part_id="+part_id+", and mechanism="+repr(mechanism)+". Defaulting to ([component], part_id, param_name) identifiers.")

        if part_id!= None:
            if (self.name, part_id, param_name) in self.parameters:
                return self.parameters[(self.name, part_id, param_name)]
            elif (type(self).__name__, part_id, param_name) in self.parameters:
                return self.parameters[(type(self).__name__, part_id, param_name)]
            elif (part_id, param_name) in self.parameters:
                return self.parameters[(part_id, param_name)]
            else:
                warn("No Parameter found with param_name="+param_name+" and part_id="+part_id+". Defaulting to (mechanism, [component], param_name)")

        if mechanism != None:
            if (mechanism.name, self.name, param_name) in self.parameters:
                return self.parameters[((mechanism.name, self.name, param_name))]
            elif (mechanism.name, param_name) in self.parameters:
                return (self.parameters[(mechanism.name, param_name)])
            elif (mechanism.type, self.name, param_name) in self.parameters:
                return self.parameters[((mechanism.name, self.name, param_name))]
            elif (mechanism.type, param_name) in self.parameters:
                return (self.parameters[(mechanism.type, param_name)])
            else:
                warn("No Parameter found with param_name="+param_name+" for mechanism "+repr(mechanism)+". Defaulting to (component, param_name) identifiers.")
            try:
                return super().get_parameter(param_name, mechanism = mechanism)
            except AttributeError:
                pass

        if (type(self).__name__, self.name, param_name) in self.parameters:
            return self.parameters[(type(self).__name__, self.name, param_name)]
        elif (self.name, param_name) in self.parameters:
            return self.parameters[(self.name, param_name)]
        elif (type(self).__name__, param_name) in self.parameters:
            return self.parameters[(type(self).__name__, param_name)]
        elif ('default', param_name) in self.parameters:
            return self.parameters[('default', param_name)]
        elif param_name in self.parameters:
            return self.parameters[param_name]
        else:
            try:
                return super().get_parameter(param_name)
            except AttributeError:
                raise AttributeError("No Valid Parameter '"+param_name+"' in "+repr(self))

    def update_species(self):
        species = []
        warn("Unsubclassed update_species called for "+repr(self))
        return species

    def update_reactions(self):
        reactions = []
        warn("Unsubclassed update_reactions called for "+repr(self))
        return reactions

    def __repr__(self):
        return type(self).__name__ + ": " + self.name

#These subclasses of Component represent different kinds of biomolecules.
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
            attributes = [],
            **keywords
    ):
        self.length = length
        Component.__init__(self = self, name = name, length = length, mechanisms=mechanisms, parameters = parameters, attributes = attributes, **keywords)

    def update_species(self):
        species = [crn.specie(self.name, type = "dna", attributes = self.attributes)]
        return species

    def update_reactions(self):
        return []

class RNA(Component):
    def __init__(
            self, name, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes = [],
            **keywords
    ):
        self.length = length
        Component.__init__(self = self, name = name, length = length, mechanisms=mechanisms, parameters = parameters, attributes = attributes, **keywords)

    def update_species(self):
        species = [crn.specie(self.name, type = "rna", attributes = self.attributes)]
        return species

    def update_reactions(self):
        return []

class Protein(Component):
    def __init__(
            self, name, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes = [],
            degredation_tag = None, **keywords
    ):
        self.length = length
        self.degredation_tag = degredation_tag
        if degredation_tag not in attributes:
            attributes.append(degredation_tag)
        Component.__init__(self = self, name = name, length = length, mechanisms=mechanisms, parameters = parameters, attributes = attributes, **keywords)

    def update_species(self):
        species = [crn.specie(self.name, type="protein", attributes=self.attributes)]
        return species

    def update_reactions(self):
        return []


class Complex(Component):
    def __init__(
            self, name,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters,
            attributes = [],
            **keywords
    ):
        Component.__init__(self = self, name = name, mechanisms=mechanisms, parameters = parameters, attributes = attributes, **keywords)

    def update_species(self):
        species = [crn.specie(self.name, type="complex", attributes=self.attributes)]
        return species

    def update_reactions(self):
        return []