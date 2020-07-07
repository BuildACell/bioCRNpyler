#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from warnings import warn
from .sbmlutil import *
import warnings
import numpy as np
from typing import List, Union, Dict



class Species(object):
    """ A formal species object for a CRN
     A Species must have a name. They may also have a material_type (such as DNA,
     RNA, Protein), and a list of attributes.
    """

    def __init__(self, name: str, material_type="", attributes=[],
                 initial_concentration=0):

        self.name = self.check_name(name)
        self.material_type = self.check_material_type(material_type)
        self.initial_concentration = initial_concentration
        if material_type == "complex":
            warn("species which are formed of two species or more should be "
                 "called using the chemical_reaction_network.ComplexSpecies "
                 "constructor for attribute inheritance purposes.")

        self.attributes = []

        if attributes is not None:
            for attribute in attributes:
                self.add_attribute(attribute)

     #Check that the string contains is alpha-numeric characters or "_" and that the first character is a letter. IF the name is a starts with a number, there must be a material type.
    def check_material_type(self, material_type):

        
        if material_type in [None, ""] and self.name[0].isnumeric():
            raise ValueError(f"species name: {self.name} contains a number as the first character and therefore requires a material_type.")
        elif material_type == None:
            return ""
        elif (material_type.replace("_", "").isalnum() and material_type.replace("_", "")[0].isalpha()) or material_type == "":
                return material_type
        else:
            raise ValueError(f"material_type {material_type} must be alpha-numeric and start with a letter.")

    
    #Check that the string contains only underscores and alpha-numeric characters
    def check_name(self, name):
        no_underscore_string = name.replace("_", "")
        if no_underscore_string.isalnum():
            return name
        else:
            raise ValueError(f"name attribute {name} must consist of letters, numbers, or underscores.")

    def __repr__(self):
        txt = ""
        if self.material_type not in ["", None]:
            txt = self.material_type + "_"

        txt += self.name

        if len(self.attributes) > 0 and self.attributes != []:
            for i in self.attributes:
                if i is not None:
                    txt += "_" + str(i)
        txt.replace("'", "")
        return txt

    def replace_species(self, species, new_species):
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if self == species:
            return new_species
        else:
            return self

    #A more powerful printing function
    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        txt = ""
        if self.material_type not in ["", None] and show_material:
            txt = self.material_type + "["

        txt += self.name

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")

        if self.material_type not in ["", None] and show_material:
            txt += "]"

        return txt

    def add_attribute(self, attribute: str):
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(), "Attribute: %s must be an alpha-numeric string" % attribute
        self.attributes.append(attribute)

    def __eq__(self, other):
        """
        Overrides the default implementation
        Two species are equivalent if they have the same name, type, and attributes
        :param other: Species instance
        :return: boolean
        """

        if isinstance(other, Species) \
                            and self.material_type == other.material_type \
                            and self.name == other.name \
                            and set(self.attributes) == set(other.attributes):
            return True
        else:
            return False
    def __gt__(self,Species2):
        return self.name > Species2.name
    def __lt__(self,Species2):
        return self.name < Species2.name

    def __hash__(self):
        return str.__hash__(repr(self))


class ComplexSpecies(Species):
    """ A special kind of species which is formed as a complex of two or more species.
        Used for attribute inheritance and storing groups of bounds Species. 
        Note taht in a ComplexSpecies, the order of the species list does not matter.
        This means that ComplexSpecies([s1, s2]) = ComplexSpecies([s2, s1]). 
        This is good for modelling order-indpendent binding complexes.
        For a case where species order matters (e.g. polymers) use OrderedComplexSpecies
    """
    def __init__(self, species: List[Union[Species,str]], name = None, material_type = "complex", attributes = None, initial_concentration = 0, **keywords):
        if len(species) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")


        self.species = []
        for s in  flatten_list(species):
            if isinstance(s, Species):
                self.species.append(s)
            elif isinstance(s, str):
                self.species.append(Species(s))
            else:
                raise ValueError("ComplexSpecies must be defined by (nested) list of Species (or subclasses thereof).")

        
        self.species_set = list(set(self.species))

        if name is not None:
            self.custom_name = True
        elif name is None:
            self.custom_name = False
            name = ""
            list.sort(self.species, key = lambda s:repr(s))
            
            list.sort(self.species_set, key = lambda s:repr(s))
            for s in self.species_set:
                count = self.species.count(s)
                if count > 1:
                    name+=f"{count}x_"
                if not (isinstance(s, ComplexSpecies) or s.material_type == ""):
                    name+=f"{s.material_type}_{s.name}_"
                else:
                    name+=f"{s.name}_"
            name = name[:-1]

        self.name = self.check_name(name)
        self.material_type = self.check_material_type(material_type)
        self.initial_concentration = initial_concentration

        if attributes is None:
            attributes = []
        for s in self.species:
            attributes += s.attributes
        attributes = list(set(attributes))

        while None in attributes:
            attributes.remove(None)

        self.attributes = attributes


    def __contains__(self,item):
        if not isinstance(item, Species):
            raise ValueError("Operator 'in' requires chemical_reaction_network.Species (or a subclass). Received: "+str(item))
        if item in self.species:
            #this is the base case
            return True
        else:
            #this is the recursive part. We want to check all
            #internal complexes for the thing we're looking for
            for content in self.species:
                if isinstance(content,ComplexSpecies) :
                    if item in content:
                        return True
            #if we got here then we've failed to find it
            return False

    #Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
    def replace_species(self, species: Species, new_species: Species):
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            if s == species:
                new_species_list.append(new_species)
            elif isinstance(s, ComplexSpecies):
                new_s = s.replace_species(species, new_species)
                new_species_list.append(new_s)
            else:
                new_species_list.append(s)

        new_name = None
        if self.custom_name == True:
            new_name = self.name
        
        return ComplexSpecies(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes)

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        txt = ""
        if self.material_type not in ["", None] and show_material:
            txt += self.material_type
        txt += "["
        for s in self.species_set:
            count = self.species.count(s)
            if count > 1:
                txt += f"{count}x_"
            txt += s.pretty_print(show_material = show_material, show_attributes = False)+":"
        txt = txt[:-1]

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")

        txt += "]"

        return txt


class Multimer(ComplexSpecies):
    """A subclass of ComplexSpecies for Complexes made entirely of the same kind of species,
    eg dimers, tetramers, etc.
    """
    def __init__(self, species, multiplicity, name = None, material_type = "complex", attributes = None, initial_concentration = 0):

        if isinstance(species, str):
            species = [Species(name = species)]
        elif not isinstance(species, Species):
            raise ValueError("Multimer must be defined by a Species (or subclasses thereof) and a multiplicity (int).")
        else:
            species = [species]

        ComplexSpecies.__init__(self, species = species*multiplicity, name = name, material_type = material_type, attributes = attributes, initial_concentration = initial_concentration)   

class OrderedComplexSpecies(ComplexSpecies):
    """ A special kind of species which is formed as a complex of two or more species.
        In OrderedComplexSpecies the order in which the complex subspecies are is defined
        denote different species, eg [s1, s2, s3] != [s1, s3, s2].
        Used for attribute inheritance and storing groups of bounds Species. 
    """

    def __init__(self, species, name = None, material_type = "ordered_complex", attributes = None, initial_concentration = 0):
        if len(species) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")


        new_species = flatten_list(species)
        self.species = []
        for s in new_species:
            if isinstance(s, Species):
                self.species.append(s)
            elif isinstance(s, str):
                self.species.append(Species(s))
            else:
                raise ValueError("OrderedComplexSpecies must be defined by (nested) list of Species (or subclasses thereof).")

        if name is not None:
            self.custom_name = True
        elif name is None:
            self.custom_name = False
            name = ""
            for s in self.species:
                if isinstance(s, str):
                    s = Species(name = s)
                if s.material_type not in ["complex", "ordered_complex", ""]:
                    name+=f"{s.material_type}_{s.name}_"
                else:
                    name+=f"{s.name}_"
            name = name[:-1]

        self.name = self.check_name(name)
        self.material_type = self.check_material_type(material_type)
        self.initial_concentration = initial_concentration

        if attributes is None:
            attributes = []
        for s in self.species:
            attributes += s.attributes
        attributes = list(set(attributes))

        while None in attributes:
            attributes.remove(None)

        self.attributes = attributes

    #Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
    def replace_species(self, species: Species, new_species: Species):
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            if s == species:
                new_species_list.append(new_species)
            elif isinstance(s, ComplexSpecies):
                new_s = s.replace_species(species, new_species)
                new_species_list.append(new_s)
            else:
                new_species_list.append(s)

        new_name = None
        if self.custom_name == True:
            new_name = self.name
        
        return OrderedComplexSpecies(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes)

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        txt = ""
        if self.material_type not in ["", None] and show_material:
            txt += self.material_type+"["
        
        txt += "["

        for s in self.species:
            txt += s.pretty_print(show_material = show_material, show_attributes = False)+":"
        txt = txt[:-1]

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")
        txt += "]"

        return txt


class Reaction(object):
    """ An abstract representation of a chemical reaction in a CRN
    A reaction has the form:
       \sum_i n_i I_i --> \sum_i m_i O_i @ rate = k
       where n_i is the count of the ith input, I_i, and m_i is the count of the
       ith output, O_i.
    If the reaction is reversible, the reverse reaction is also included:
       \sum_i m_i O_i  --> \sum_i n_i I_i @ rate = k_rev
    """
    def __init__(self, inputs, outputs, k = 0, input_coefs = None,
                 output_coefs = None, k_rev = 0, propensity_type = "massaction",
                 rate_formula = None, propensity_params = None):

        if len(inputs) == 0 and len(outputs) == 0:
            warn("Reaction Inputs and Outputs both contain 0 Species.")

        if k != 0 and propensity_params is not None and "k" not in propensity_params:
            propensity_params["k"] = k
        elif k == 0 and propensity_params is not None and "k" in propensity_params:
            k = propensity_params["k"]
        elif k != 0 and propensity_params is not None and k != propensity_params['k']:
            print("k=", k, "propensity_params[k]", propensity_params["k"], "propensity_params", propensity_params)
            raise ValueError("Inconsistent rate constants: propensity_params['k'] != k.")

        if propensity_type == "massaction" and propensity_params is not None:
            warn("ValueWarning: propensity_params dictionary passed into a "
                 "massaction propensity. Massaction propensities do not "
                 "require a param dictionary.")
        elif propensity_type != "massaction" and propensity_params is None:
            raise ValueError("Non-massaction propensities require a propensity_params dictionary passed to the propensity_params keyword.")
        elif propensity_type != "massaction" and k_rev != 0:
            raise ValueError("Invalid reversible reaction for propensity "
                             f"type = {propensity_type}. Only massaction "
                             "propensities support the reversible rate k_rev. "
                             "Consider creating two seperate reactions "
                             "instead.")
        elif propensity_type == "hillpositive":
            if not ("k" in propensity_params and "s1" in propensity_params and "K" in propensity_params \
                    and "n" in propensity_params):
                raise ValueError("hillpositive propensities, p(s1; k, K, n) "
                        "= k*s1^n/(s1^n + K), require the following "
                        "propensity_params: "
                        "'k':rate constant (float)"
                        "'s1':species (chemical_reaction_network.species), "
                        "'n':cooperativity(float), "
                        "and 'K':dissociationc constant (float).")
        elif propensity_type == "hillnegative":
            if not ("k" in propensity_params and "s1" in propensity_params and "K" in propensity_params \
                    and "n" in propensity_params):
                raise ValueError("hillnegative propensities, "
                        "p(s1; k, K, n) = k*1/(s1^n + K), require "
                        "the following propensity_params:"
                        "'k':rate constant (float)"
                        "'s1':species (chemical_reaction_network.species), "
                        "'n':cooperativity(float), "
                        "and 'K':dissociationc constant (float)")
        elif propensity_type == "proportionalhillpositive":
            if not ("k" in propensity_params and "s1" in propensity_params and "d" in propensity_params \
                    and "K" in propensity_params \
                    and "n" in propensity_params):
                raise ValueError("proportionalhillpositive propensities, "
                    "p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K), require the "
                    "following propensity_params: "
                    "'k':rate constant (float)"
                    "'s1':species (chemical_reaction_network.species), "
                    "'d':species (chemical_reaction_network.species), "
                    "'n':cooperativity(float), "
                    "and 'K':dissociationc onstant (float)")
        elif propensity_type == "proportionalhillnegative":
            if not ("k" in propensity_params and "s1" in propensity_params and "d" in propensity_params \
                    and "K" in propensity_params \
                    and "n" in propensity_params):
                raise ValueError("proportionalhillnegative propensities, "
                    "p(s1, d; k, K, n) = k*d/(s1^n + K), require the "
                    "following propensity_params: "
                    "'k':rate constant (float)"
                    "'s1':species (chemical_reaction_network.species), "
                    "'d':species (chemical_reaction_network.species), "
                    "'n':cooperativity(float), "
                    "and 'K':dissociationc onstant (float)")
        elif propensity_type == "general":
            if "rate" not in propensity_params:
                raise ValueError("general propensities, p(s) = k * f(s), "
                    "require the propensity_params: "
                    "'rate':f(s) where f(s) is an SBML compatable function "
                    "of arbitrary species, "
                    "s (use repr(chemical_reaction_network.species) to get "
                    "the proper text representation of a species name).")
        elif propensity_type != "massaction":
            raise ValueError(f"Unknown propensity type: {propensity_type}.")
        self.propensity_type = propensity_type
        self.propensity_params = propensity_params

        # Check that inputs and outputs only contain species
        #if inputs or outputs is a nested list, flatten that list
        new_inputs = []
        for s in flatten_list(inputs):
            if isinstance(s, Species):
                new_inputs.append(s)
            else:
                raise ValueError(f"A non-species object was used as a species: {s}!")
        inputs = new_inputs

        new_outputs = []
        for s in flatten_list(outputs):
            if isinstance(s, Species):
                new_outputs.append(s)
            else:
                raise ValueError(f"A non-species object was used as a species: {s}!")
        outputs = new_outputs

        #OLD CHECK
        # Check that inputs and outputs only contain species
        #if any(not isinstance(s, Species) for s in inputs + outputs):
        #    raise ValueError("A non-species object was used as a species.")

        # internal representation of a reaction

        #self.inputs and self.outputs should be ordered lists.
        self.inputs = []
        for s in inputs:
            if s not in self.inputs:
                self.inputs.append(s)
        self.outputs = []
        for s in outputs:
            if s not in self.outputs:
                self.outputs.append(s)

        #self.input_coefs[i] is the number of self.inputs[i] into the reaction
        self.input_coefs = None
        #self.output coefs is analogous to above
        self.output_coefs = None

        # Check that rates are valid
        if k <= 0:
            raise ValueError(f"Reaction rate <= 0: k={k}")
        else:
            self.k = k
        if k_rev > 0:
            self.reversible = True
            self.k_r = k_rev
        else:
            self.k_r = 0
            self.reversible = False

        # TODO input coefficients should be stored with the species a dictionary (same for the output )
        # Set input coefficients
        if input_coefs is None:
            self.input_coefs = [inputs.count(s) for s in self.inputs]
        elif input_coefs is not None and len(input_coefs) == len(self.inputs):
            self.input_coefs = input_coefs
        elif len(input_coefs) == len(inputs) \
             and len(self.inputs) != len(inputs):
            raise ValueError("Input species and input_coefs contain "
                             "contradictory counts.")
        else:
            raise ValueError(f"len(input_coefs) ({len(input_coefs)}) doesn't "
                             f"match len(self.inputs) ({len(self.inputs)}).")

        # Set Output Coefs
        if output_coefs is None:
            self.output_coefs = [outputs.count(s) for s in self.outputs]
        elif output_coefs is not None \
             and len(output_coefs) == len(self.outputs):
            self.output_coefs = output_coefs
        elif len(output_coefs) == len(outputs) \
             and len(self.outputs) != len(outputs):
            raise ValueError(f"Output species ({self.outputs}) and output_coefs ({output_coefs}) contain "
                             "contradictory counts.")
        else:
            raise ValueError(f"len(output_coefs) ({len(output_coefs)}) doesn't "
                             f"match len(self.outputs) ({len(self.outputs)}).")


    #Replaces species with new_species in the entire CRN
    def replace_species(self, species: Species, new_species: Species):
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_inputs = []
        for s in self.inputs:
            if s == species:
                new_inputs.append(new_species)
            elif isinstance(s, ComplexSpecies):
                new_s = s.replace_species(species, new_species)
                new_inputs.append(new_s)
            else:
                new_inputs.append(s)
        self.inputs = new_inputs

        new_outputs = []
        for s in self.outputs:
            if s == species:
                new_outputs.append(new_species)
            elif isinstance(s, ComplexSpecies):
                new_s = s.replace_species(species, new_species)
                new_outputs.append(new_s)
            else:
                new_outputs.append(s)
        self.outputs = new_outputs

        if self.propensity_params is not None:
            new_params = {}
            for key in self.propensity_params:
                if isinstance(self.propensity_params[key], ComplexSpecies):
                    new_params[key] = self.propensity_params[key].replace_species(species, new_species)
                elif isinstance(self.propensity_params[key], Species) and self.propensity_params[key] == species:
                    new_params[key] = new_species
                else:
                    new_params[key] = self.propensity_params[key]

        new_r = Reaction(inputs = self.inputs, outputs = self.outputs, input_coefs = self.input_coefs, output_coefs = self.output_coefs, propensity_type = self.propensity_type, propensity_params = self.propensity_params, k = self.k, k_rev = self.k_r)
        return new_r

    #Helper function to print the text of a rate function
    def rate_func_text(self, pretty_print = False,  show_material = True, show_attributes = True, **kwargs):
        tab = (" " * 8)
        txt = ""
        if self.propensity_type == "massaction":
            input_func_args = ""
            input_prod = f"{self.k}"
            for i in range(len(self.inputs)):
                if pretty_print:
                    sin = self.inputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
                else:
                    sin = repr(self.inputs[i])

                input_func_args += f"{sin}"

                if self.input_coefs[i] > 1:
                    input_prod+=f"*{sin}^{self.input_coefs[i]}"
                else:
                    input_prod+=f"*{sin}"

                if i < len(self.inputs)-1:
                    input_func_args += ","

            if len(self.inputs) > 0:
                input_func_args = "("+input_func_args+")"
            txt += f"massaction: k_f{input_func_args}={input_prod}"

            if self.reversible:
                output_func_args = ""
                output_prod = f"{self.k_r}"
               
                for i in range(len(self.outputs)):
                    if pretty_print:
                        sout = self.outputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
                    else:
                        sout = repr(self.outputs[i])

                    output_func_args += f"{sout}"

                    if self.output_coefs[i] > 1:
                        output_prod+=f"*{sout}^{self.output_coefs[i]}"
                    else:
                        output_prod+=f"*{sout}"

                    if i < len(self.outputs)-1:
                        output_func_args += ","

                if len(self.outputs) > 0:
                    output_func_args = "("+output_func_args+")"
                txt += f" k_r{output_func_args}={output_prod}"

        elif self.propensity_type == "hillpositive":
            if pretty_print:
                s1 = self.propensity_params["s1"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                s1 = repr(self.propensity_params["s1"])

            kd = str(self.propensity_params["K"])
            n = str(self.propensity_params["n"])
            txt += f"hillpositive: k({s1})={self.k}*{s1}^{n}/({kd}+{s1}^{n})"
        elif self.propensity_type == "hillnegative":

            if pretty_print:
                s1 = self.propensity_params["s1"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                s1 = repr(self.propensity_params["s1"])

            kd = str(self.propensity_params["K"])
            n = str(self.propensity_params["n"])
            txt += f"hillnegative: k({s1})={self.k}*1/({kd}+{s1}^{n})"
        elif self.propensity_type == "proportionalhillpositive":
            if pretty_print:
                s1 = self.propensity_params["s1"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
                s2 = self.propensity_params["d"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                s1 = repr(self.propensity_params["s1"])
                s2 = repr(self.propensity_params["d"])

            kd = str(self.propensity_params["K"])
            n = str(self.propensity_params["n"])
            txt += (f"proportionalhillpositive: k({s1}, "
                   f"{s2})={self.k}*{s2}*{s1}^{n}/({kd}+{s1}^{n})")
        elif self.propensity_type == "proportionalhillnegative":
            if pretty_print:
                s1 = self.propensity_params["s1"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
                s2 = self.propensity_params["d"].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                s1 = repr(self.propensity_params["s1"])
                s2 = repr(self.propensity_params["d"])

            kd = str(self.propensity_params["K"])
            n = str(self.propensity_params["n"])
            txt += (f"proportionalhillnegative: k({s1}, "
                   f"{s2})={self.k}*{s2}/({kd}+{s1}^{n})")
        elif self.propensity_type == "general":
            eq = self.propensity_params["rate"]
            txt += f"general: k(x)={self.k}*{eq}"
        else:
            raise ValueError("Unknown Propensity Type: "
                             f"{self.propensity_type}.")

        return txt

    def __repr__(self, **kwargs):
        tab = (" " * 8)
        txt = ""
        for i in range(len(self.inputs)):
            if self.input_coefs[i] > 1:
                txt += str(self.input_coefs[i]) + " " + str(self.inputs[i])
            else:
                txt += str(self.inputs[i])
            if i < len(self.inputs) - 1:
                txt += " + "
        if self.reversible:
            txt += " <--> "
        else:
            txt += " --> "
        for i in range(len(self.outputs)):
            if self.output_coefs[i] > 1:
                txt += str(self.output_coefs[i]) + " " + str(self.outputs[i])
            else:
                txt += str(self.outputs[i])
            if i < len(self.outputs) - 1:
                txt += " + "
        txt += tab
        txt += self.rate_func_text(**kwargs)
        
        return txt

    def pretty_print(self, show_rates = True, show_material = True, show_attributes = True, **kwargs):
        tab = (" " * 8)
        txt = ""
        for i in range(len(self.inputs)):
            if self.input_coefs[i] > 1:
                txt += str(self.input_coefs[i]) + " " + self.inputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                txt += self.inputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            if i < len(self.inputs) - 1:
                txt += " + "
        if self.reversible:
            txt += " <--> "
        else:
            txt += " --> "
        for i in range(len(self.outputs)):
            if self.output_coefs[i] > 1:
                txt += str(self.output_coefs[i]) + " " + self.outputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            else:
                txt += self.outputs[i].pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs)
            if i < len(self.outputs) - 1:
                txt += " + "
        txt += tab
        if show_rates:
            if self.reversible:
                rate_txt = self.rate_func_text(pretty_print = True, show_material = show_material, show_attributes = show_attributes, **kwargs)
                rate_txt = rate_txt.replace(" k_r", "\n"+tab+"k_r")
                txt += "\n"+tab+rate_txt
            else:
                txt += "\n"+tab+self.rate_func_text(pretty_print = True, show_material = show_material, show_attributes = show_attributes, **kwargs)

        return txt


    def __eq__(self, other):
        """Overrides the default implementation.
           Two reactions are equivalent if they have the same inputs, outputs,
           and rates."""
        complexes_equal = Reaction.complex_set_equality(self.inputs,
                                                    self.input_coefs,
                                                    other.inputs,
                                                    other.input_coefs) \
                           and Reaction.complex_set_equality(self.outputs,
                                                         self.output_coefs,
                                                         other.outputs,
                                                         other.output_coefs)
        rates_equal = (other.k == self.k and other.k_r == self.k_r)
        propensity_types_equal = (self.propensity_type == other.propensity_type)

        # must both be reactions with the same rates and numbers of inputs and
        # outputs.
        if not isinstance(other, Reaction):
            return False
        if complexes_equal and rates_equal and propensity_types_equal:
            return True
        elif complexes_equal and propensity_types_equal:
            #warn("Two reactions with the same inputs and outputs but different "
                 #"rates are formally different, but may be undesired:"
                 #f"{repr(self)} and {repr(other)}.")
            return False

        # If the reactions are reversible inverses of eachother, one's forward
        # reaction could be the other's reverse
        elif self.reversible and other.reversible:
            reverse_complex_equal = Reaction.complex_set_equality(self.inputs,
                                                            self.input_coefs,
                                                            other.outputs,
                                                            other.output_coefs)\
                        and Reaction.complex_set_equality(self.outputs,
                                                      self.output_coefs,
                                                      other.inputs,
                                                      other.input_coefs)
            reverse_rates_equal = (other.k == self.k_r and other.k_r == self.k)
            if reverse_complex_equal and reverse_rates_equal:
                return True
            elif reverse_complex_equal:
                warn("Two reversible reactions with the same inputs and outputs"
                    " (reversed) but different rates are formally equal, but "
                    f"may be undesired:{repr(self)} and {repr(other)}")
                return True
            else:
                return False
        else:
            return False

    @staticmethod
    def complex_set_equality(c1, c1_coefs, c2, c2_coefs):
        """Checks to see if two formal complexes (reaction input or output sets) are equal."""
        if len(c1) != len(c2):
            return False
        else:
            for i in range(len(c1)):
                s1 = c1[i]
                coef1 = c1_coefs[i]
                if s1 not in c2 or coef1 != c2_coefs[c2.index(s1)]:
                    return False
        return True

    def pyrepr(self):
        if self.reversible:
            return [
                ([repr(i) for i in self.inputs], self.input_coefs,
                 [repr(i) for i in self.outputs], self.output_coefs,
                 self.k),
                ([repr(i) for i in self.outputs], self.output_coefs,
                 [repr(i) for i in self.inputs], self.input_coefs,
                 self.k_r)]
        else:
            return [([repr(i) for i in self.inputs], self.input_coefs,
                     [repr(i) for i in self.outputs],
                     self.output_coefs, self.k)]


class ChemicalReactionNetwork(object):
    """ A chemical reaction network is a container of species and reactions
    chemical reaction networks can be compiled into SBML or represented
    conveniently as python tuple objects.
    reaction types:
       mass action: standard mass action semantics where the propensity of a
                reaction is given by deterministic propensity =
                        k \Prod_{inputs i} [S_i]^a_i
               stochastic propensity =
                        k \Prod_{inputs i} (S_i)!/(S_i - a_i)!
               where a_i is the spectrometric coefficient of species i
    """
    def __init__(self, species: List[Species], reactions: List[Reaction], warnings = False):
        self.species, self.reactions = ChemicalReactionNetwork.check_crn_validity(reactions, species, warnings=warnings)

        # TODO check whether we need this data structure
        self.species2index = {}
        for i in range(len(self.species)):
            self.species2index[str(self.species[i])] = i

    def add_species(self, species, warnings = False):
        if not isinstance(species, list):
            species = [species]
        self.species, self.reactions = ChemicalReactionNetwork.check_crn_validity(self.reactions, self.species+species, warnings=warnings)

    def add_reactions(self, reactions, warnings = False):
        if not isinstance(reactions, list):
            reactions = [reactions]
        self.species, self.reactions = ChemicalReactionNetwork.check_crn_validity(self.reactions+reactions, self.species, warnings=warnings)

    @staticmethod
    def check_crn_validity(reactions: List[Reaction], species: List[Species], warnings = False):
        # Check to make sure species are valid and only have a count of 1
        checked_species = []
        if not all(isinstance(s, Species) for s in species):
            print(f"A non-species object was used as a species: {species}!")
            raise ValueError("A non-species object was used as a species!")

        for s in species:
            if species.count(s) > 1:
                pass
                #warn("Species "+str(s)+" duplicated in CRN definition.
                # Duplicates have been removed.")
            if s not in checked_species:
                checked_species.append(s)

        # Check to make sure reactions are valid meaning:
        #   only have a count of 1
        #   all species in the inputs/outputs are also in the species list
        checked_reactions = []

        if not all(isinstance(r, Reaction) for r in reactions):
            raise ValueError("A non-reaction object was used as a reaction!")

        for r in reactions:
            if reactions.count(r) > 1:
                pass

            checked_reactions.append(r)
            #if r not in checked_reactions:
            #    checked_reactions.append(r)

            for s in r.inputs:
                if s not in checked_species and warnings:
                    warn(f"Reaction {repr(r)} contains a species {repr(s)} "
                         "which is not in the CRN.")

            for s in r.outputs:
                if s not in checked_species and warnings:
                    warn(f"Reaction {repr(r)} contains a species {repr(s)} "
                         "which is not in the CRN.")

        return checked_species, checked_reactions

    def __repr__(self):
        txt = "Species = "
        for s in self.species:
            txt += repr(s) + ", "
        txt = txt[:-2] + '\n'
        txt += "Reactions = [\n"

        for r in self.reactions:
            txt += "\t" + repr(r) + "\n"
        txt += "]"
        return txt

    def pretty_print(self, show_rates = True, show_material = True, show_attributes = True, **kwargs):
        txt = f"Species ({len(self.species)}) = "+"{"
        for sind in range(len(self.species)):
            s = self.species[sind]
            txt += f"{sind}. "+s.pretty_print(show_material = show_material, show_attributes = show_attributes, **kwargs) + ", "
        txt = txt[:-2] + '}\n'
        txt += f"Reactions ({len(self.reactions)}) = [\n"

        for rind in range(len(self.reactions)):
            r = self.reactions[rind]
            txt += f"{rind}. " + r.pretty_print(show_rates = show_rates, show_material = show_material, show_attributes = show_attributes, **kwargs) + "\n"
        txt += "]"
        return txt

    def pyrepr(self):
        reactions = []
        for r in self.reactions:
            reactions += r.pyrepr()
        species = [str(s) for s in self.species]
        return species, reactions

    # TODO check whether we need this method
    def species_index(self, species: Species):
        if len(self.species2index) != len(self.species):
            self.species2index = {}
            for i in range(len(self.species)):
                self.species2index[str(self.species[i])] = i
        return self.species2index[str(species)]

    def initial_condition_vector(self, init_cond_dict: Dict[str,float]):
        x0 = [0.0] * len(self.species)
        for idx, s in enumerate(self.species):
            if s in init_cond_dict:
                x0[idx] = init_cond_dict[s]
        return x0

    def get_all_species_containing(self, species: Species, return_as_strings = False):
        """Returns all species (complexes and otherwise) containing a given species
           (or string).
        """
        return_list = []
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        for s in self.species:
            if species == s or (isinstance(s, ComplexSpecies) and species in s.species):
                if return_as_strings:
                    return_list.append(repr(s))
                else:
                    return_list.append(s)
        return return_list

    #Replaces species with new_species in the entire CRN
    def replace_species(self, species: Species, new_species: Species):
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            if s == species:
                new_species_list.append(new_species)
            elif isinstance(s, ComplexSpecies):
                new_s = s.replace_species(species, new_species)
                new_species_list.append(new_s)
            else:
                new_species_list.append(s)

        new_reaction_list = []

        for r in self.reactions:
            new_r = r.replace_species(species, new_species)
            new_reaction_list.append(new_r)


        return ChemicalReactionNetwork(new_species_list, new_reaction_list)


    def generate_sbml_model(self, stochastic_model=False, **keywords):
        document, model = create_sbml_model(**keywords)

        for s in self.species:

            add_species(model=model, compartment=model.getCompartment(0),
                    species=s, initial_concentration=s.initial_concentration)

        rxn_count = 0
        for r in self.reactions:
            rxn_id = "r" + str(rxn_count)
            add_reaction(model, r.inputs, r.input_coefs, r.outputs,
                         r.output_coefs, rxn_id, r.k,
                         stochastic = stochastic_model,
                         propensity_type=r.propensity_type,
                         propensity_params = r.propensity_params)
            rxn_count += 1

            if r.reversible and r.propensity_type == "massaction":
                add_reaction(model, r.outputs, r.output_coefs, r.inputs,
                             r.input_coefs, rxn_id, r.k_r,
                             stochastic=stochastic_model,
                             propensity_type=r.propensity_type)
                rxn_count += 1

        if document.getNumErrors():
            warn('SBML model generated has errors. Use document.getErrorLog() to print all errors.')
        return document, model

    def write_sbml_file(self, file_name=None, **keywords):
        document, _ = self.generate_sbml_model(**keywords)
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, 'w') as f:
            f.write(sbml_string)
        return True

    def create_bioscrape_model(self):
        from bioscrape.types import Model

        species_list = []
        initial_condition_dict = {}
        for s in self.species:
            species_list.append(repr(s))
            if s.initial_concentration is None:
                initial_condition_dict[repr(s)] = 0
            else:
                initial_condition_dict[repr(s)] = s.initial_concentration

        reaction_list = []
        reaction_counter = 0
        rate_list = []
        for rxn in self.reactions:

            reactants = []
            for i in range(len(rxn.inputs)):
                reactants += [repr(rxn.inputs[i])]*int(rxn.input_coefs[i])
            products = []
            for i in range(len(rxn.outputs)):
                products += [repr(rxn.outputs[i])]*int(rxn.output_coefs[i])

            prop_type = rxn.propensity_type
            if rxn.propensity_params is None:
                prop_params = {}
            else:
                prop_params = {}
                for k in rxn.propensity_params:
                    v = rxn.propensity_params[k]
                    if isinstance(v, Species):
                        prop_params[k] = repr(v)
                    elif isinstance(v, str):
                        prop_params[k] = v
                    else:
                        prop_params[k] = float(v)


            prop_params['propensity_type'] = rxn.propensity_type
            prop_params['k'] = rxn.k

            reaction_list.append((reactants, products, prop_type,
                                  dict(prop_params)))

            if rxn.reversible and rxn.propensity_type == "massaction":
                prop_params['k'] = rxn.k_r
                reaction_list.append((products, reactants, prop_type,
                                      dict(prop_params)))
            elif rxn.reversible:
                raise ValueError("Only massaction irreversible reactions are "
                                 "supported for automatic bioscrape simulation."
                                 " Consider creating two seperate reactions.")
        model = Model(species = species_list, reactions = reaction_list,
                      initial_condition_dict = initial_condition_dict)
        return model

    def simulate_with_bioscrape(self, timepoints, initial_condition_dict = {},
                                stochastic = False, return_dataframe = True,
                                safe = False, **kwargs):

        """Simulate CRN model with bioscrape (https://github.com/biocircuits/bioscrape).
        Returns the data for all species as Pandas dataframe.
        """
        result = None
        try:
            from bioscrape.simulator import py_simulate_model
            m = self.create_bioscrape_model()
            m.set_species(initial_condition_dict)
            if not stochastic and safe:
                safe = False
            result = py_simulate_model(timepoints, Model = m,
                                        stochastic = stochastic,
                                        return_dataframe = return_dataframe,
                                        safe = safe)
        except ModuleNotFoundError:
            warnings.warn('bioscrape was not found, please install bioscrape')

        return result

    def simulate_with_bioscrape_via_sbml(self, timepoints, file = None,
                initial_condition_dict = {}, return_dataframe = True,
                stochastic = False, **kwargs):

        """Simulate CRN model with bioscrape via writing a SBML file temporarily.
        [Bioscrape on GitHub](https://github.com/biocircuits/bioscrape).

        Returns the data for all species as Pandas dataframe.
        """
        result = None
        m = None
        try:
            from bioscrape.simulator import py_simulate_model
            if file is None:
                self.write_sbml_file(file_name ="temp_sbml_file.xml")
                file_name = "temp_sbml_file.xml"
            elif isinstance(file, str):
                file_name = file
            else:
                file_name = file.name

            if 'sbml_warnings' in kwargs:
                sbml_warnings = kwargs.get('sbml_warnings')
            else:
                sbml_warnings = False
            m = bioscrape.types.Model(sbml_filename = file_name, sbml_warnings = sbml_warnings)
            # m.write_bioscrape_xml('temp_bs'+ file_name + '.xml') # Uncomment if you want a bioscrape XML written as well.
            m.set_species(initial_condition_dict)
            result = py_simulate_model(timepoints, Model = m,
                                                stochastic = stochastic,
                                                return_dataframe = return_dataframe)
        except ModuleNotFoundError:
            warnings.warn('bioscrape was not found, please install bioscrape')

        return result, m

    def runsim_roadrunner(self, timepoints, filename, species_to_plot = []):

        """ To simulate using roadrunner

        Arguments:
        timepoints: The array of time points to run the simulation for. 
        filename: Name of the SBML file to simulate

        Returns the results array as returned by RoadRunner.

        Refer to the libRoadRunner simulator library documentation 
        for details on simulation results: (http://libroadrunner.org/)[http://libroadrunner.org/]
        NOTE : Needs roadrunner package installed to simulate.
        """
        res_ar = None
        try:
            import roadrunner

            rr = roadrunner.RoadRunner(filename)
            result = rr.simulate(timepoints[0],timepoints[-1],len(timepoints))
            res_ar = np.array(result)
        except ModuleNotFoundError:
            warnings.warn('libroadrunner was not found, please install libroadrunner')
        return res_ar

#Helper function to flatten lists
def flatten_list(in_list):
    out_list = []
    for element in in_list:
        if isinstance(element, list):
            out_list += flatten_list(element)
        else:
            out_list += [element]
    return out_list
