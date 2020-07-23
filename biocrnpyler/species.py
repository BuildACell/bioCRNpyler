#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .sbmlutil import *
import warnings
import numpy as np
from typing import List, Union, Dict
import itertools
import copy


class Species(object):
    """ A formal species object for a CRN
     A Species must have a name. They may also have a material_type (such as DNA,
     RNA, Protein), and a list of attributes.
    """

    def __init__(self, name: str, material_type="", attributes: Union[List,None] = None,
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
            if not isinstance(attributes,list):
                attributes = list(attributes)
            for attribute in attributes:
                self.add_attribute(attribute)

    def check_material_type(self, material_type):
        """Check that the string contains is alpha-numeric characters or "_" and that the first character is a letter.
        Material types cannot contain double underscores "__" or end with an underscore

        If the name is a starts with a number, there must be a material type.
        """
        if material_type in [None, ""] and self.name[0].isnumeric():
            raise ValueError(f"species name: {self.name} contains a number as the first character and therefore requires a material_type.")
        elif material_type is None:
            return ""
        elif (material_type.replace("_", "").isalnum()  
            and material_type.replace("_", "")[0].isalpha() and "__" not in material_type and material_type[len(material_type)-1] != "_") or material_type == "":

                return material_type
        else:
            raise ValueError(f"material_type {material_type} must be alpha-numeric, start with a letter, cannot start or end with an underscore and cannot contain a double underscore.")

    def check_name(self, name):
        """
        Check that the string contains only underscores and alpha-numeric characters
        names cannot contain double underscores "__" or end with an underscore
        """
        no_underscore_string = name.replace("_", "")
        if no_underscore_string.isalnum() and "__" not in name and name[len(name)-1] != "_":
            return name
        else:
            raise ValueError(f"name attribute {name} must consist of letters, numbers, or underscores, cannot start or end with an underscore or contain a double underscore.")

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

    def get_species(self, **kwargs):
        """Used in some recursive calls where ComplexSpecies returns a list and Species will return just themselves (in a list)
        """
        return [self]

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        """
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
        """Adds attributes to a Species
        """
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(), "Attribute: %s must be an alpha-numeric string" % attribute
        self.attributes.append(attribute)

    def __eq__(self, other):
        """Two species are equivalent if they have the same name, type, and attributes

        :param other: Species instance
        :return: boolean
        """
        if other.__class__ is self.__class__:
            return (self.material_type, self.name, set(self.attributes)) \
                   == (other.material_type, other.name, set(other.attributes))
        else:
            return False

    def __gt__(self, Species2):
        return self.name > Species2.name

    def __lt__(self, Species2):
        return self.name < Species2.name

    def __hash__(self):
        return str.__hash__(repr(self))

    def __contains__(self, item):
        if item == self:
            return True
        return False

    def __rmul__(self, other):
        if isinstance(other, int):
            if other <= 0:
                raise ValueError(f'Stoichiometry must be positive! We got {other}!')
            return WeightedSpecies(species=self, stoichiometry=other)
        else:
            raise TypeError(f'Species can be multiplied by integer only! We got {type(other)}')

    def __add__(self, other):
        if isinstance(other, Species):
            return [WeightedSpecies(species=other),
                    WeightedSpecies(species=self)]
        else:
            if isinstance(other, WeightedSpecies):
                if other.species == self:
                    return WeightedSpecies(species=other.species, stoichiometry=other.stoichiometry+1)
                else:
                    return [other, WeightedSpecies(species=self)]

        raise NotImplementedError(f'Either Species or WeightedSpecies can be added together! We got {type(other)}')

    @staticmethod
    def flatten_list(in_list) -> List:
        """Helper function to flatten lists
        """
        out_list = []
        if not isinstance(in_list,list):
            out_list.append(in_list)
        else:
            for element in in_list:
                if isinstance(element, list):
                    out_list += Species.flatten_list(element)
                else:
                    out_list += [element]
        return out_list


class WeightedSpecies:
    def __init__(self, species: Species, stoichiometry:int=1):
        """Container object for a all types of species and its stoichiometry
        """
        self.species: Species = species
        self.stoichiometry: int = stoichiometry

    @property
    def stoichiometry(self):
        return self._stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        if not isinstance(new_stoichiometry, int) or new_stoichiometry <= 0:
            raise ValueError(f'Stoichiometry must be positive integer! We got {new_stoichiometry}!')
        self._stoichiometry = new_stoichiometry

    def pretty_print(self, **kwargs):
        return f'{self.stoichiometry if self.stoichiometry > 1 else ""}{self.species.pretty_print(**kwargs)}'

    def replace_species(self, *args, **kwargs):
        return self.species.replace_species(*args, **kwargs)

    @staticmethod
    def _count_weighted_species(weighted_species: List):
        """Helper function merge the same species in a list with different stoichiometry

        >>> s1 = Species(name='a')
        >>> ws1 = WeightedSpecies(species=s1, stoichiometry=2)
        >>> ws2 = WeightedSpecies(species=s1, stoichiometry=5)
        >>> ws_list = [ws1, ws2]
        >>> freq_dict = WeightedSpecies._count_weighted_species(ws_list)
        >>> len(freq_dict)
        1

        :param weighted_species: list of weighted_species
        :return: unique list of weighted_species, i.e. set(weighted_species)
        """
        # convert to set doesn't work because we need only species equality
        unique_species = []
        for w_species in weighted_species:
            if not any(w_species.species == u_s.species for u_s in unique_species):
                unique_species.append(w_species)

        freq_dict = dict(zip(unique_species, [0]*len(unique_species)))
        for w_species in weighted_species:
            for key in freq_dict:
                if key.species == w_species.species:
                    freq_dict[key] += w_species.stoichiometry

        return freq_dict

    def __eq__(self, other):
        if other.__class__ is self.__class__:
            return (other.species, other.stoichiometry) == (self.species, self.stoichiometry)
        return False

    def __hash__(self):
        return hash(self.species)+hash(self.stoichiometry)


class ComplexSpecies(Species):
    """ A special kind of species which is formed as a complex of two or more species.

        Used for attribute inheritance and storing groups of bounds Species. 
        Note taht in a ComplexSpecies, the order of the species list does not matter.
        This means that ComplexSpecies([s1, s2]) = ComplexSpecies([s2, s1]). 
        This is good for modelling order-indpendent binding complexes.
        For a case where species order matters (e.g. polymers) use OrderedComplexSpecies
    """
    def __init__(self, species: List[Union[Species,str]], name: Union[str,None] = None, material_type="complex", attributes=None, initial_concentration = 0, **keywords):
        if len(species) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")

        self.species = []
        for s in Species.flatten_list(species):
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
        """
        Returns a list of species inside the ComplexSpecies
        """
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

    def replace_species(self, species: Species, new_species: Species):
        """Replaces species with new_species in the entire Complex Species.

        Acts recursively on nested ComplexSpecies
        Does not act in place - returns a new ComplexSpecies.
        """
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self.custom_name is True:
            new_name = self.name
        
        return ComplexSpecies(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes)

    def get_species(self, recursive = False):
        """
        Returns all species in the ComplexSpecies. If recursive = True, returns species inside internal ComplexSpecies recursively as well.
        """
        if not recursive:
            species = [self]
        else:
            species = []
            for s in self.species:
                species += s.get_species(recursive = True)

        return species

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        """

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

        new_species = Species.flatten_list(species)
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

    def replace_species(self, species: Species, new_species: Species):
        """
        Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
        """
        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self.custom_name == True:
            new_name = self.name
        
        return OrderedComplexSpecies(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes)

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
        """
        A more powerful printing function.
        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        """

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
