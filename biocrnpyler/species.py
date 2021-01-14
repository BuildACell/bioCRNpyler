#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import warnings
from typing import List, Union
from .compartments import Compartment
from .polymer import OrderedMonomer, OrderedPolymer, MonomerCollection


class Species(OrderedMonomer):

    """ A formal species object for a CRN
     A Species must have a name. They may also have a material_type (such as DNA,
     RNA, Protein), and a list of attributes.
    """

    def __init__(self, name: str, material_type="", attributes: Union[List, None] = None,
                 compartment=None, **keywords):
        OrderedMonomer.__init__(self, **keywords)

        self.name = name
        self.material_type = material_type
        self._attributes = []  # Set this to avoid errors
        self.attributes = attributes
        self.compartment = compartment

    @property
    def attributes(self):
        if(not hasattr(self, "_attributes")):
            self._attributes = []
        return self._attributes

    @attributes.setter
    def attributes(self, attributes):
        if(not hasattr(self, "_attributes")):
            self._attributes = []
        if attributes is not None:
            if not isinstance(attributes, list):
                attributes = list(attributes)
            for attribute in attributes:
                self.add_attribute(attribute)
        elif attributes is None:
            self._attributes = []

    def remove_attribute(self, attribute: str):
        """
        removes an attribute from a Species
        """
        if(not hasattr(self, "_attributes")):
            self._attributes = []
            return
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(
        ), "Attribute: %s must be an alpha-numeric string" % attribute
        if attribute in self.attributes:
            new_attrib = []
            for attrib in self._attributes:
                if(attrib == attribute):
                    pass
                else:
                    new_attrib += [attrib]
            self._attributes = new_attrib

    def add_attribute(self, attribute: str):
        """
        Adds attributes to a Species
        """
        if(not hasattr(self, "_attributes")):
            self._attributes = []
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(
        ), "Attribute: %s must be an alpha-numeric string" % attribute
        if attribute not in self.attributes:
            self._attributes.append(attribute)

    @property
    def name(self):
        if self._name is None:
            return ""
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        if name is None:
            raise TypeError("Name must be a string.")
        else:
            self._name = self._check_name(name)

    @property
    def compartment(self):
        return self._compartment

    @compartment.setter
    def compartment(self, compartment):
        if compartment is None:
            self._compartment = Compartment(name="default")
        else:
            if isinstance(compartment, str):
                self._compartment = Compartment(
                    name=self._check_name(compartment))
            elif isinstance(compartment, Compartment):
                self._compartment = compartment

    # Use OrderedMonomers getter
    direction = property(OrderedMonomer.direction.__get__)

    @direction.setter
    def direction(self, direction):
        """
        This is inheritted from OrderedMonomer.
        A species with direction will use it as an attribute as well.
        This is overwritten to make direction an attribute
        """

        self._direction = direction
        if direction is not None:
            self.add_attribute(direction)

    def remove(self):
        """
        Added functionality to remove direction as an attribute.
        """
        if self.direction is not None:
            self.remove_attribute(self.direction)
        OrderedMonomer.remove(self)  # call the OrderedMonomer function

    # Note: this is used because properties can't be overwritten without setters being overwritten in subclasses.
    def _check_name(self, name):
        """
        Check that the string contains only underscores and alpha-numeric characters or is None. 
        Additionally cannot end in "_" or contain double "__", also cannot start with a number
        """
        if name is None:
            return name
        elif isinstance(name, str):
            no_underscore_string = name.replace("_", "")
            if no_underscore_string.isalnum() and "__" not in name and name[len(name)-1] != "_" and name[0].isalnum():
                return name
            else:
                raise ValueError(
                    f"name attribute {name} must consist of letters, numbers, or underscores and cannot contain double underscores or begin/end with a special character.")
        else:
            raise TypeError("Name must be a string.")

    @property
    def material_type(self):
        return self._material_type

    @material_type.setter
    def material_type(self, material_type: str):
        """ 
        Check that the string contains is alpha-numeric characters or "_" and that the first character is a letter. 
        If the name is a starts with a number, there must be a material type.
        """
        if material_type in [None, ""] and self.name[0].isnumeric():
            raise ValueError(
                f"species name: {self.name} contains a number as the first character and therefore requires a material_type.")
        elif material_type == None:
            self._material_type = None
        elif (material_type.replace("_", "").isalnum() and material_type.replace("_", "")[0].isalpha() and "__" not in material_type and material_type[len(material_type)-1] != "_") or material_type == "":
            self._material_type = material_type
        else:
            raise ValueError(
                f"material_type {material_type} must be alpha-numeric and start with a letter.")

    def __repr__(self):
        txt = ""
        if self.material_type not in ["", None]:
            txt = self.material_type + "_"

        txt += self.name

        if len(self.attributes) > 0 and self.attributes != []:
            for i in self.attributes:
                if i is not None:
                    txt += "_" + str(i)
        if self.compartment.name != 'default':
            # Only add a compartment name if it is not the default one.
            txt += "_" + self.compartment.name
        txt.replace("'", "")
        return txt

    def replace_species(self, species, new_species):
        if not isinstance(species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        if self == species:
            return new_species
        else:
            return self

    def get_species(self, **kwargs):
        """
        Used in some recursive calls where ComplexSpecies returns a list and Species will return just themselves (in a list)
        """
        return [self]

    def pretty_print(self, show_material=True, show_compartment=False, show_attributes=True,
                     show_initial_condition=False, **kwargs):
        """
        #A more powerful printing function. 
        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        """
        txt = ""
        if self.material_type not in ["", None] and show_material:
            txt = self.material_type + "["

        txt += self.name

        if self.compartment not in ["default", None] and show_compartment:
            txt += " in " + self.compartment.name + "."

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")
        if(hasattr(self, "direction") and self.direction is not None):
            txt += "-"+self.direction
        if self.material_type not in ["", None] and show_material:
            txt += "]"

        return txt

    def __eq__(self, other):
        """
        Overrides the default implementation
        Two species are equivalent if they have the same name, type, and attributes
        :param other: Species instance
        :return: boolean
        """

        #Note: "==" equality between parents and children can result in loops, so string equality is used
        if isinstance(other, Species) \
                and self.material_type == other.material_type \
                and self.name == other.name \
                and set(self.attributes) == set(other.attributes)\
                and str(self.parent) == str(other.parent)\
                and self.compartment == other.compartment\
                and self.position == other.position:
            return True
        else:
            return False

    def __gt__(self, Species2):
        return self.name > Species2.name

    def __lt__(self, Species2):
        return self.name < Species2.name

    def __hash__(self):
        return str.__hash__(repr(self))

    def __contains__(self, item):
        return item in self.get_species()

    @staticmethod
    def flatten_list(in_list) -> List:
        """Helper function to flatten lists
        """
        out_list = []
        if not isinstance(in_list, list):
            out_list.append(in_list)
        else:
            for element in in_list:
                if isinstance(element, list):
                    out_list += Species.flatten_list(element)
                elif(element is None):
                    pass
                else:
                    out_list += [element]
        return out_list


class WeightedSpecies:
    def __init__(self, species: Species, stoichiometry: int = 1):
        """Container object for a all types of species and its stoichiometry
        """
        self.species: Species = species
        self.stoichiometry: int = stoichiometry

    @property
    def stoichiometry(self):
        return self._stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        if new_stoichiometry <= 0:
            raise ValueError(
                f'Stoichiometry must be positive integer! We got {new_stoichiometry}!')
        self._stoichiometry = int(new_stoichiometry)

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
    """ 
    !!!ComplexSpecies and OrderedComplexSpecies should ALWAYS be created with the Complex function!!!

        A special kind of species which is formed as a complex of two or more species.
        Used for attribute inheritance and storing groups of bounds Species. 
        Note taht in a ComplexSpecies, the order of the species list does not matter.
        This means that ComplexSpecies([s1, s2]) = ComplexSpecies([s2, s1]). 
        This is good for modelling order-indpendent binding complexes.
        For a case where species order matters (e.g. polymers) use OrderedComplexSpecies
    """

    def __init__(self, species: List[Union[Species, str]], name: Union[str, None] = None, material_type="complex",
                 attributes=None, compartment=None, **keywords):

        # A little check to enforce use of Complex() to create ComplexSpecies
        if "called_from_complex" not in keywords or not keywords["called_from_complex"]:
            warnings.warn(
                "ComplexSpecies should be created using the Complex([List of Species]) function, not directly!")

        # Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")
        self.species = species

        # call super class
        Species.__init__(self, name=name, material_type=material_type, attributes=attributes,
                         compartment=compartment)

    def __repr__(self):
        """
        ComplexSpecies add an additional "_" onto the end of their string representation
        This ensures that some edge cases are differentiated.
        """
        txt = Species.__repr__(self)
        txt += "_"
        return txt

    @property
    def name(self):
        if self._name is None:
            name = ""
            for s in self.species_set:
                count = self.species.count(s)
                name += str(s)+"_"
                if count > 1:
                    name += f"{count}x_"
            name = name[:-1]
            return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    def __contains__(self, item):
        """
        Returns a list of species inside the ComplexSpecies
        """
        if not isinstance(item, Species):
            raise ValueError(
                "Operator 'in' requires chemical_reaction_network.Species (or a subclass). Received: "+str(item))
        if item in self.species:
            # this is the base case
            return True
        else:
            # this is the recursive part. We want to check all
            # internal complexes for the thing we're looking for
            for content in self.species:
                if isinstance(content, ComplexSpecies):
                    if item in content:
                        return True
            # if we got here then we've failed to find it
            return False

    @property
    def species_set(self):
        species_set = list(set(self.species))
        list.sort(species_set, key=lambda s: repr(s))
        return species_set

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if not isinstance(species, list):
            raise TypeError(f"species must be a list: recieved {species}.")
        species = Species.flatten_list(species)
        if not all(isinstance(s, Species) for s in species):
            raise TypeError(
                f"recieved a non-species as a member of the list species: {species}.")
        else:
            list.sort(species, key=lambda s: repr(s))
            self._species = species

    def replace_species(self, species: Species, new_species: Species):
        """
        Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
        Does not act in place - returns a new ComplexSpecies. 
        """
        if not isinstance(species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self._name is not None:
            new_name = self.name

        return Complex(species=new_species_list, name=new_name, material_type=self.material_type, attributes=self.attributes)

    def get_species(self, recursive=False):
        """
        Returns all species in the ComplexSpecies. If recursive = True, returns species inside internal ComplexSpecies recursively as well.
        """
        if not recursive:
            species = [self]
        else:
            species = [self]
            for s in self.species:
                species += s.get_species(recursive=True)

        return species

    def pretty_print(self, show_material=True, show_compartment=False, show_attributes=True,
                     show_initial_condition=False, **kwargs):
        """
        A more powerful printing function. 
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
            txt += s.pretty_print(show_material=show_material,
                                  show_attributes=False)+":"
        txt = txt[:-1]

        if self.compartment not in ["default", None] and show_compartment:
            txt += " in " + self.compartment.name + "."

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")
        if(hasattr(self, "direction") and self.direction is not None):
            txt += "-"+self.direction
        txt += "]"

        return txt


class OrderedComplexSpecies(ComplexSpecies):
    """ 
    !!!ComplexSpecies and OrderedComplexSpecies should ALWAYS be created with the Complex function!!!

    A special kind of species which is formed as a complex of two or more species.
    In OrderedComplexSpecies the order in which the complex subspecies are is defined
    denote different species, eg [s1, s2, s3] != [s1, s3, s2].
    Used for attribute inheritance and storing groups of bounds Species. 
    """

    def __init__(self, species, name=None, material_type="ordered_complex", attributes=None,
                 compartment=None, **keywords):
        # Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")
        self.species = species

        # Call the Species superclass constructor
        Species.__init__(self, name=name, material_type=material_type, attributes=attributes,
                         compartment=compartment)

    @property
    def name(self):
        if self._name is None:
            name = ""
            for s in self.species:
                if isinstance(s, str):
                    s = Species(name=s)
                if s.material_type not in [""]:
                    name += f"{s.material_type}_{s.name}_"
                else:
                    name += f"{s.name}_"
            name = name[:-1]
            return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if not isinstance(species, list):
            raise TypeError(f"species must be a list: recieved {species}.")
        species = Species.flatten_list(species)
        if not all(isinstance(s, Species) for s in species):
            raise TypeError(
                f"recieved a non-species as a member of the list species: {species}.")
        else:
            self._species = species

    def replace_species(self, species: Species, new_species: Species):
        """
        Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
        """
        if not isinstance(species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError(
                'species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self._name is not None:
            new_name = self.name

        return Complex(species=new_species_list, name=new_name, material_type=self.material_type, attributes=self.attributes, ordered=True)

    def pretty_print(self, show_material=True, show_compartment=False, show_attributes=True,
                     show_initial_condition=False, **kwargs):
        """
        A more powerful printing function. 
        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        """

        txt = ""
        if self.material_type not in ["", None] and show_material:
            txt += self.material_type

        txt += "["

        for s in self.species:
            txt += s.pretty_print(show_material=show_material,
                                  show_attributes=False)+":"
        txt = txt[:-1]

        if self.compartment not in ["default", None] and show_compartment:
            txt += " in " + self.compartment.name + "."

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")
        txt += "]"

        return txt


class OrderedPolymerSpecies(OrderedComplexSpecies, OrderedPolymer):
    """
    A class to represent OrderedPolymers which can also participate in chemical reactions.
    OrderedPolymerSpecies is made up of Species (which are also OrderedMonomers).

    The Species inside an OrderedPolymerSpecies are meant to model multiple binding sites and/or
    functional regions. ComplexSpecies can be formed inside an OrderedPolymer by passing
    the internal Species at a specific location.

    When used as an input to a reaction, OrderedPolymerSpecies can be passed 
    or one if its internal Species (eg a Species with Species.parent = OrderedPolymerSpecies)
    can also be used to produce the same reaction. This allows flexibility in the arguments to
    different Mechanisms. Sometimes, it is convenient to pass in the OrderedPolymerSpecies,
    sometimes it is convenient to pass an internal Species. Both will work from the point of view
    of any Mechanism.
    """
    default_material = "ordered_polymer"

    def __init__(self, species, name=None, base_species=None, material_type=default_material,
                 compartment=None, attributes=None, circular=False):

        self.material_type = material_type
        self.compartment = compartment
        self.parent = None
        self.position = None
        self.direction = None

        self.circular = circular

        if attributes is None:
            self.attributes = []
        else:
            self.attributes = attributes

        self._name = OrderedComplexSpecies._check_name(self, name)

        self.material_type = material_type
        #self.species = []
        monomers = []
        for specie in species:
            if isinstance(specie, Species) and isinstance(specie, OrderedMonomer):
                monomers += [specie]
            elif (isinstance(specie, tuple) or isinstance(specie, list)) and (isinstance(specie[0], Species) and isinstance(specie[0], OrderedMonomer)):
                monomers += [specie]
            elif isinstance(specie, OrderedPolymer):
                raise NotImplementedError(
                    f"OrderedPolymer cannot be used as a Monomer at this time.")
            else:
                raise ValueError(
                    "{} should be a Species or list [Species, 'direction']".format(specie))
                # only species are acceptable

        OrderedPolymer.__init__(self, monomers)
        self.material_type = material_type

        if(base_species is None):
            self.base_species = Species("NA", material_type=material_type)
            self.base_species._name = self.name #Bipass check until base-species is removed in Andrey's PR
        elif(isinstance(base_species, Species)):
            self.base_species = base_species
        else:
            raise TypeError("base_species is of type "+type(base_species) +
                            " which is not acceptable. Use Species or str")

    @classmethod
    def from_polymer_species(cls, ops, replace_dict, **keywords):
        """
        Created a new OrderedPolymerSpecies with certain monomers replaced based upon replace_dict:

        inputs: replace_dict {monomer index --> new Species}
        outputs: OrderedPolymerSpecies
        """
        monomers = []
        for i in range(len(ops.polymer)):
            if i in replace_dict:
                monomers.append(replace_dict[i])
            else:
                monomer = copy.copy(ops[i])
                direction = monomer.direction
                monomer.remove()
                monomers.append([monomer, direction])

        #Set keywords
        if "circular" not in keywords:
            keywords["circular"] = ops.circular
        if "material_type" not in keywords:
            keywords["material_type"] = ops.default_material
        if "compartment" not in keywords:
            keywords["compartment"] = ops.compartment
        if "attributes" not in keywords:
            keywords["attributes"] = ops.attributes

        return cls(monomers, **keywords) #Produces a new OrderedPolymerSpecies

    @property
    def species_set(self):
        return set(self.polymer)

    @property
    def species(self):
        return self.polymer

    def get_species_list(self):
        return self.polymer

    @property
    def circular(self):
        if("circular" in self.attributes):
            return True
        else:
            return False

    @circular.setter
    def circular(self, value: bool):
        if(value):
            self.add_attribute("circular")
        else:
            self.remove_attribute("circular")

    def set_species_list(self, spec_tuple: tuple):
        OrderedPolymer.__init__(self, spec_tuple)

    @property
    def name(self):
        if self._name is None:
            name = ""
        else:
            name = self._name
        outlst = []

        for monomer in self.polymer:
            assert(isinstance(monomer, Species))
            pname = monomer.name
            pdir = None
            if(hasattr(monomer, "direction")):
                pdir = monomer.direction
            if(pdir is not None):
                outlst += [pname+"_"+str(pdir)[0]]
            else:
                outlst += [pname]
        if(self.circular):
            outlst += ["o"]
        name = '_'.join(outlst)
        return name

    def __hash__(self):
        ophash = OrderedPolymer.__hash__(self)
        ophash += hash(self.circular)+hash(self.base_species) + \
            hash(self.name)+hash(self.material_type)
        return ophash

    def replace(self, position, part, direction=None):
        # TODO only change the name if the part we are replacing is actually different
        mydir = direction
        if((mydir is None) and (part.direction is not None)):
            mydir = part.direction
        if(part == self.polymer[position] and self.polymer[position].direction == mydir):
            # in this case we are replacing a part with the same thing, so do nothing
            # but it could be true that the reference changes? That shouldnt be
            pass
        else:
            OrderedPolymer.replace(
                self, position=position, part=part, direction=mydir)

    def __contains__(self, item):
        for part in self.species:
            part = copy.copy(part) #Only compare things which aren't None
            if item.parent is None:
                part.parent = None
            if item.direction is None:
                part.direciton = None
            if item.position is None:
                part.position = None
            if part == item:
                return True
        return False


class PolymerConformation(Species, MonomerCollection):
    """
    This class stores a set of PolymerSpecies and a set of connections between them in the form of ComplexSpecies containing Monomers inside the PolymerSpecies.

    The main function of this class is to provide a unique name to each conformation. The name is given by:

        conformation__[PolymerSpecies 1]_..._[PolymerSpecies N]_[ComplexSpecies_1 parent Polymer indices]_[ComplexSpecies_1]..._[ComplexSpecies_M]__
    
    where the list of PolymerSpecies and ComplexSpecies are in alphabetical order. 
    The ComplexSpecies parent Polymer indices notes which Polymers each Species in the ComplexSpecies comes from, with 'n' used for None.

    In general, users should not produce PolymerConformations directly. The Complex function will automatically produce these
    when a complex is formed involving Multiple OrderedMonomers contained within one or more PolymerSpecies.

    In effect, this can be thought of as a data structure for a hypergraph. The monomers of the PolymerSpecies are
    vertices and ComplexSpecies form edges that connect an arbitrary number of vertices (potentially including 
    other Species as well). Note that this class allows for multiple edges between the same sets of vertices.
    """

    def __init__(self, complexes, material_type = "conformation", name = None, **keywords):
        """
        complexes: a list of ComplexSpecies each of which must contain Monomers from the OrderedPolymerSpecies in the conformation
        """
        Species.__init__(self, name = name, material_type = material_type, **keywords)

        self.complexes = complexes

    @classmethod
    def from_polymer_conformation(cls, pcs, complexes, **keywords):
        """
        This function produces a new PolymerConformation from previously existing PolymerConformations and new Complexes.

        pcs: a list of PolymerConformations
        complexes: a list of complexes to add to the polymer conformation
        """

        if not isinstance(pcs, list) or not any([isinstance(pc, PolymerConformation) for pc in pcs]):
            raise TypeError(f"pcs must be a list of PolymerConformations. Recieved {pcs}.")

        #generate a list of all complexes
        for pc in pcs:
            complexes += pc.complexes

        return cls(complexes, **keywords)

    @classmethod
    def from_polymer_replacement(cls, pc, old_polymers, new_polymers, **keywords):
        """
        This function produces a PolymerConformation from a previously existing PolymerConformation by replacing old_polymers with new_polymers
        
        pc: the PolymerConformation to replace polymers from.
        old_polymers: a list of PolymerSpecies instances. These must be the same instances stored inside pc or an error is thrown.
        new_polymers: a list of new PolymerSpecies instances to replace each of the old_polymers. Must be the same length as old_polymers.
        """
        if not isinstance(pc, PolymerConformation):
            raise TypeError(f"pc must be a PolymerConformation. Recieved {pc}.")

        if not isinstance(new_polymers, list) and all([isinstance(p, OrderedPolymerSpecies) for p in new_polymers]):
            raise TypeError(f"new_polymers must be a list of OrderedPolymerSpecies. Recieved: {new_polymers}")

        if not isinstance(old_polymers, list) and all([isinstance(p, OrderedPolymerSpecies) for p in old_polymers]):
            raise TypeError(f"new_polymers must be a list of OrderedPolymerSpecies. Recieved: {old_polymers}")

        if not all([any([p is pp for pp in pc.polymers]) for p in old_polymers]):
            raise ValueError("All OrderedPolymerSpecies in old_polymers must be contained (as instances, not string equivalents) in pc.polymers.")

        if len(old_polymers) != len(new_polymers):
            raise ValueError("old_polymers and new_polymers must be the same length.")

        

        complexes = pc.complexes
        new_complexes = []

        #Create a set of new ComplexeSpecies with the correct monomers inside
        for c in complexes:
            pol_inds = pc.get_polymer_indices(c)
            species = []
            for ci, s in enumerate(c.species):
                pi = pol_inds[ci]
                #If s is part of a polymer:
                if pi is not None:
                    #add the new_polymers[pi][s.position] to the new species
                    species.append(new_polymers[pi][s.position])
                #Otherwise, the species is not part of a polymer
                else:
                    species.append(s)

            #Create the new Complex
            if isinstance(c, OrderedComplexSpecies): #is the complex ordered?
                new_complexes += [OrderedComplexSpecies(species, material_type = c.material_type, attributes = c.attributes, comparment = c.compartment)]
            elif isinstance(c, ComplexSpecies):
                new_complexes += [ComplexSpecies(species, material_type = c.material_type, attributes = c.attributes, comparment = c.compartment)]
            else:
                raise TypeError(f"Invalid object found in PolymerConformation.complexes: {c}.")

        #Set keywords
        if "material_type" not in keywords:
            keywords["material_type"] = pc.material_type
        if "compartment" not in keywords:
            keywords["compartment"] = pc.compartment
        if "attributes" not in keywords:
            keywords["attributes"] = pc.attributes

        return cls(new_complexes, **keywords)

    #To be a valid MonomerCollection
    @property
    def monomers(self):
        return tuple(self._species)

    @property
    def polymers(self):
        return self._polymers

    @property
    def complexes(self):
        return self._complexes


    @complexes.setter
    def complexes(self, complexes):
        """
        This setter copies the complexes and the polymers they connect into the PolymerConformation.
        This is done in such a way to preserve references between parents and children without relying on hash functions (in case two polymers are identical).
        """
        self._complexes = []
        self._polymers = []
        if not isinstance(complexes, list) or len(complexes) < 1 or not all([isinstance(c, ComplexSpecies) for c in complexes]):
            raise TypeError(f"complexes must be a list containing ComplexSpecies. Recieved {complexes}.")

        if any([sum([cc==c and cc.species == c.species for cc in complexes]) > 1 for c in complexes]):
            raise ValueError("Complexes contains two identical complexes.")

        #Sort the complexes by their name, and the ids of the polymers inside of them to get a unique ordering for identically named Complexes and Polymers
        #this will produce a unique ordering of the internal polymers as well
        sort_func = lambda c: (c, tuple(id(s.parent) for s in c.species))
        complexes.sort(key = sort_func)


        polymers = [] #Original polymers stored here
        copied_polymers = [] #Polymers are copied here
        copied_complexes = [] #Complexes are copied here

        #Find all the polymers in the complexes. "is" is used instead of equality to deal with the possibility of
        #multiple instances (copies) of the same polymer being bound together.
        for c in complexes:
            c_copy = copy.deepcopy(c)
            complex_contains_polymer = False
            for i, s in enumerate(c.species):
                parent = s.parent

                if parent is not None and isinstance(parent, OrderedPolymerSpecies):
                    complex_contains_polymer = True
                    #parent has not been copied
                    if not any([parent is p for p in polymers]):
                        polymers.append(parent)
                        parent_copy = copy.deepcopy(parent)
                        parent_copy.parent = self #set the polymer's parent
                        copied_polymers.append(parent_copy)
                        c_copy.species[i] = parent_copy[s.position] #set s to be the new version from the deep copy
                    #parent has already been copied
                    else:
                        index = [parent is p for p in polymers].index(True) #get parent index
                        parent_copy = copied_polymers[index]
                        c_copy.species[i] = parent_copy[s.position] #set s to be the new version from the deep copy

                #This edgecase should not occur.
                elif parent is not None and not isinstance(parent, OrderedPolymerSpecies):
                    raise ValueError(f"Species {s} found inside complex {c} with a parent {parent} which is not an OrderedPolymerSpecies.")

            if not complex_contains_polymer:
                raise ValueError(f"Complex {c} does not contain any Species inside of Polymers.")

            #Set the parent of the ComplexSpecies to self
            c_copy.parent = self
            copied_complexes.append(c_copy)

        
        #Sort the polymers by their name
        #polymer_sort_func = lambda p: (str(p), id(p))
        copied_polymers.sort()
        self._polymers = copied_polymers

        #Sort the complexes ordered by their polymer indices then name
        complex_sort_func = lambda c: (tuple(self.get_polymer_indices(c)), c)
        copied_complexes.sort(key = complex_sort_func) 
        self._complexes = copied_complexes
    
    def get_polymer_indices(self, c):
        #Takes a complex and returns a list of indices to the polymers that complex contains monomers from
        #this complex should be in self.complexes. 
        #If a Species in the Complex is not in a Polymer, None is added to the list.
        indices = []
        for j, s in enumerate(c.species):
            parent = s.parent
            for i, p in enumerate(self.polymers):
                if parent is p:
                    indices.append(i)
            if len(indices) <= j:
                indices.append(None) 

        return indices

    def get_polymer_positions(self, c, polymer_ind):
        #Takes a complex and the index of a polymer in the conformation and returns a list of positions that ComplexSpecies is bound at
        p = self.polymers[polymer_ind]
        positions = []
        for s in c.species():
            if s.parent is p:
                positions.append(s.position)
        return positions

    def get_polymer(self, p):
        #Takes a polymer and returns a matching instance inside self.polymers (or None)
        p = copy.copy(p)
        p.parent = self
        if p in self.polymers:
            i = self.polymers.index(p)
            return self.polymers[i]
        else:
            return None

    def get_complex(self, c):
        #takes a ComplexSpecies and returns the instance of this Species in the PolymerConformation (or None)
        c = copy.copy(c)
        c.parent = self
        if c in self.complexes:
            i = self.complexes.index(c)
            return self.complexes[i]
        else:
            return None

    def remove_complex(self, c):
        #Removes a copy of this complex from PolymerConformation if possible
        if c in self.complexes:
            self.complexes.remove(c)
        else:
            raise ValueError(f"Complex {c} not in PolymerConformation {self}.")

    @property
    def name(self):
        if self._name is None:
            name = ""
            for p in self.polymers:
                name += "_"+str(p)
            for c in self.complexes:
                parent_inds = self.get_polymer_indices(c)
                name += "_"
                for ind in parent_inds:
                    if ind is None:
                        name+="n"
                    else:
                        name+=f"p{ind}"
                name += "_"+str(c)
            return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    def __repr__(self):
        """
        PolymerConformations add an additional "_" onto the end of their string representation
        This ensures that some edge cases are differentiated.
        """
        txt = Species.__repr__(self)
        txt += "_"
        return txt


class Complex:
    """
     !!!ComplexSpecies and OrderedComplexSpecies should ALWAYS be created with the Complex function!!!

    Complex is not a class that gets instantiated - it creates ComplexSpecies and OrderedComplexSpecies.
    The Logic encoded in the __new__ function is used to insert these classes into the binding sites of
    OrderedPolymerSpecies.

    arguments:
    species: a list of species to put into ComplexSpecies or OrderedComplexSpecies

    keywords:
    ordered: whether to produce an OrderedComplexSpecies (default = False)
    """
    def __new__(cls, *args, **keywords):
        """
        This function effectively produces the instance of the correct Species Class based upon the arguments passed in.

        Cases: Here species refer to the Species in the Species list passed into the construct.
        1.  No Species have parents. 
            Produces: an ComplexSpecies or an OrderedComplexSepcies
        2.  A single Species S has a parent which is an OrderedPolymerSpecies with no parent. 
            Produces: an OrderedPolymerSpecies with a ComplexSpecies or OrderedComplexSpecies containing S in S's location in the OrderedPolymerSpecies.
        3.  A single Species S has a parent which is an OrderedPolymerSpecies with a PolymerConformation as a parent.
            Produces: a (Ordered)ComplexSpecies in a PolymerConformation with in the OrderedPolymerSpecies inside the PolymerConformation
        4.  Multiple Species S have parents which are OrderedPolymerSpecies without parents.
            Produces: a (Ordered)ComplexSpecies containing all S inside a PolymerConformation
        5.  One or More Species S have parents which are OrderedPolymerSpecies with parents and/or PolymerConformations
            Produces: a (Ordered)ComplexSpecies containing all S inside a PolymerConformation which merges all PolymerComformation Complexes.
        """
        species = []
        # below is extracting the "species" keywords from the args
        keywarg = None
        if("species" in keywords):
            keywarg = True
            species = keywords["species"]
            del keywords["species"]
        elif(len(args) >= 1):
            keywarg = False
            species = args[0]
            args = args[1:]

        if not isinstance(species, list):
            raise TypeError(f"First argument to Complex (or species keyword), must be a list of Species; recieved {species}.")

        # Check whether ot make a ComplexSpecies or OrderedComplexSpecies
        if "ordered" in keywords:
            ordered = keywords["ordered"]
            del keywords["ordered"]
            ComplexClass = OrderedComplexSpecies
        else:
            ordered = False
            ComplexClass = ComplexSpecies

        #Use to supress errors in ComplexSpecies and OrderedComplexSpecies
        keywords["called_from_complex"] = True

        parent_species = []  # parent_species is a list of OrderedPolymerSpecies and/or PolymerConformations
        bindlocs = []  # bindloc is the location a Species is bound to parent_species
        insertlocs = [] #insertloc is the location of the species inside the Complex.species list (important to maintain order in OrderedComplexSpecies)
        child_species = [] #the species with parents
        other_species = []  # Other species in the Complex

        # Below cycle through species and see if one has a parent. If it does, that means the species is
        # in an OrderedPolymerSpecies and the Complex should be formed around it.
        for i, specie in enumerate(species):
            if(hasattr(specie, "parent") and (specie.parent is not None)):
                #This adds to a growing list of parents, which will be placed in an PolymerConformation
                # It is very important to deepcopy here because the underlying OrderedPolymerSpecies or PolymerConformation will be modified.
                parent_species.append(specie.parent)
                bindlocs.append(specie.position)
                insertlocs.append(i)
                child_species.append(specie)
            else:
                other_species += [specie]

        #Case 1: If no OrderedPolymerSpecies is found, just call the regular constructor.
        if len(parent_species) == 0:
            return ComplexClass(species, *args, **keywords)

        #Case 2 - 3: the Complex is being formed inside an OrderedPolymerSpecies (and only a single species has a parent).
        elif len(parent_species) == 1 and isinstance(parent_species[0], OrderedPolymerSpecies):
            parent_species = parent_species[0]
            bindloc = bindlocs[0]

            #Ensure the monomer in the Complex has no parent
            #new_species = other_species
            #monomer = copy.copy(child_species[0]) #Copy the monomer
            #monomer.remove()
            #new_species.insert(insertlocs[0], monomer) #Insert the monomer back into the Species

            # Create an OrderedcomplexSepcies or ComplexSpecies
            new_complex = ComplexClass(species, *args, **keywords)

            #Create a new OrderedPolymerSpecied which is copied from the parent with the new complex replacing bindloc (inheriting the same direction).
            new_polymer_species = OrderedPolymerSpecies.from_polymer_species(parent_species, {bindloc:[new_complex, child_species[0].direction]})

            #Case 2: OrderedPolymerSpecies has no parent
            if parent_species.parent is None:
                return new_polymer_species[bindloc]

            #Case 3: OrderedPolymerSpecies is inside a PolymerConformation
            elif isinstance(parent_species.parent, PolymerConformation):

                #create a new PolymerConformation that replaces the appropriate monomer in the polymer
                new_pc = PolymerConformation.from_polymer_replacement(parent_species.parent, [parent_species], [new_polymer_species])
                return new_pc.get_polymer(new_polymer_species)[bindloc] #return the newly created Monomer attached to the new OrderedPolymerSpecies inside the new PolymerConformation
            else:
                #This error should never occur
                raise TypeError(f"Unknown parent type {type(parent_species)} recieved for {parent_species} .parent {parent_species.parent}.")

        #Case 4-5: In the following cases, multiple species may have parents
        else:
            #Case 4: Multiple OrderedPolymerSpecies are being Complexed Together (and they aren't already part of PolymerConformations)
            if all([isinstance(p, OrderedPolymerSpecies) for p in parent_species]) and all([p.parent is None for p in parent_species]):

                #Create a Complex and PolymerConformation
                new_complex = ComplexClass(species, *args, **keywords)
                new_pc = PolymerConformation([new_complex])
                return new_pc.get_complex(new_complex)

            #Case 5: Multiple species in one more more PolymerConformations are being Complexed Together
            else:
                pcs = []
                merged_species = other_species #these Species will go inside the ComplexSpecies later

                #Cycle through Species with parents
                insert_loc_offset = 0
                for i, p in enumerate(parent_species):
                    s = child_species[i]

                    #if the parent is a OrderedPolymerSpecies
                    if isinstance(p, OrderedPolymerSpecies):
                        merged_species.insert(insertlocs[i]+insert_loc_offset, s)

                        #if the Polymer is already in a Conformation...
                        if p.parent is not None and not any([p.parent is pp for pp in pcs]):
                            pcs.append(p.parent)

                    #if the parent is a PolymerConformation and child is a ComplexSpecies
                    elif isinstance(p, PolymerConformation) and isinstance(s, ComplexSpecies) and not isinstance(s, OrderedPolymerSpecies):
                        #Store all the unique PolymerConformations
                        if not any([p is pp for pp in pcs]):
                            p = copy.copy(p)
                            p.remove_complex(s)
                            pcs.append(p)

                        #Merge the species lists
                        merged_species = merged_species[:insertlocs[i]+insert_loc_offset]+s.species+merged_species[insertlocs[i]+insert_loc_offset:]
                        insert_loc_offset += len(s.species)-1 #this takes care of ordering offsets during the merge
                    else:
                        raise TypeError(f"Cannot form a complex from {species}. Invalid Parent Species {p} for child {s}.")

                # Create a Complex and merged PolymerConformation
                new_complex = ComplexClass(merged_species, *args, **keywords)
                new_pc = PolymerConformation.from_polymer_conformation(pcs, [new_complex])
                return new_pc.get_complex(new_complex)

