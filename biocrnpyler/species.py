#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .sbmlutil import *
import warnings
import numpy as np
from typing import List, Union, Dict
import itertools
import copy
from .polymer import OrderedPolymer, OrderedMonomer


class Species(OrderedMonomer):

    """ A formal species object for a CRN
     A Species must have a name. They may also have a material_type (such as DNA,
     RNA, Protein), and a list of attributes.
    """

    def __init__(self, name: str, material_type="", attributes: Union[List,None] = None,
                 initial_concentration=0, **keywords):
        OrderedMonomer.__init__(self,**keywords)

        self.name = name
        self.material_type = material_type
        self.initial_concentration = initial_concentration
        self._attributes = [] #Set this to avoid errors
        self.attributes = attributes

    @property
    def attributes(self):
        if(not hasattr(self,"_attributes")):
            self._attributes = []
        return self._attributes

    @attributes.setter
    def attributes(self, attributes):
        if(not hasattr(self,"_attributes")):
            self._attributes = []
        if attributes is not None:
            if not isinstance(attributes,list):
                attributes = list(attributes)
            for attribute in attributes:
                self.add_attribute(attribute)
        elif attributes is None:
            self._attributes = []
    def remove_attribute(self,attribute:str):
        """
        removes an attribute from a Species
        """
        if(not hasattr(self,"_attributes")):
            self._attributes = []
            return
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(), "Attribute: %s must be an alpha-numeric string" % attribute
        if attribute in self.attributes:
            new_attrib = []
            for attrib in self._attributes:
                if(attrib==attribute):
                    pass
                else:
                    new_attrib += [attrib]
            self._attributes = new_attrib

    def add_attribute(self, attribute: str):
        """
        Adds attributes to a Species
        """
        if(not hasattr(self,"_attributes")):
            self._attributes = []
        assert isinstance(attribute, str) and attribute is not None and attribute.isalnum(), "Attribute: %s must be an alpha-numeric string" % attribute
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

    
    #Use OrderedMonomers getter
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
            self.attributes.remove(self.direction)
        OrderedMonomer.remove(self) #call the OrderedMonomer function

    #Note: this is used because properties can't be overwritten without setters being overwritten in subclasses.
    def _check_name(self, name):
        """
        Check that the string contains only underscores and alpha-numeric characters or is None
        """
        if name is None:
            return name
        elif isinstance(name, str):
            no_underscore_string = name.replace("_", "")
            if no_underscore_string.isalnum():
                return name
            else:
                raise ValueError(f"name attribute {name} must consist of letters, numbers, or underscores.")
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
            raise ValueError(f"species name: {self.name} contains a number as the first character and therefore requires a material_type.")
        elif material_type == None:
            self._material_type = None
        elif (material_type.replace("_", "").isalnum() and material_type.replace("_", "")[0].isalpha()) or material_type == "":
            self._material_type = material_type
        else:
            raise ValueError(f"material_type {material_type} must be alpha-numeric and start with a letter.")
    
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
        """
        Used in some recursive calls where ComplexSpecies returns a list and Species will return just themselves (in a list)
        """
        return [self]

    
    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
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

        if len(self.attributes) > 0 and self.attributes != [] and show_attributes:
            txt += "("
            for i in self.attributes:
                if i is not None:
                    txt += str(i)+", "
            txt = txt[:-2]+")"

        txt.replace("'", "")
        if(hasattr(self,"direction") and self.direction is not None):
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

        if isinstance(other, Species) \
                            and self.material_type == other.material_type \
                            and self.name == other.name \
                            and set(self.attributes) == set(other.attributes)\
                            and self.parent == other.parent\
                            and self.position == other.position:
            return True
        else:
            return False
    def __gt__(self,Species2):
        return self.name > Species2.name
    def __lt__(self,Species2):
        return self.name < Species2.name

    def __hash__(self):
        return str.__hash__(repr(self))

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
        if new_stoichiometry <= 0:
            raise ValueError(f'Stoichiometry must be positive integer! We got {new_stoichiometry}!')
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
    def __new__(cls,*args,**keywords):
        species = []
        #below is extracting the "species" keywords from the args
        keywarg = None
        if("species" in keywords):
            keywarg = True
            species = keywords["species"]
            del keywords["species"]
        elif(len(args)>=1):
            keywarg = False
            species = args[0]
            args = args[1:]

        #Check whether ot make a ComplexSpecies or OrderedComplexSpecies
        if "ordered" in keywords:
            ordered = keywords["ordered"]
            del keywords["ordered"]
        else:
            ordered = False

        valent_complex = None #valent_complex is an OrderedPolymer Species
        bindloc = None #bindloc is the location a Species is bound to valent_complex
        other_species = [] #Other species in the Complex

        #Below cycle through species and see if one has a parent. If it does, that means the species is
        #in an OrderedPolymerSpecies and the Complex should be formed around it.
        for specie in species:
            if(hasattr(specie,"parent") and (specie.parent is not None)):
                if(valent_complex is None):
                    #It is very important to deepcopy here because the underlying OrderedPolymerSpecies will be modified.
                    valent_complex = copy.deepcopy(specie.parent)
                    bindloc = specie.position
                else:
                    #If valent_complex has already been found - it means there are two OrderedPolymer
                    #or two Species in the same OrderedPolymer. This kind of binding has not been implemented.
                    raise NotImplementedError("binding together two OrderedPolymerSpecies or two species in the same OrderedPolymerSpecies!")
            else:
                other_species += [specie]

        if len(other_species) == 0:
            raise ValueError("Trying to create a Complex from a single species!")

        #If no OrderedPolymerSpecies is found, just call the regular constructor.
        if(valent_complex is None):
            if ordered:
                #this creates an OrderedComplexSpecies
                #pass in all the args and keywords appropriately
                keywords["called_from_complex"] = True
                return OrderedComplexSpecies(species,*args,**keywords)
            else:
                #this creates aComplexSpecies
                #pass in all the args and keywords appropriately
                keywords["called_from_complex"] = True
                return ComplexSpecies(species,*args,**keywords)

        #This means the Complex is being formed inside an OrderedPolymerSpecies
        else:
            #this is the species around which the complex is being formed
            prev_species = valent_complex[bindloc]
            prev_direction = copy.deepcopy(valent_complex[bindloc].direction)

            #combine what was in the OrderedMonomer with the new stuff in the list
            new_species = other_species+[prev_species]

            #Create an OrderedcomplexSepcies
            if ordered:
                keywords["called_from_complex"] = True
                new_complex = OrderedComplexSpecies(new_species,*args,**keywords)
            #Create a ComplexSpecies
            else:
                keywords["called_from_complex"] = True
                new_complex = ComplexSpecies(new_species,*args,**keywords)
            #now we replace the monomer inside the parent polymer
            valent_complex.replace(bindloc,new_complex,prev_direction)

            return valent_complex[bindloc]

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
    def __init__(self, species: List[Union[Species,str]], name: Union[str,None] = None, material_type = "complex", attributes = None, initial_concentration = 0, **keywords):
        
        #A little check to enforce use of Complex() to create ComplexSpecies
        if "called_from_complex" not in keywords or not keywords["called_from_complex"]:
            warnings.warn("ComplexSpecies should be created using the Complex([List of Species]) function, not directly!")

        #Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")
        self.species = species 

        #call super class
        Species.__init__(self, name = name, material_type = material_type, attributes = attributes, initial_concentration = initial_concentration)

    @property
    def name(self):
        if self._name is None:
            name = ""
            for s in self.species_set:
                count = self.species.count(s)
                if count > 1:
                    name+=f"{count}x_"
                if not (isinstance(s, ComplexSpecies) or s.material_type == ""):
                    name+=f"{s.material_type}_{s.name}_"
                else:
                    name+=f"{s.name}_"
            name = name[:-1]
            return name
        else:
            return self._name
    
    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

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

    @property
    def species_set(self):
        species_set = list(set(self.species))
        list.sort(species_set, key = lambda s:repr(s))
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
             raise TypeError(f"recieved a non-species as a member of the list species: {species}.")
        else:
            list.sort(species, key = lambda s:repr(s))
            self._species = species

    def replace_species(self, species: Species, new_species: Species):
        """
        Replaces species with new_species in the entire Complex Species. Acts recursively on nested ComplexSpecies
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
        if self._name is not None:
            new_name = self.name
        
        return Complex(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes)

    
    def get_species(self, recursive = False):
        """
        Returns all species in the ComplexSpecies. If recursive = True, returns species inside internal ComplexSpecies recursively as well.
        """
        if not recursive:
            species = [self]
        else:
            species = [self]
            for s in self.species:
                species += s.get_species(recursive = True)

        return species


    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
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
    def __init__(self, species, multiplicity, name = None, material_type = "complex", attributes = None, initial_concentration = 0, **keywords):

        if "called_from_complex" not in keywords or not keywords["called_from_complex"]:
            warnings.warn("OrderedComplexSpecies should be created from the Complex([List of Species], ordered = True) function, not directly!")

        if isinstance(species, str):
            species = [Species(name = species)]
        elif not isinstance(species, Species):
            raise ValueError("Multimer must be defined by a Species (or subclasses thereof) and a multiplicity (int).")
        else:
            species = [species]

        ComplexSpecies.__init__(self, species = species*multiplicity, name = name, material_type = material_type, attributes = attributes, initial_concentration = initial_concentration, **keywords)   

class OrderedComplexSpecies(ComplexSpecies):
    """ 
    !!!ComplexSpecies and OrderedComplexSpecies should ALWAYS be created with the Complex function!!!
    
    A special kind of species which is formed as a complex of two or more species.
    In OrderedComplexSpecies the order in which the complex subspecies are is defined
    denote different species, eg [s1, s2, s3] != [s1, s3, s2].
    Used for attribute inheritance and storing groups of bounds Species. 
    """

    def __init__(self, species, name = None, material_type = "ordered_complex", attributes = None, initial_concentration = 0, **keywords):
        #Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError("chemical_reaction_network.complex requires 2 "
                             "or more species in its constructor.")
        self.species = species

        #Call the Species superclass constructor
        Species.__init__(self, name = name, material_type = material_type, attributes = attributes, initial_concentration = initial_concentration)


    
    @property
    def name(self):
        if self._name is None:
            name = ""
            for s in self.species:
                if isinstance(s, str):
                    s = Species(name = s)
                if s.material_type not in ["complex", "ordered_complex", ""]:
                    name+=f"{s.material_type}_{s.name}_"
                else:
                    name+=f"{s.name}_"
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
             raise TypeError(f"recieved a non-species as a member of the list species: {species}.")
        else:
            self._species = species

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
        if self._name is not None:
            new_name = self.name
        
        return Complex(species = new_species_list, name = new_name, material_type = self.material_type, attributes = self.attributes, ordered = True)

    def pretty_print(self, show_material = True, show_attributes = True, **kwargs):
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

class OrderedPolymerSpecies(OrderedComplexSpecies,OrderedPolymer):
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
    def __init__(self,species, name=None, base_species = None, material_type = "ordered_polymer", \
                             attributes = None, initial_concentration = 0,circular = False):

        self.material_type = material_type
        self.parent=None
        self.position=None
        self.direction=None
        
        self.initial_concentration = initial_concentration
        self.circular = circular

        if attributes is None:
            self.attributes = []
        else:
            self.attributes = attributes

        self._name = OrderedComplexSpecies._check_name(self,name)

        self.material_type = material_type
        #self.species = []
        monomers = []
        for specie in species:
            if isinstance(specie,Species) and isinstance(specie, OrderedMonomer):
                monomers += [specie]
            elif (isinstance(specie, tuple) or isinstance(specie, list)) and (isinstance(specie[0],Species) and isinstance(specie[0], OrderedMonomer)):
                monomers += [specie]
            elif isinstance(specie, OrderedPolymer):
                raise NotImplementedError(f"OrderedPolymer cannot be used as a Monomer at this time.")
            else:
                raise ValueError("{} should be a Species or list [Species, 'direction']".format(specie))
                #only species are acceptable
        
        OrderedPolymer.__init__(self, monomers)
        self.material_type = material_type
        
        if(base_species is None):
            self.base_species = Species(self.name, material_type = material_type)
        elif(isinstance(base_species,Species)):
            self.base_species = base_species
        else:
            raise TypeError("base_species is of type "+type(base_species)+" which is not acceptable. Use Species or str")
        
    @property
    def species_set(self):
        return set(self._polymer)
    @property
    def species(self):
        return self._polymer

    def get_species_list(self):
        return self._polymer
    
    @property
    def circular(self):
        if("circular" in self.attributes):
            return True
        else:
            return False

    @circular.setter
    def circular(self,value:bool):
        if(value):
            self.add_attribute("circular")
        else:
            self.remove_attribute("circular")

    def set_species_list(self,spec_tuple:tuple):
        OrderedPolymer.__init__(self,spec_tuple)

    @property
    def name(self):
        if self._name is None:
            name = ""
        else:
            name = self._name
        outlst = []

        for monomer in self._polymer:
            assert(isinstance(monomer,Species))
            pname = monomer.name
            pdir = None
            if(hasattr(monomer,"direction")):
                pdir =monomer.direction
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
        ophash += hash(self.circular)+hash(self.base_species)+hash(self.name)+hash(self.material_type)
        return ophash
    

    def replace(self,position,part,direction=None):
        #TODO only change the name if the part we are replacing is actually different
        mydir = direction
        if((mydir is None) and (part.direction is not None)):
            mydir = part.direction
        if(part == self._polymer[position] and self._polymer[position].direction == mydir):
            #in this case we are replacing a part with the same thing, so do nothing
            #but it could be true that the reference changes? That shouldnt be
            pass
        else:
            OrderedPolymer.replace(self,position=position,part=part,direction=mydir)
            #print("replacing")
            #print([a.data for a in self._polymer])
            #self.name = self.make_name()

    def __contains__(self,item):
        for part in self.species:
            if((part==item ) or (item == part.data) or (item in part)):
                return True
        return False