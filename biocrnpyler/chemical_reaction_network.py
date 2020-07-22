#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .sbmlutil import *
from .propensities import Propensity, MassAction
import warnings
import numpy as np
from typing import List, Union, Dict
from dataclasses import dataclass
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

        If the name is a starts with a number, there must be a material type.
        """
        if material_type in [None, ""] and self.name[0].isnumeric():
            raise ValueError(f"species name: {self.name} contains a number as the first character and therefore requires a material_type.")
        elif material_type is None:
            return ""
        elif (material_type.replace("_", "").isalnum() and material_type.replace("_", "")[0].isalpha()) or material_type == "":
                return material_type
        else:
            raise ValueError(f"material_type {material_type} must be alpha-numeric and start with a letter.")

    def check_name(self, name):
        """
        Check that the string contains only underscores and alpha-numeric characters
        """
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


@dataclass(frozen=True, eq=True)
class WeightedSpecies:
    """Container object for a species and its stoichiometry
    """
    species: Species
    stoichiometry: int = 1

    def pretty_print(self, **kwargs):
        return f'{self.stoichiometry if self.stoichiometry > 1 else ""}{self.species.pretty_print(**kwargs)}'

    def replace_species(self, *args, **kwargs):
        return self.species.replace_species(*args, **kwargs)

    @staticmethod
    def _count_weighted_species(weighted_species):
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


class Reaction(object):
    """ An abstract representation of a chemical reaction in a CRN
    A reaction has the form:
    .. math::
       \sum_i n_i I_i --> \sum_i m_i O_i @ rate = k
       where n_i is the count of the ith input, I_i, and m_i is the count of the
       ith output, O_i.
    If the reaction is reversible, the reverse reaction is also included:
    .. math::
       \sum_i m_i O_i  --> \sum_i n_i I_i @ rate = k_rev
    """
    def __init__(self, inputs: Union[List[Species], List[WeightedSpecies]],
                 outputs: Union[List[Species], List[WeightedSpecies]],
                 propensity_type: Propensity):

        if len(inputs) == 0 and len(outputs) == 0:
            warn("Reaction Inputs and Outputs both contain 0 Species.")

        self.inputs = Species.flatten_list(inputs)
        self.outputs = Species.flatten_list(outputs)
        self.propensity_type = propensity_type

    @property
    def propensity_type(self) -> Propensity:
        return self._propensity_type

    @propensity_type.setter
    def propensity_type(self, new_propensity_type: Propensity):
        """ Replace the propensity type associated with the reaction object

        :param new_propensity_type: Valid propensity type
        """
        if not Propensity.is_valid_propensity(new_propensity_type):
            raise ValueError(f'unknown propensity type: {new_propensity_type} '
                             f'({type(new_propensity_type)})!')

        self._propensity_type = new_propensity_type

    @classmethod
    def from_massaction(cls, inputs: Union[List[Species], List[WeightedSpecies]],
                        outputs: Union[List[Species], List[WeightedSpecies]],
                        k_forward: float, k_reverse: float = None):
        """ Initialize a Reaction object with mass action kinetics
        :param inputs:
        :param outputs:
        :param k_forward:
        :param k_reverse:
        :return: Reaction object
        """
        mak = MassAction(k_forward=k_forward, k_reverse=k_reverse)

        return cls(inputs=inputs, outputs=outputs, propensity_type=mak)

    @property
    def is_reversible(self) -> bool:
        return self.propensity_type.is_reversible

    @property
    def inputs(self) -> List[WeightedSpecies]:
        return self._input_complexes

    @inputs.setter
    def inputs(self, new_input_complexes: List[WeightedSpecies]):
        self._input_complexes = Reaction._check_and_convert_complex_list(complexes=new_input_complexes)

    @property
    def outputs(self) -> List[WeightedSpecies]:
        return self._output_complexes

    @outputs.setter
    def outputs(self, new_output_complexes: List[WeightedSpecies]):
        self._output_complexes = Reaction._check_and_convert_complex_list(complexes=new_output_complexes)

    @staticmethod
    def _check_and_convert_complex_list(complexes: Union[List[Species], List[WeightedSpecies]]) -> List[WeightedSpecies]:
        if all([isinstance(one_complex, Species) for one_complex in complexes]):
            # we wrap each Species object to WeightedSpecies
            complexes = [WeightedSpecies(species=species) for species in complexes]
        else:
            if not all([isinstance(one_complex, WeightedSpecies) for one_complex in complexes]):
                raise TypeError(f'inputs must be list of Species or list of ChemicalComplexes!')

        # filter out duplicates and adjust stoichiometry
        out_list = []
        # Create a dictionary of unique species and their stoichiometry count
        stoichiometry_count = WeightedSpecies._count_weighted_species(complexes)
        for one_complex, stoichiometry in stoichiometry_count.items():
            new_complex = WeightedSpecies(species=one_complex.species, stoichiometry=stoichiometry)
            out_list.append(new_complex)

        return out_list

    @property
    def k_forward(self):
        return self.propensity_type.k_forward

    @property
    def k_reverse(self):
        return self.propensity_type.k_reverse

    def replace_species(self, species: Species, new_species: Species):
        """Replaces species with new_species in the reaction
        :param species: species to be replaced
        :param new_species: the new species the old species is replaced with
        :return: a new Reaction instance
        """
        if not isinstance(species, Species) or not isinstance(new_species, Species):
            raise ValueError('both species and new_species argument must be an instance of Species!')

        new_inputs = []
        for s in self.inputs:
            new_s = s.replace_species(species, new_species)
            new_inputs.append(new_s)

        new_outputs = []
        for s in self.outputs:
            new_s = s.replace_species(species, new_species)
            new_outputs.append(new_s)

        # get a shallow copy of the parameters and species, so we can replace some of them
        propensity_type_dict = copy.copy(self.propensity_type.propensity_dict)
        for key, prop_species in propensity_type_dict['species'].items():
            propensity_type_dict['species'][key] = prop_species.replace_species(species, new_species)

        new_propensity_type = self.propensity_type.from_dict(propensity_type_dict)
        print(new_propensity_type.propensity_dict)

        new_r = Reaction(inputs=new_inputs, outputs=new_outputs, propensity_type=new_propensity_type)
        return new_r

    def __repr__(self):
        """Helper function to print the text of a rate function"""
        return self.pretty_print(show_rates=False, show_material=True, show_attributes=False)

    def pretty_print(self, show_rates=True, show_material=True, show_attributes=True, **kwargs):

        kwargs['show_rates'] = show_rates
        kwargs['show_material'] = show_material
        kwargs['show_attributes'] = show_attributes

        txt = '+'.join([s.pretty_print(**kwargs) for s in self.inputs])

        if self.is_reversible:
            txt += " <--> "
        else:
            txt += " --> "

        txt += '+'.join([s.pretty_print(**kwargs) for s in self.outputs])
        if show_rates:
            txt += f'{self.propensity_type.pretty_print(**kwargs)}'

        return txt

    def __eq__(self, other):
        """Two reactions are equivalent if they have the same inputs, outputs,
           and propensity."""
        if not isinstance(other, Reaction):
            raise TypeError(f'Only reactions can be compared with reaction! We got {type(other)}.')

        if len(self.inputs) != len(other.inputs) or len(self.outputs) != len(other.outputs):
            return False

        return (set(self.inputs), set(self.outputs), self.propensity_type) == (set(other.inputs), set(other.outputs), other.propensity_type)

    def __contains__(self, item: Species):
        """It checks whether a species is part of a reaction.

         it checks the input and output lists as well as the propensity type for the species

        :param item: a Species instance
        :return: bool
        :exception NotImplementedError for non-Species objects
        """
        if isinstance(item, Species):
            if item in self.inputs \
                    or item in self.outputs \
                    or item in self.propensity_type.species:
                return True
        else:
            raise NotImplementedError
        return False

    @property
    def species(self) -> List[Species]:
        """returns a list of species in the reactions collected from the inputs
        and outputs and the propensity (e.g. Hill kinetics has species in it)

        :return: list of species in the reactions
        """
        in_part = []
        for s in self.inputs:
            in_part.extend(Species.flatten_list(s.species))
        out_part = []
        for s in self.outputs:
            out_part.extend(Species.flatten_list(s.species))

        return list(itertools.chain(in_part, out_part, self.propensity_type.species))


class ChemicalReactionNetwork(object):
    """ A chemical reaction network is a container of species and reactions
    chemical reaction networks can be compiled into SBML or represented
    conveniently as python tuple objects.
    reaction types:
       mass action: standard mass action semantics where the propensity of a
                reaction is given by deterministic propensity =
       .. math::
                k \Prod_{inputs i} [S_i]^a_i
               stochastic propensity =
        .. math::
                k \Prod_{inputs i} (S_i)!/(S_i - a_i)!
               where a_i is the spectrometric coefficient of species i
    """
    def __init__(self, species: List[Species], reactions: List[Reaction], show_warnings=False):
        self.reactions, self.species = ChemicalReactionNetwork.check_crn_validity(reactions, species, show_warnings=show_warnings)

    def add_species(self, species, show_warnings=False):
        if not isinstance(species, list):
            species = [species]
        self.reactions, self.species = ChemicalReactionNetwork.check_crn_validity(self.reactions, list(itertools.chain(self.species, species)), show_warnings=show_warnings)

    def add_reactions(self, reactions, show_warnings=False):
        if not isinstance(reactions, list):
            reactions = [reactions]
        self.reactions, self.species = ChemicalReactionNetwork.check_crn_validity(list(itertools.chain(self.reactions, reactions)), self.species, show_warnings=show_warnings)

    @staticmethod
    def check_crn_validity(reactions: List[Reaction], species: List[Species], show_warnings=True):

        if not all(isinstance(r, Reaction) for r in reactions):
            raise ValueError("A non-reaction object was used as a reaction!")

        if not all(isinstance(s, Species) for s in species):
            raise ValueError("A non-species object was used as a species!")

        for r in reactions:
            if reactions.count(r) > 1 and show_warnings:
                warn(f"Reaction {r} may be duplicated in CRN definitions. "
                     f"Duplicates have NOT been removed.")

        for s in species:
            if species.count(s) > 1 and show_warnings:
                warn(f"Species {s} is duplicated in the CRN definition. "
                     f"Duplicates have NOT been removed.")

        # check that all species in the reactions are also in the species list and vice versa
        unique_species = set(species)
        all_species_in_reactions = set(Species.flatten_list([r.species for r in reactions]))
        if unique_species != all_species_in_reactions:
            species_without_reactions = unique_species - all_species_in_reactions
            if species_without_reactions and show_warnings:
                warn(f'These Species {list(species_without_reactions)} are not part of any reactions in the CRN!')
            unlisted_reactions = all_species_in_reactions - unique_species
            if unlisted_reactions and show_warnings:
                warn(f'These Species {list(unlisted_reactions)} are not listed in the Species list, but part of the reactions!')

        return reactions, species

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
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed
        show_rates toggles whether reaction rate functions are printed
        """

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

    def initial_condition_vector(self, init_cond_dict: Union[Dict[str, float], Dict[Species, float]]):
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

    def replace_species(self, species: Species, new_species: Species):
        """Replaces species with new_species in the entire CRN.

        Does not act in place: returns a new CRN.
        """

        if not isinstance(species, Species):
            raise ValueError('species argument must be an instance of Species!')

        if not isinstance(new_species, Species):
            raise ValueError('species argument must be an instance of Species!')

        new_species_list = []
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_reaction_list = []
        for r in self.reactions:
            new_r = r.replace_species(species, new_species)
            new_reaction_list.append(new_r)

        return ChemicalReactionNetwork(new_species_list, new_reaction_list)

    def generate_sbml_model(self, stochastic_model=False, **keywords):
        """Creates an new SBML model and populates with the species and
        reactions in the ChemicalReactionNetwork object

        :param stochastic_model: whether the model is stochastic
        :param keywords: extra keywords pass onto create_sbml_model()
        :return: tuple: (document,model) SBML objects
        """
        document, model = create_sbml_model(**keywords)

        add_all_species(model=model, species=self.species)

        add_all_reactions(model=model, reactions=self.reactions,
                          stochastic=stochastic_model)

        if document.getNumErrors():
            warn('SBML model generated has errors. Use document.getErrorLog() to print all errors.')
        return document, model

    def write_sbml_file(self, file_name=None, **keywords) -> bool:
        """Writes CRN to an SBML file
        """
        document, _ = self.generate_sbml_model(**keywords)
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, 'w') as f:
            f.write(sbml_string)
        return True

    def create_bioscrape_model(self):
        """Creates a Bioscrape Model of the CRN directly.
        """
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

            if rxn.is_reversible and rxn.propensity_type == "massaction":
                prop_params['k'] = rxn.k_r
                reaction_list.append((products, reactants, prop_type,
                                      dict(prop_params)))
            elif rxn.is_reversible:
                raise ValueError("Only massaction irreversible reactions are "
                                 "supported for automatic bioscrape simulation."
                                 " Consider creating two seperate reactions.")
        model = Model(species = species_list, reactions = reaction_list,
                      initial_condition_dict = initial_condition_dict)
        return model

    def simulate_with_bioscrape(self, timepoints, initial_condition_dict=None,
                                stochastic = False, return_dataframe = True,
                                safe = False, **kwargs):

        """Simulate CRN model with bioscrape (https://github.com/biocircuits/bioscrape).
        Returns the data for all species as Pandas dataframe.
        """
        result = None
        try:
            from bioscrape.simulator import py_simulate_model
            m = self.create_bioscrape_model()
            if not initial_condition_dict:
                initial_condition_dict = {}
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
                initial_condition_dict = None, return_dataframe = True,
                stochastic = False, **kwargs):

        """Simulate CRN model with bioscrape via writing a SBML file temporarily.
        [Bioscrape on GitHub](https://github.com/biocircuits/bioscrape).

        Returns the data for all species as Pandas dataframe.
        """
        result = None
        m = None
        try:
            from bioscrape.simulator import py_simulate_model
            from bioscrape.types import Model
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
            m = Model(sbml_filename = file_name, sbml_warnings = sbml_warnings)
            # m.write_bioscrape_xml('temp_bs'+ file_name + '.xml') # Uncomment if you want a bioscrape XML written as well.
            m.set_species(initial_condition_dict)
            result = py_simulate_model(timepoints, Model = m,
                                                stochastic = stochastic,
                                                return_dataframe = return_dataframe)
        except ModuleNotFoundError:
            warnings.warn('bioscrape was not found, please install bioscrape')

        return result, m

    def runsim_roadrunner(self, timepoints, filename, species_to_plot = None):
        """To simulate using roadrunner.
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
