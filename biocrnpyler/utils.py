################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 7/31/2020
#       Copyright (c) 2020, Build-A-Cell. All rights reserved.
#       See LICENSE file in the project root directory for details.
#
################################################################
import itertools as it
import numbers

from typing import Dict, Union
from .species import WeightedSpecies, Species
from .parameter import Parameter


def all_comb(input_list):
    out_list = []
    for i in range(1,len(input_list)+1):
        out_list += it.combinations(input_list,i)
    return out_list


def rev_dir(dir):
    reversedict = {"forward":"reverse","reverse":"forward"}
    return reversedict[dir]

def recursive_parent(s):
    #Recursively goes through Species and gets the top level parent
    if hasattr(s, "parent") and s.parent is not None:
        return recursive_parent(s.parent)
    else:
        return s

def remove_bindloc(spec_list):
        """go through every species on a list and remove any "bindloc" attributes. This is used
        to convert monomers with a parent polymer into the correct species after combinatorial binding
        in things like DNAassembly and RNAassembly."""
        if not isinstance(spec_list, list):
            spec_list = [spec_list]

        out_sp_list = []
        for s in spec_list:
            #go through the species and replace species with their parents, recursively

            if isinstance(s, WeightedSpecies):
                parent = recursive_parent(s.species)
                s.species = parent
                out_sp_list.append(s)
            else:
                parent = recursive_parent(s)
                out_sp_list.append(parent)
                
        return out_sp_list

#Converts a parameter to its Value
def parameter_to_value(p):
    if isinstance(p, Parameter):
        return p.value
    else:
        return p

#Converts a dictionary of Species (or strings) --> Parameters (or Numbers) to Strings --> Numbers
def process_initial_concentration_dict(initial_concentration_dict: Dict[Union[str, Species], Union[numbers.Real, Parameter]]) -> Dict[str, numbers.Real]:
    return {str(key): parameter_to_value(value) for (key, value) in initial_concentration_dict.items()}
