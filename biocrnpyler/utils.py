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


def remove_bindloc(spec_list):
        """go through every species on a list and remove any "bindloc" attributes"""
        #spec_list2 = copy.copy(spec_list)
        out_sp_list = []
        for specie in spec_list:
            #go through the species and remove the "bindloc" attribute
            #I don't care about the binding now that I am done generating species
            if(type(specie)==WeightedSpecies):
                spec2 = specie.species
                if(hasattr(spec2,"parent") and (spec2.parent is not None)):
                    specie.species = spec2.parent
            if(hasattr(specie,"parent") and (specie.parent is not None)):
                out_sp_list += [specie.parent]
            else:
                out_sp_list+= [specie]
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
