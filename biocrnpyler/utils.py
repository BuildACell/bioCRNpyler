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
        """go through every species on a list and remove any "bindloc" attributes. This is used
        to convert monomers with a parent polymer into the correct species after combinatorial binding
        in things like DNAassembly and RNAassembly."""
        
        out_sp_list = []
        for specie in spec_list:
            #waited species case used in reactions
            if(isinstance(specie,WeightedSpecies)):
                spec2 = specie.species
                #replace the species (representing a binding location) with its parent
                if(hasattr(spec2,"parent") and (spec2.parent is not None)):
                    specie.species = spec2.parent
            #OrderedMonomerSpecies is inside an OrderedPolymerSpecies
            if(hasattr(specie,"parent") and (specie.parent is not None)):
                #replace the Species with its parent
                out_sp_list += [specie.parent]
            #Standard Species (without parents) are not effected
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
