################################################################
#
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.
#
################################################################
import itertools as it
import numbers

from typing import Dict, Union
from ..core.species import WeightedSpecies, Species
from ..core.parameter import Parameter
from warnings import warn


def all_comb(input_list):
    out_list = []
    for i in range(1, len(input_list) + 1):
        out_list += it.combinations(input_list, i)
    return out_list


def rev_dir(dir):
    reversedict = {"forward": "reverse", "reverse": "forward"}
    return reversedict[dir]


def recursive_parent(s):
    # Recursively goes through Species and gets the top level parent
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
        # go through the species and replace species with their parents, recursively
        if isinstance(s, WeightedSpecies):
            parent = recursive_parent(s.species)
            s.species = parent
            out_sp_list.append(s)
        else:
            parent = recursive_parent(s)
            out_sp_list.append(parent)

    return out_sp_list


def parameter_to_value(p):
    if isinstance(p, Parameter):
        return p.value
    else:
        return p


def process_initial_concentration_dict(
    initial_concentration_dict: Dict[
        Union[str, Species], Union[numbers.Real, Parameter]
    ],
) -> Dict[str, numbers.Real]:
    processed_x0_dict = {
        str(key): parameter_to_value(value)
        for (key, value) in initial_concentration_dict.items()
    }
    return processed_x0_dict


def combine_dictionaries(dict1, dict2):
    """append lists that share the same key, and add new keys
    WARNING: this only works if the dictionaries have values that are lists"""
    outdict = dict1
    for key in dict2:
        if key in outdict:
            assert isinstance(dict2[key], list)
            assert isinstance(outdict[key], list)
            outdict[key] += dict2[key]
        else:
            outdict[key] = dict2[key]
    return outdict


def member_dictionary_search(member, dictionary):
    """searches through a dictionary for keys relevant to the given data member.
    Order of returning:
    repr
    name
    material_type
    propensity name
    propensity partid
    propensity mechanism
    """
    if dictionary is None:
        return None
    elif repr(member) in dictionary:
        return dictionary[repr(member)]
    elif hasattr(member, "name") and member.name in dictionary:
        return dictionary[member.name]
    elif (
        hasattr(member, "integrase")
        and hasattr(member, "site_type")
        and (str(member.integrase.name), member.site_type) in dictionary
    ):
        return dictionary[(str(member.integrase.name), member.site_type)]
    elif hasattr(member, "site_type") and member.site_type in dictionary:
        return dictionary[member.site_type]
    elif hasattr(member, "assembly"):
        typename = type(member).__name__
        if typename in dictionary:
            return dictionary[typename]
    elif (
        hasattr(member, "material_type")
        and hasattr(member, "attributes")
        and (member.material_type, tuple(member.attributes)) in dictionary
    ):
        return dictionary[(member.material_type, tuple(member.attributes))]
    elif hasattr(member, "material_type") and member.material_type in dictionary:
        return dictionary[member.material_type]
    elif hasattr(member, "attributes") and tuple(member.attributes) in dictionary:
        return dictionary[tuple(member.attributes)]
    elif hasattr(member, "propensity_type"):
        out_value = None
        try:
            for k, p in member.propensity_type.propensity_dict["parameters"].items():
                if hasattr(p, "search_key"):
                    mech_str = repr(p.search_key.mechanism).strip("'\"")
                    partid_str = repr(p.search_key.part_id).strip("'\"")
                    name_str = repr(p.search_key.name).strip("'\"")
                    cur_value = out_value
                    if name_str in dictionary:
                        # name of the mechanism that made the reaction
                        cur_value = dictionary[name_str]
                    if partid_str in dictionary:
                        # partid used to make the reaction
                        cur_value = dictionary[partid_str]
                    if mech_str in dictionary:
                        # the type of mechanism that made the reaction
                        cur_value = dictionary[mech_str]
                    if out_value is not None and cur_value != out_value:
                        warn(
                            f"dictionary search output was {out_value} but now it will be {cur_value}"
                        )
                    out_value = cur_value
        except KeyError:
            pass
        return out_value
    return None
