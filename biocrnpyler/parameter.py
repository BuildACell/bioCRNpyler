# parameter.py - parameter processing
# RMM, 19 Aug 2018
#
# This file contains the Parameter class that is used for representing
# parameters, as well as utility functions for manipulating
# parameters.
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.
"""
Parameter processing
--------------------

Parameters are used to set the rate constants for reactions, as well
as the initial concentrations of extract and buffer.

The following general guidelines hold for usage of parameters in the
`txtl` toolbox:

1. Parameter values for all pre-defined components (including extracts
   and buffers) should be contained in a config file, formatted as a
   CSV file with filename that matches the component name.  The
   configuration file will be searched for in the current directory,
   then in the directory of the component, and finally in the default
   TX-TL component library [not yet implemented].

2. Component parameter values are local to the reactions set up by the
   component, allowing reuse of the same name (eg, promoter binding
   strenth) across multiple parameters.  Extract parameter values are
   stored in the containing mixture and create global parameters in
   the SBML model for the system.  These global parameters can be
   accessed by using the extract parameter name within an txtl
   reaction.

3. Default configuration files can be overridden at the time that a
   component is defined by using the `config_file` argument.  If a
   pathname is not provided, the file will be searched for in the same
   order as default configuration files.  Examples:

     promoter = Promoter('pconst')
     promoter = Promoter('pconst', 'pconst.csv')
     promoter = Promoter('pconst', './configs/pconst.csv')

4. Individual parameter values can be overriden at the time that a
   component is defined using named arguments in the construction of
   the component or using the `parameters` argument (passed as a dict
   object, described below).  Examples:

     myptet = Promoter('ptet', RNAPbound_F=20)
     myptet = Promoter('ptet', parameters={'RNAPbound_F':20})
     ptet_param = {
       'RNAPbound_F' : 20,
       'RNAPbound_R' : Parameter('my_param', 'Numeric', 20)
     }
     myptet = Promoter('ptet', ptet_param)

5. Parameter dictionaries should have the parameter name as the
   dictionary key and the value can be one of three different object
   types:
     * A parameter object (from the Parameter class)
     * A floating point number
     * A string (representing a global parameter name)

6. Parameter objects can be of type 'Numeric', 'Expression', or
   'Global' [not implemented].  Numeric parameters are passed directly
   to libsbml.  Expression parameters are evaluated using the current
   dictionary, augmented by the following variables:

     * RNA_LENGTH - number of basepairs in the RNA sequence
     * AA_LENGTH = number of amino acides in the protein sequence

7. Parameters for all possible mechanisms that are defined for a
   component should be defined using the config file or parameter
   argument when the component is created.

"""

import csv
import os
import sys
import re
from warnings import warn


# Takes a list L and string s and returns a list of flexible variations of all
# the strings in S
def get_flexible_string_list_index(L, S):
    s_list = []
    for s in S:
        s_list += [s, s.replace(" ", "_"), s.capitalize(), s.casefold(),
                   s.casefold().replace(" ", "_")]
    ind = None
    for t in s_list:
        try:
            indt = L.index(t)
            if ind != None:
                raise ValueError("List contains multiple elements that are too "
                                 f"similar: '{str(L[ind])}', '{str(L[indt])}'")
            else:
                ind = indt
        except ValueError:
            pass
    return ind

#Duplicates in the new files overwrite current things in parameter dictionary
def create_parameter_dictionary(parameters, parameter_file):
    if not isinstance(parameters, dict) and parameters!= None:
        raise ValueError("parameter keyword must be a dictionary "
                         "{param_key:val}")
    if not isinstance(parameter_file, str) \
       and not isinstance(parameter_file, list) \
       and parameter_file != None:
        raise ValueError("parameter_file keyword must be a valid filepath "
                         "string or list of such strings.")

    if isinstance(parameter_file, list):
        file_list = parameter_file
    elif parameter_file != None:
        file_list = [parameter_file]
    elif parameter_file == None:
        file_list = []

    param_dict = {}
    if parameters != None:
        for k in parameters:
            param_dict[k] = parameters[k]

    for fname in file_list:

        f = open(fname)
        if fname[-4:] in [".txt",".tsv"]:
            L0 = f.readline().replace("\n", "").split("\t")
        else:
            L0 = f.readline().replace("\n", "").split(",")

        try:
            pID_ind = get_flexible_string_list_index(L0, ["part_id", "part"])
            if pID_ind == None:
                warn(f"file {fname} contains no part ID column. Please add a "
                     "column the name 'part_id' or 'part'.")
        except ValueError:
            raise ValueError("'part_id', 'part', or a similar string appears "
                            f"multiple times in the top line of {fname}. "
                            "keyword list = {L0}.")
        try:
            mech_ind = get_flexible_string_list_index(L0, ["mechanism",
                                                           "mechanism_id"])
            if mech_ind == None:
                warn(f"file {fname} contains no mechanism column. Please add a"
                     "column the name 'mechanism' or 'mechanism_id'.")

        except ValueError:
            raise ValueError("'mechanism', 'mechanism_id' or a similar string "
                             "appears multiple times in the top line of "
                             f"{fname}. keyword list = {L0}")
        try:
            param_ind = get_flexible_string_list_index(L0,
                        ["param_name", "parameter_name", "parameter", "param"])
            if param_ind == None:
                ValueError(f"file {fname} contains no parameter name column. "
                           "Please add a column the name 'param', 'parameter', "
                           "'param_name', or 'parameter_name'.")
        except ValueError:
            raise ValueError("'param', 'parameter', 'parameter_name', "
                             "'param_name', or a similar string appears "
                            f"multiple times in the top line of {fname}. "
                            "keyword list = {L0}")
        try:
            val_ind = get_flexible_string_list_index(L0,
                            ["val", "value", "param_val", "parameter_value"])
            if val_ind == None:
                raise ValueError(f"file {fname} contains no parameter value "
                                 "column. Please add a column the name 'val', "
                                 "'value', 'param_val', or 'parameter_value'.")
        except ValueError:
            raise ValueError("'val', 'value', 'param_val', 'parameter_value', "
                             "or a similar string appears multiple times in "
                            f"the top line of {fname}. keyword list = {L0}")

        for line in f:
            if fname[-4:] in [".txt", ".tsv"]:
                L = line.replace("\n", "").split("\t")
            else:
                L = line.replace("\n", "").split(",")

            if pID_ind != None:
                id = L[pID_ind]
            else:
                id = ""

            if mech_ind != None:
                mech = L[mech_ind]
            else:
                mech = ""

            param = L[param_ind]
            val = L[val_ind]

            try:
                if param == "":
                    pass
                elif mech == "" and id == "":
                    param_dict[param] = float(val)
                elif mech == "":
                    param_dict[(id, param)] = float(val)
                elif id == "":
                    param_dict[(mech, param)] = float(val)
                else:
                    param_dict[(mech, id, param)] = float(val)
            except ValueError:
                raise ValueError("Unable to parse parameter file line = "+line)
        f.close()
    return param_dict


class Parameter:
    """Parameter value (reaction rates)"""

    def __init__(self, name, param_type, value, comment="", debug=False):
        self.name = name.strip()
        self.param_type = param_type.strip()
        self.comment = comment.strip()

        # Set the value of the parameter
        if debug: print("%s [%s] = %s" % (self.name, self.param_type, value))
        if param_type.strip() == 'Numeric':
            self.value = float(value)  # store as float
        elif param_type.strip() == 'Expression':
            self.value = value  # store as string
        else:
            raise ValueError("can't parse value of parameter %s" % name)

    def get_value(self):
        return float(self.value)


def load_config(filename, extension=".csv", debug=False):
    # Find the configuration file
    # ! TODO: update this to search along a path (in pathutil)
    module_path = os.path.dirname(sys.modules[__name__].__file__)

    # Look for the config file in a list of paths
    csvfile = None
    for path in (module_path + "/components/", module_path + "/config/"):
        try:
            # ! TODO: add extension if not present
            filepath = path + filename
            csvfile = open(filepath)
            break
        except:
            continue

    # If we didn't find the file, return None
    if csvfile == None: return None

    # Open up the CSV file for reaching
    csvreader = csv.reader(csvfile)
    params = {}
    for row in csvreader:
        # Get rid of extraneous spaces
        for i in range(len(row)): row[i] = row[i].strip()

        # Skip blank lines (and malformed lines)
        if len(row) < 3 or row[0] == "": continue

        # Create a new parameter object to keep track of this row
        # ! TODO: this should be done in a better (and more pythonic) way
        if len(row) >= 4:
            #                 name  param_type value  comment
            param = Parameter(row[0], row[1], row[2], row[3])
        else:
            param = Parameter(row[0], row[1], row[2], "")

        # Name simplification for backward compatibility with MATLAB code
        param.name = re.sub("_Forward$", "_F", param.name)
        param.name = re.sub("_Reverse$", "_R", param.name)
        param.name = re.sub("_ic$", "_IC", param.name)
        param.name = re.sub("_Concentration$", "_IC", param.name)

        # Set up as dictionary for easy access
        params[param.name] = param

    csvfile.close()  # Close the file now that we are done
    return params  # Return the parameters we read


# Process parameter input
def get_parameters(config_file, custom, default={}, **keywords):
    # Start with the default parameters values (if given)
    parameters = default.copy() if default != None else {}

    # Now load parameters from the configuration file (if available)
    if config_file != None:
        config_parameters = load_config(config_file)
        if config_parameters != None:
            parameters.update(config_parameters)
        else:
            warn("get_parameters: couldn't find file %s" % config_file)

    # Override any parameters given as a parameter dictionary
    if custom != None:
        for key, value in custom.items():
            parameters[key] = _to_parameter(key, value)

    # Finally, check to see if any of the keywords are parameter names
    for key, value in keywords.items():
        if key in parameters.keys():
            parameters[key] = _to_parameter(key, value)

    # All done!
    return parameters


# Update any missing parameter values in a parameter dictionary
def update_missing(existing_dict, default_dict):
    """Fill in missing parameter values with defaults

    This takes a parameter dictionary, looks to see if a specified set
    of keys are missing, and, if so, fills in the keys with default
    values.  This function is useful when you have an existing
    parameter list and need to make sure that a certain set of default
    keys are present.  If they keys are already present in the
    parameter dictionary, they will not be replaced (unlike a regular
    dictionary update).

    """
    for key, value in default_dict.items():
        if key not in existing_dict.keys():
            existing_dict[key] = _to_parameter(key, value)


# Update any existing parameter values in a parameter dictionary
def update_existing(existing_dict, custom_dict):
    """Update existing parameter values with new values

    This takes a parameter dictionary, looks to see if a specified set
    of keys exist, and, if so, fills in the keys with new values.
    This function is useful when you have an existing parameter list
    that you want to override with parameter defintions from a longer
    list.

    """
    for key, value in custom_dict.items():
        if key in existing_dict.keys():
            existing_dict[key] = _to_parameter(key, value)


def eval_parameter(component, name, assignments={}):
    parameters = component.parameters
    if name not in parameters.keys() or parameters[name] == None:
        # Couldn't find the parmaeter
        return None
    param = parameters[name]

    # See if we already have a floating point number
    # ! TODO: decide if we need this; can just use the evaluation below?
    if isinstance(param.value, (float, int)): return float(param.value)

    # Evaluate the expression
    return float(eval(param.value, assignments))

# Convert a value input to a parameter object
def _to_parameter(key, value):
    if isinstance(value, Parameter):
        return value
    elif isinstance(value, (float, int)):
        return Parameter(key, 'Numeric', value)
    elif isinstance(value, str):
        return Parameter(key, 'Global', value)
    else:
        ValueError('Unknown parameter type')
