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
#### Parameter Value Defaulting:
Not all parameters need to have the required headings. 
    The only two required columns are "param_val" and "param_name". 
    BioCRNpyler uses a form of parameter name defaulting discussed below to find default parameters 
    if no exact match is in the config file. This makes it easy to set default parameters for things 
    like "ku" and "ktx" to quickly build models.

#### Parameters inside BioCRNpyler:
Inside of bioCRNpyler, parameters are stored as a dictionary key value pair: 
    (mechanism_name, part_id, param_name) --> param_val. If that particular parameter key cannot be found, 
    the software will default to the following keys: 
    (mechanism_type, part_id, param_name) >> (part_id, param_name) >> 
    (mechanism_name, param_name) >> (mechanism_type, param_name) >>
       (param_name) and give a warning. 
    As a note, mechanism_name refers to the .name variable of a Mechanism. mechanism_type refers to the .type variable of a Mechanism. 
    Either of these can be used as a mechanism_id. This allows for models to be constructed easily using default parameter values and 
    for parameters to be shared between different Mechanisms and/or Components.

#### Initial Conditions are also Parameters
The initial condition of any Species (or Component) will also be looked up as a parameters automatically.
    Initial conditions can be customized in through the custom_initial_condition keyword in the Mixture constructor.
    custom_initial_conditions will take precedent to parameter initial conditions.

    During compilation, Mixture.set_initial_condition() checks for parameters for all species in the following order:

    # First checks if (mixture.name, repr(species) is in the self.custom_initial_condition_dict
    # Then checks if (repr(species) is in the self.custom_initial_condition_dict
    # Then checks if (mixture.name, component.name) is in the self.custom_initial_condition_dictionary
    # Then checks if (component.name) is in the self.custom_initial_condition_dictionary

    # Then checks if (mixture.name, repr(species) is in the parameter dictionary
    # Then checks if repr(species) is in the parameter dictionary
    # Then checks if (mixture.name, component.name) is in the parameter dictionary
    # Then checks if component.name is in the parameter dictionary
    # Then defaults to 0
"""

import csv
from warnings import warn
from typing import List, Dict, Union
import numbers
import re


class Parameter(object):
    def __init__(self, parameter_name: str, parameter_value: Union[str, numbers.Real]):
        """A class for representing parameters in general. Only the below subclasses are ever used.

        :param parameter_name: is the name of the parameter
        :param parameter_value: is the value of the parameter
        """
        self.parameter_name = parameter_name
        self.parameter_value = parameter_value

    @property
    def parameter_name(self):
        return self._parameter_name

    @parameter_name.setter
    def parameter_name(self, new_parameter_name):
        if not isinstance(new_parameter_name, str):
            raise ValueError(f"parameter_name must be a string: received {type(new_parameter_name)}.")
        if not re.search('^[a-z]+', new_parameter_name, re.IGNORECASE):
            raise ValueError(f'parameter_name should be at least one character and cannot start with a number!')

        self._parameter_name = new_parameter_name

    @property
    def parameter_value(self):
        return self._parameter_value

    @parameter_value.setter
    def parameter_value(self, new_parameter_value):
        if not (isinstance(new_parameter_value, numbers.Real) or isinstance(new_parameter_value, str)):
            raise ValueError(f"parameter_value must be a float or int: received {type(new_parameter_value)}.")
        if isinstance(new_parameter_value, str):
            print(re.search('(^[1-9]+/[1-9]+)|(^[1-9]+e[0-9]+)|(^[0-9])', new_parameter_value, re.I))
            print(re.search('[a-d-f-z]', new_parameter_value, re.I))
            if re.search('[a-d-f-z]', new_parameter_value, re.I) \
                    or re.search('(^[1-9]+/[1-9]+)|(^[1-9]+e[0-9]+)|(^[0-9])', new_parameter_value, re.I) is None:
                raise ValueError('No valid parameter value! Accepted formats: 1.00 or 1e4 or 2/5 ')

            self._parameter_value = Parameter._convert_rational(new_parameter_value)
        else:
            self._parameter_value = new_parameter_value

    @staticmethod
    def _convert_rational(p_value):
        if '/' in p_value:
            nom, denom = p_value.split('/')
            return float(nom)/float(denom)
        else:
            return float(p_value)

#A class for representing parameters in a parameter stored the ParameterDatabase
# parameter_keys is a dictionary {key:value} of keys for looking up the parameter
# parameter_info is a dictionary {key:value} of additional information about the parameter. 
#     For example: additional columns in the parameter file or the parameter file name.
class ParameterEntry(Parameter):
    def __init__(self, parameter_name, parameter_value, parameter_keys = None, parameter_info = None):
        Parameter.__init__(self, parameter_name, parameter_value)

        
        self.key_tuple = (None, None, parameter_name)

        if parameter_keys is None:
            self.keys = {}
        elif isinstance(parameter_keys, dict):
            self.keys = dict(parameter_keys)

            if "mechanism" in self.keys and "part_id" in self.keys:
                self.key_tuple = (self.keys["mechanism"], self.keys["part_id"], parameter_name)
            elif "part_id" in self.keys:
                self.key_tuple = (None, self.keys["part_id"], parameter_name)
            elif "mechanism" in self.keys:
                self.key_tuple = (self.keys["mechanism"], None, parameter_name)
        else:
            raise ValueError(f"parameter_keys must be a None or a dictionary: received {parameter_keys}.")

        if parameter_info is None:
            self.info = {}
        elif isinstance(parameter_info, dict):
            self.info = dict(parameter_info)
        else:
            raise ValueError(f"parameter_info must be None or a dictionary: received {parameter_info}.")

#A class for representing parameters used in the Model
#  search_key is a tuple searched for to find the parameter, eg (mech_id, part_id, param_name), :
#  found_key is the tuple used after defaulting to find the parameter eg (param_name)
class ModelParameter(ParameterEntry):
    def __init__(self, parameter_name, parameter_value, search_key, found_key, parameter_keys = None, parameter_info = None):

        ParameterEntry.__init__(self, parameter_name, parameter_value, parameter_keys = parameter_keys, parameter_info = parameter_info)

        if not isinstance(search_key, tuple):
            raise ValueError(f"search_key must be a tuple: received {search_key}.")
        else:
            self.search_key = search_key

        if not isinstance(found_key, tuple):
            raise ValueError(f"found_key must be a tuple: recieved {found_key}.")
        else:
            self.found_key = found_key


#A class for storing parameters in Components and Mixtures
class ParameterDatabase():
    def __init__(self, parameter_dictionary = None, parameter_file = None, overwrite_parameters = False):

        self.parameters = {} #create an emtpy dictionary to get parameters.

        if isinstance(parameter_file, str):
            self.load_parameters_from_file(parameter_file, overwrite_parameters = overwrite_parameters)
        elif isinstance(parameter_file, list):
            for p in parameter_file:
                if isinstance(p, str):
                    self.load_parameters_from_file(p, overwrite_parameters = overwrite_parameters)
                else:
                    raise ValueError("parameter_file must be a string or list of strings representing file names (and paths).")
        elif parameter_file is not None:
            raise ValueError("parameter_file must be a string representing a file name (and path).")
            
        if isinstance(parameter_dictionary, dict):
            self.load_parameters_from_dictionary(parameter_dictionary, overwrite_parameters = overwrite_parameters)
        elif parameter_dictionary is not None:
            raise ValueError("parameter_dictionary must be None or a dictionary!")

    #To check if a key or ParameterEntry is in a the ParameterDatabase
    def __contains__(self, val):
        if isinstance(val, ParameterEntry):
            key = val.key_tuple
            if key in self.parameters and self.parameters[key] == val:
                return True
            else:
                return False
        else:
            if isinstance(val, str):
                key = (None, None, val)
            else:
                key = val

            if key in self.parameters:
                return True
            else:
                return False


    #Ability to loop through parameters eg
    # for entry in ParameterDatabase: ...
    def __iter__(self):
        self.keys = list(self.parameters.keys())
        self.current_key_ind = 0
        return self
    def __next__(self):
        if self.current_key_ind < len(self.keys):
            key = self.keys[self.current_key_ind]
            entry = self.parameters[key]
            self.current_key_ind += 1
            return entry
        else:
            raise StopIteration

    #Gets a parameter from the database
    #Only returns exact matches.
    def __getitem__(self, key):
        param = None
        if isinstance(key, str):
            if (None, None, key) in self.parameters:
                param = self.parameters[(None, None, key)]
        elif isinstance(key, tuple) and len(key) == 3:
            if key in self.parameters:
                param = self.parameters[key]
        else:
            raise ValueError(f"Invalid parameter key {key}. Key must be be a tuple ('mechanism', 'part_id', 'parameter_name') or a string 'parameter_name'.")

        if param == None:
            raise ValueError(f"Parameter {key} not in ParameterDatabase.")
        else:
            return param

    #Sets a parameter in the databases - useful for quickly changing parameters, but add_parameter is recommended.
    def __setitem__(self, key, value):

        if isinstance(key, str):
            parameter_name = key
            parameter_keys = None
            key = (None, None, parameter_name)
        elif isinstance(key, tuple) and len(key) == 3:
            parameter_name = key[2]
            parameter_keys = {"mechanism":key[0], "part_id":key[1]}
        else:
            raise ValueError(f"Invalid parameter key {key}. Key must be be a tuple ('mechanism', 'part_id', 'parameter_name') or a string 'parameter_name'.")        
        
        if isinstance(value, ParameterEntry):
            if key != value.key_tuple:
                raise ValueError(f"Parameter Key does not match: ParameterDatabase key {key} is not the same as ParameterEntry Key {value.key_tuple}.")
            self.parameters[key] = value
        else:
            self.add_parameter(parameter_name, value, parameter_keys = parameter_keys, parameter_origin = "Set Manually", overwrite_parameters = True)
            self.add_parameter(parameter_name, value, parameter_origin = "Set Manually", overwrite_parameters = True)
    
    #Adds a parameter to the database with appropriate metadata
    def add_parameter(self, parameter_name, parameter_value, parameter_origin = None, parameter_keys = None, parameter_info = None, overwrite_parameters = False):

        #Put parameter origin into parameter_info
        if parameter_info is None:
            parameter_info = {}
        elif not isinstance(parameter_info, dict):
            raise ValueError("parameter_info must be None or a dictionary!")

        if parameter_origin is not None:
            parameter_info["parameter origin"] = parameter_origin

        #Create ParameterEntry
        param = ParameterEntry(parameter_name, parameter_value, parameter_keys = parameter_keys, parameter_info = parameter_info)

        #Determine correct parameter key
        if parameter_keys is None or (isinstance(parameter_keys, dict) and len(parameter_keys) == 0):
            key = (None, None, parameter_name)
        elif "mechanism" in parameter_keys and "part_id" in parameter_keys:
            key = (parameter_keys["mechanism"], parameter_keys["part_id"], parameter_name)
        elif "mechanism" in parameter_keys:
            key = (parameter_keys["mechanism"], None, parameter_name)
        elif "part_id" in parameter_keys:
            key = (None, parameter_keys["part_id"], parameter_name)
        else:
            raise ValueError(f"Invalid Parameter Keys {parameter_keys} for parameter ({parameter_name} = {parameter_value}): expected 'mechanism' and/or 'part_id'.")

        #Update parameter dictionary
        if key in self.parameters and not overwrite_parameters:
            raise ValueError(f"Duplicate parameter detected. Parameter with key = {key} is already in the ParameterDatabase. To Overwrite existing parameters, use overwrite_parameters = True.")
        else:
            self.parameters[key] = param

    #Loads Parameters from a dictionary
    def load_parameters_from_dictionary(self, parameter_dictionary, overwrite_parameters = False):
        for k in parameter_dictionary:
            if isinstance(k, str):
                self.add_parameter(k, parameter_dictionary[k], parameter_origin = "parameter_dictionary", overwrite_parameters = overwrite_parameters)

            elif isinstance(k ,tuple) and len(k) == 3:
                self.add_parameter(k[2], parameter_dictionary[k], parameter_keys = {"part_id":k[1], "mechanism":k[0]}, parameter_origin = "parameter_dictionary", overwrite_parameters = overwrite_parameters)
            else:
                raise ValueError(f"Invalid parameter key {k}. parameter_dictionary keys must be a 3-tuple or a single string: ('mechanism', 'part_id', 'param_name') OR 'param_name'. Note 'mechanism' and 'part_id' can be None.")

    #Loads parameters from another ParameterDatabase
    def load_parameters_from_database(self, parameter_database, overwrite_parameters = False):

        if not isinstance(parameter_database, ParameterDatabase):
            raise TypeError(f"paramater_database must be a ParamaterDatabase: recievied {parameter_database}.")

        for k in parameter_database:
            if k not in self.parameters or overwrite_parameters:
                self.parameters[k.key_tuple] = parameter_database[k.key_tuple]
            else:
                raise ValueError(f"Duplicate parameter detected. Parameter with key = {key} is already in the ParameterDatabase. To Overwrite existing parameters, use overwrite_parameters = True.")


    def load_parameters_from_file(self, filename, overwrite_parameters = False):

        #Figure out the format of the parameter file
        with open(filename) as f:
            file_type = filename.split(".")[-1]
            if file_type in ["tsv", "txt"]:
                delimiter = '\t'
            elif file_type in ["csv"]:
                delimiter = ","
            else:
                raise ValueError("Parameter files must be tab-seperated (.tsv or .txt) or comma-seperated (.csv) files.")

            csvreader = csv.DictReader(f, delimiter=delimiter)
            #Used for flexible column headings
            accepted_field_names = {
                'mechanism': ['mechanism', 'mechanism_id'],
                'param_name': ["parameter_name", "parameter", "param", "param_name"],
                'part_id': ['part_id', 'part'],
                'param_val': ["val", "value", "param_val", "parameter_value"]
            }

            field_names = self._get_field_names(csvreader.fieldnames, accepted_field_names)

            #Determine which columns are in the CSV
            if field_names['param_name'] is None:
                warn('No param name column was found, could not load parameter')
                return param_dict
            if field_names['mechanism'] is None:
                no_mechism_column = True
            else:
                no_mechism_column = False

            if field_names['part_id'] is None:
                no_part_id_column = True
            else:
                no_part_id_column = False

            #Load all parameters
            for row in csvreader:
                # TODO what about integers? float might cause numerical drift in simulations, e.g. cooperativity=2.001
                param_value = float(row[field_names['param_val']])
                field_columns = [field_names['param_name'], field_names['part_id'], field_names['mechanism'], field_names['param_val']]
                parameter_info = {k:row[k] for k in row if k not in field_columns}
                # TODO test all these cases!

                #Case 1: No Param Name so skip the row
                if row[field_names['param_name']] is None or len(row[field_names['param_name']]) == 0:
                    pass

                #Case 2: Just a Param Name
                elif no_mechism_column and no_part_id_column:
                    param_name = row[field_names['param_name']]
                    self.add_parameter(param_name, param_value, parameter_origin = filename, 
                        parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                #Case 3: Part_id and Param Name
                elif no_mechism_column and no_part_id_column is False:
                    param_name = row[field_names['param_name']]
                    part_id = row[field_names['part_id']]

                    if part_id is not None and len(part_id) > 0:
                        self.add_parameter(param_name, param_value, parameter_keys = {"part_id":part_id}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info , overwrite_parameters = overwrite_parameters)

                #Case 4: mechanism and param name
                elif no_part_id_column and no_mechism_column is False:
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if mech_name is not None and len(mech_name) > 0:
                        self.add_parameter(param_name, param_value, parameter_keys = {"mechanism":mech_name}, parameter_origin = filename, 
                        parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                #Case 5: mechanism, part_id, and param name
                else:
                    part_id = row[field_names['part_id']]
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if part_id is not None and len(part_id) > 0 and mech_name is not None and len(mech_name) >0:
                        self.add_parameter(param_name, param_value, parameter_keys = {"part_id":part_id, "mechanism":mech_name}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    
                    elif part_id is not None and len(part_id) > 0:    
                        self.add_parameter(param_name, param_value, parameter_keys = {"part_id":part_id}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                    elif mech_name is not None and len(mech_name) >0:
                        self.add_parameter(param_name, param_value, parameter_keys = {"mechanism":mech_name}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)


    @staticmethod
    def _get_field_names(field_names: List[str], accepted_field_names: Dict[str, List[str]]) -> Dict[str, str]:
        """ Searches through valid field names and finds the currently used one. It builds a dictionary of currently
            used field names
        :param field_names: list of field names (columns) found in the csv file
        :param accepted_field_names: dictionary of possible field names and their valid aliases
        :return: dictionary of currently used field names (aliases)
        """
        if not isinstance(field_names, list):
            raise ValueError('field_names must be a list of strings')
        if isinstance(field_names, list) and len(field_names) == 0:
            raise ValueError('field_names cannot be empty list!')
        if not isinstance(accepted_field_names, dict):
            raise ValueError('accepted_field_names must be a dictionary')
        if isinstance(accepted_field_names, dict) and len(accepted_field_names) == 0:
            raise ValueError('accepted_field_names cannot be empty dictionary')

        return_field_names = dict.fromkeys(accepted_field_names.keys())
        for accepted_name in accepted_field_names:
            # try to find an possible accepted names in the field_names using a generator
            try:
                loc_gen = (idx for idx, name in enumerate(accepted_field_names[accepted_name]) if name in field_names)
                loc_idx = next(loc_gen)
            except StopIteration:
                # we have reached the end of the possible names
                return_field_names[accepted_name] = None
                warn(f"parameter file contains no {accepted_name} column! Please add a "
                    f"column named {accepted_field_names[accepted_name]}.")
            else:
                return_field_names[accepted_name] = accepted_field_names[accepted_name][loc_idx]

        return return_field_names

    #Searches the database for the best matching parameter.
    def find_parameter(self, mechanism, part_id, param_name, parameter_warnings = False):

        found_entry = None
        warning_txt = None

        if isinstance(mechanism, str):
            mech_name = mechanism
            mech_type = mechanism
        elif hasattr(mechanism, "name") and hasattr(mechanism, "mechanism_type"):
            mech_name = mechanism.name
            mech_type = mechanism.mechanism_type
        elif mechanism is not None:
            raise ValueError(f"mechanism keyword must be or string or have name and mechanism_type attributes: recievied {mechanism}.")
        # Ideally parameters can be found
        # (mechanism.name/mechanism_type, part_id, param_name) --> val
        if part_id is not None and mechanism is not None:
            if mechanism is not None and (mechanism.name, part_id, param_name) in self.parameters:
                found_entry = self.parameters[(mech_name, part_id, param_name)]
                found_key = (mech_name, part_id, param_name)
            elif mechanism is not None and (mech_type, part_id, param_name) in self.parameters:
                found_entry = self.parameters[(mech_type, part_id, param_name)]
                found_key = (mech_type, part_id, param_name)

        # Next try (part_id, param_name) --> val
        if part_id is not None and found_entry is None:
            if (None, part_id, param_name) in self.parameters:
                found_entry = self.parameters[(None, part_id, param_name)]
                found_key = (None, part_id, param_name)

                warning_txt = ("No Parameter found with "
                    f"param_name={param_name} and part_id={part_id} and "
                    f"mechanism={repr(mechanism)}. Parameter found under "
                    f"the key (part_id, param_name)=({part_id}, "
                    f"{param_name}).")
        # Next try (Mechanism.name/mechanism_type, param_name) --> val
        if mechanism is not None and found_entry is None:
            if (mechanism.name, None, param_name) in self.parameters:
                found_entry = self.parameters[(mech_name, None, param_name)]
                found_key = (mech_name, None, param_name)
                warning_txt = ("No Parameter found with "
                    f"param_name={param_name} and part_id={part_id} and "
                    f"mechanism={repr(mechanism)}. Parameter found under "
                    f"the key (mechanism.name, "
                    f"param_name)=({mech_name}, {param_name})")

            elif (mech_type, None, param_name) in self.parameters:
                found_entry = self.parameters[(mech_type, None, param_name)]
                found_key = (mech_type, None, param_name)
                warning_txt = ("No Parameter found with "
                    f"param_name={param_name} and part_id={part_id} and "
                    f"mechanism={repr(mechanism)}. Parameter found under "
                    f"the key (mechanism.name, "
                    f"param_name)=({mech_type}, {param_name})")
        # Finally try (param_name) --> return val
        if (None, None, param_name) in self.parameters and found_entry is None:
            found_entry = self.parameters[None, None, param_name]
            found_key = (None, None, param_name)
            warning_txt = (f"No parameter found with "
                f"param_name={param_name} and part_id={part_id} "
                f"and mechanism={repr(mechanism)}. Parameter "
                f"found under the key param_name={param_name}")

        if warning_txt is not None and parameter_warnings and found_entry is not None:
            warn(warning_txt)

        if found_entry is None:
            return None
        else:
            return_param = ModelParameter(found_entry.parameter_name,found_entry.value, (mechanism, part_id, param_name), found_key,
                parameter_keys = found_entry.keys, parameter_info = found_entry.info)
            return return_param

