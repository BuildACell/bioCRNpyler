# parameter.py - parameter processing
# RMM, 19 Aug 2018
#
# This file contains the Parameter class that is used for representing
# parameters, as well as utility functions for manipulating
# parameters.
#
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
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
    Initial conditions can be customized in through the custom_initial_conditions keyword in the Mixture constructor.
    custom_initial_conditions will take precedent to parameter initial conditions.

#### Units are read directly read from the column labeled "units" in the parameter file. 

"""

import csv
import numbers
import re
from collections import namedtuple  # Used for the parameter keys
from typing import Dict, List, Union
from warnings import warn

ParameterKey = namedtuple('ParameterKey', 'mechanism part_id name')  # This could later be extended


class Parameter(object):
    def __init__(self, parameter_name: str, parameter_value: Union[str, numbers.Real], unit=None):
        """A class for representing parameters in general. Only the below subclasses are ever used.

        :param parameter_name: is the name of the parameter
        :param parameter_value: is the value of the parameter
        :param unit: is the unit of the parameter or a species
        """
        self.parameter_name = parameter_name
        self.value = parameter_value
        self.unit = unit

    @property
    def parameter_name(self) -> str:
        return self._parameter_name

    @parameter_name.setter
    def parameter_name(self, new_parameter_name: str):
        if not isinstance(new_parameter_name, str):
            raise ValueError(f"parameter_name must be a string: received {type(new_parameter_name)}.")
        if not re.search('^[a-z]+', new_parameter_name, re.IGNORECASE):
            raise ValueError(f'parameter_name should be at least one character and cannot start with a number!')

        self._parameter_name = new_parameter_name

    @property
    def value(self) -> numbers.Real:
        return self._value

    @value.setter
    def value(self, new_parameter_value: Union[str, numbers.Real]):
        if not (isinstance(new_parameter_value, numbers.Real) or isinstance(new_parameter_value, str)):
            raise ValueError(f"parameter_value must be a float or int: received {type(new_parameter_value)}.")
        if isinstance(new_parameter_value, str):
            if re.search('[a-df-z]', new_parameter_value, re.I) \
                    or re.search('(^[1-9]+/[1-9]+)|(^[1-9]+e-?[0-9]+)|(^.?[0-9])', new_parameter_value, re.I) is None:
                raise ValueError(f'No valid parameter value! Accepted formats: 1.00 or 1e4 or 2/5, we got {new_parameter_value} ')

            self._value = Parameter._convert_rational(new_parameter_value)
        else:
            self._value = new_parameter_value

    @property
    def unit(self) -> str:
        return self._unit
    
    @unit.setter
    def unit(self, new_unit: str):
        if new_unit is None:
            self._unit = ""
        elif type(new_unit) is not str:
            raise ValueError(f"All units must be strings. Recieved {new_unit}.")
        else:
            self._unit = new_unit

    @staticmethod
    def _convert_rational(p_value: str) -> numbers.Real:
        if '/' in p_value:
            nom, denom = p_value.split('/')
            return float(nom)/float(denom)
        else:
            return float(p_value)
    def __eq__(self, other):
        if isinstance(other,Parameter):
            return self.value == other.value
        else:
            try:
                return self.value == float(other)
            except TypeError:
                raise TypeError(f"Cannot compare parameter {self} with {other}.")
        
    def __str__(self):
        return f"Parameter {self.parameter_name} = {self.value}"
    def __hash__(self):
        return hash(self._parameter_name)+hash(self._value)+hash(self._unit)


class ParameterEntry(Parameter):
    """A class for representing parameters in a parameter stored the ParameterDatabase.

     parameter_keys is a dictionary {key:value} or named_tuple (type ParameterKey) of keys for looking up the parameter
     parameter_info is a dictionary {key:value} of additional information about the parameter. 
         For example: additional columns in the parameter file or the parameter file name.
    """
    def __init__(self, parameter_name: str, parameter_value: Union[str,numbers.Real], parameter_key=None, parameter_info=None, **keywords):
        Parameter.__init__(self, parameter_name, parameter_value, **keywords)

        self.parameter_key = parameter_key
        self.parameter_info = parameter_info

    # Helper function to create ParameterKeys
    @staticmethod
    def create_parameter_key(new_key: Union[Dict, ParameterKey, str], parameter_name=None) -> ParameterKey:
        # New Key can be a named_tuple
        if isinstance(new_key, dict):
            new_key = dict(new_key)
            if parameter_name is not None:
                new_key["name"] = parameter_name
            for k in ParameterKey._fields:
                if k not in new_key:
                    new_key[k] = None
            return ParameterKey(**new_key) #automatically unpack the keywords
        elif isinstance(new_key, ParameterKey):
            return new_key
        elif isinstance(new_key, tuple) and len(list(new_key)) == len(ParameterKey._fields):
            # make a dictionary assuming correct ordering
            keywords = {ParameterKey._fields[i]:new_key[i] for i in range(len(ParameterKey._fields))}
            return ParameterKey(**keywords) #automatically unpack the keywords
        elif isinstance(new_key, str):
            return ParameterKey(mechanism = None, part_id = None, name = new_key)
        elif new_key is None and parameter_name is not None:
            return ParameterKey(mechanism = None, part_id = None, name = parameter_name)
        else:
            raise ValueError(f"parameter_key must be None, a dictionary, a ParameterKey, a {len(ParameterKey._fields)}-tuple, or a string (parameter name): received {new_key}.")

    @property
    def parameter_key(self) -> ParameterKey:
        return self._parameter_key

    @parameter_key.setter
    def parameter_key(self, parameter_key: Union[Dict, ParameterKey, str]):
        self._parameter_key = self.create_parameter_key(parameter_key, self.parameter_name)

    @property
    def parameter_info(self) -> Dict:
        return self._parameter_info

    @parameter_info.setter
    def parameter_info(self, parameter_info: Dict):
        if parameter_info is None:
            self._parameter_info = {}
        elif isinstance(parameter_info, dict):
            self._parameter_info = dict(parameter_info)

            #Update the units attribute, if necessary
            if "unit" in parameter_info and self.unit != "" and self.unit != parameter_info["unit"]:
                raise ValueError(f"Recieved multiple parameter units through constructor {self.unit} and parameter_info dictionary {parameter_info['unit']}.")
            elif "unit" in parameter_info:
                self.unit = parameter_info["unit"]

        else:
            raise ValueError(f"parameter_info must be None or a dictionary: received {parameter_info}.")

    def get_sbml_id(self):
        sbml_id = self.parameter_key.name+"_"
        if self.parameter_key.part_id is not None:
            sbml_id += self.parameter_key.part_id
        sbml_id += "_"
        if self.parameter_key.mechanism is not None:
            sbml_id += self.parameter_key.mechanism
        return sbml_id

    def __str__(self):
        return f"ParameterEntry({self.parameter_key}) = {self.value}"


class ModelParameter(ParameterEntry):
    """A class for representing parameters used in the Model.

      search_key is a tuple searched for to find the parameter, eg (mech_id, part_id, param_name), :
      found_key is the tuple used after defaulting to find the parameter eg (param_name)
    """
    def __init__(self, parameter_name: str, parameter_value: Union[str, numbers.Real], search_key, found_key, unit=None, parameter_key=None, parameter_info=None, **keywords):

        ParameterEntry.__init__(self, parameter_name, parameter_value, unit=unit, parameter_key=parameter_key, parameter_info=parameter_info, **keywords)
        self.search_key = search_key
        self.found_key = found_key

    @property
    def search_key(self):
        return self._search_key

    @search_key.setter
    def search_key(self, search_key):
        self._search_key = self.create_parameter_key(search_key, self.parameter_name)

    @property
    def found_key(self):
        return self._found_key

    @found_key.setter
    def found_key(self, found_key):
        self._found_key = self.create_parameter_key(found_key, self.parameter_name)

    def __str__(self):
        return f"ModelParameter({self.parameter_key}) = {self.value}\tsearch_key={self.search_key}"


class ParameterDatabase(object):
    def __init__(self, parameter_dictionary=None, parameter_file=None, overwrite_parameters=False):
        """A class for storing parameters in Components and Mixtures.

        :param parameter_dictionary:
        :param parameter_file:
        :param overwrite_parameters: whether to overwrite existing entries in the parameter database
        """

        self.parameters = {} #create an emtpy dictionary to get parameters.

        if isinstance(parameter_file, str):
            self.load_parameters_from_file(parameter_file, overwrite_parameters = overwrite_parameters)
        elif isinstance(parameter_file, list):
            for p in parameter_file:
                if isinstance(p, str):
                    self.load_parameters_from_file(p, overwrite_parameters = overwrite_parameters)
                else:
                    raise ValueError("parameter_file must be a string or list of strings representing file names and paths.")
        elif parameter_file is not None:
            raise ValueError("parameter_file must be a string representing a file name and path.")
            
        if isinstance(parameter_dictionary, dict):
            self.load_parameters_from_dictionary(parameter_dictionary, overwrite_parameters = overwrite_parameters)
        elif parameter_dictionary is not None:
            raise ValueError("parameter_dictionary must be None or a dictionary!")

    # To check if a key or ParameterEntry is in a the ParameterDatabase
    def __contains__(self, val):
        if isinstance(val, ParameterEntry):
            key = val.parameter_key
            if key in self.parameters and self.parameters[key] == val:
                return True
            else:
                return False
        else:
            try:
                key = ParameterEntry.create_parameter_key(val)
                return key in self.parameters
            except ValueError:
                return False

    # Ability to loop through parameters eg
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

    # Length method
    def __len__(self):
        return len(self.parameters)

    # Gets a parameter from the database
    # Only returns exact matches.
    def __getitem__(self, key):
        param_key = ParameterEntry.create_parameter_key(key)
        return self.parameters[param_key]

    # Sets a parameter in the databases - useful for quickly changing parameters, but add_parameter is recommended.
    def __setitem__(self, parameter_key, value):

        key = ParameterEntry.create_parameter_key(parameter_key)
        
        if isinstance(value, ParameterEntry):
            if key != value.parameter_key:
                raise ValueError(f"Parameter Key does not match: ParameterDatabase key {key} is not the same as ParameterEntry Key {value.parameter_key}.")
            self.parameters[key] = value
        else:
            self.add_parameter(key.name, value, parameter_key = key, parameter_origin = "Set Manually", overwrite_parameters = True)

    def __str__(self):
        txt = "ParameterDatabase:"
        param_txt = "\n".join([repr(p) for p in self.parameters])
        return txt+param_txt

    def add_parameter(self, parameter_name: str, parameter_value: Union[str,numbers.Real], parameter_origin = None, parameter_key = None, parameter_info = None, overwrite_parameters = False):
        """Adds a parameter to the database with appropriate metadata

        :param parameter_name: the name of the parameter
        :param parameter_value: the value of the parameter
        :param parameter_origin:
        :param parameter_key:
        :param parameter_info:
        :param overwrite_parameters: whether to overwrite existing entries in the parameter database
        :return:
        """

        # Put parameter origin into parameter_info
        if parameter_info is None:
            parameter_info = {}
        if "parameter origin" not in parameter_info:
            parameter_info["parameter origin"] = parameter_origin

        # Create ParameterEntry
        param = ParameterEntry(parameter_name, parameter_value, parameter_key = parameter_key, parameter_info = parameter_info)
        key = param.parameter_key
        
        # Update parameter dictionary
        if key in self.parameters and not overwrite_parameters:
            raise ValueError(f"Duplicate parameter detected. Parameter with key = {key} is already in the ParameterDatabase. To Overwrite existing parameters, use overwrite_parameters = True.")
        else:
            self.parameters[key] = param

    def load_parameters_from_dictionary(self, parameter_dictionary: Dict[ParameterKey, Union[str,numbers.Real]], overwrite_parameters=False) -> None:
        """Loads Parameters from a parameter dictionary.

        :param parameter_dictionary: Dictionary with keys ParameterKey types and values with real numbers
        :param overwrite_parameters: whether to overwrite existing entries in the parameter database
        """
        for k in parameter_dictionary:
            key = ParameterEntry.create_parameter_key(k)
            self.add_parameter(key.name, parameter_dictionary[k], parameter_key = {"part_id":key.part_id, "mechanism":key.mechanism}, parameter_origin = "parameter_dictionary", overwrite_parameters = overwrite_parameters)

    def load_parameters_from_database(self, parameter_database, overwrite_parameters=False) -> None:
        """Loads parameters from another ParameterDatabase.

        :param parameter_database: instance of another ParameterDatabase
        :param overwrite_parameters:  whether to overwrite existing entries in the parameter database
        """

        if not isinstance(parameter_database, ParameterDatabase):
            raise TypeError(f"paramater_database must be a ParamaterDatabase: recievied {parameter_database}.")

        for k in parameter_database:
            if k not in self.parameters or overwrite_parameters:
                self.parameters[k.parameter_key] = parameter_database[k.parameter_key]
            else:
                raise ValueError(f"Duplicate parameter detected. Parameter with key = {k} is already in the ParameterDatabase. To Overwrite existing parameters, use overwrite_parameters = True.")

    def load_parameters_from_file(self, filename: str, overwrite_parameters=False) -> None:
        """Loads parameters from a file to the ParameterDatabase.

        Parameter files must be tab-separated (.tsv or .txt) or comma-separated (.csv) files!
        :param filename: name of the file (with valid file path)
        :param overwrite_parameters: whether to overwrite existing entries in the parameter database
        """

        # Figure out the format of the parameter file from the file extension
        with open(filename) as f:
            file_type = filename.split(".")[-1]
            if file_type in ["tsv", "txt"]:
                delimiter = '\t'
            elif file_type in ["csv"]:
                delimiter = ","
            else:
                raise ValueError("Parameter files must be tab-seperated (.tsv or .txt) or comma-seperated (.csv) files.")

            csvreader = csv.DictReader(f, delimiter=delimiter)
            # Used for flexible column headings
            accepted_field_names = {
                'mechanism': ['mechanism', 'mechanism_id'],
                'param_name': ["parameter_name", "parameter", "param", "param_name"],
                'part_id': ['part_id', 'part'],
                'param_val': ["val", "value", "param_val", "parameter_value"],
                'unit':['unit', 'units']
            }

            field_names = self._get_field_names(csvreader.fieldnames, accepted_field_names)

            # Determine which columns are in the CSV
            if field_names['param_name'] is None:
                warn('No param_name column was found, could not load parameter!')
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
                param_value = row[field_names['param_val']]
                field_columns = [field_names['param_name'], field_names['part_id'], field_names['mechanism'], field_names['param_val']]
                parameter_info = {k:row[k] for k in row if k not in field_columns}
                # TODO test all these cases!

                # Case 1: No Param Name so skip the row
                if row[field_names['param_name']] is None or len(row[field_names['param_name']]) == 0:
                    pass

                # Case 2: Just a Param Name
                elif no_mechism_column and no_part_id_column:
                    param_name = row[field_names['param_name']]
                    self.add_parameter(param_name, param_value, parameter_origin = filename, 
                        parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                # Case 3: Part_id and Param Name
                elif no_mechism_column and no_part_id_column is False:
                    param_name = row[field_names['param_name']]
                    part_id = row[field_names['part_id']]

                    if part_id is not None and len(part_id) > 0:
                        self.add_parameter(param_name, param_value, parameter_key = {"part_id":part_id}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info , overwrite_parameters = overwrite_parameters)

                # Case 4: mechanism and param name
                elif no_part_id_column and no_mechism_column is False:
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if mech_name is not None and len(mech_name) > 0:
                        self.add_parameter(param_name, param_value, parameter_key = {"mechanism":mech_name}, parameter_origin = filename, 
                        parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                # Case 5: mechanism, part_id, and param name
                else:
                    part_id = row[field_names['part_id']]
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if part_id is not None and len(part_id) > 0 and mech_name is not None and len(mech_name) >0:
                        self.add_parameter(param_name, param_value, parameter_key = {"part_id":part_id, "mechanism":mech_name}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    
                    elif part_id is not None and len(part_id) > 0:    
                        self.add_parameter(param_name, param_value, parameter_key = {"part_id":part_id}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

                    elif mech_name is not None and len(mech_name) >0:
                        self.add_parameter(param_name, param_value, parameter_key = {"mechanism":mech_name}, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)
                    else:
                        self.add_parameter(param_name, param_value, parameter_origin = filename, 
                            parameter_info = parameter_info, overwrite_parameters = overwrite_parameters)

    @staticmethod
    def _get_field_names(field_names: List[str], accepted_field_names: Dict[str, List[str]]) -> Dict[str, str]:
        """Searches through valid field names and finds the currently used one. It builds a dictionary of currently
            used field names.

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

    def find_parameter(self, mechanism, part_id, param_name):
        """Searches the database for the best matching parameter. 
        
        Parameter defaulting hierarchy:
        (mechanism_name, part_id, param_name) --> param_val. If that particular parameter key cannot be found, 
        the software will default to the following keys: 
        (mechanism_type, part_id, param_name) >> (part_id, param_name) >> 
        (mechanism_name, param_name) >> (mechanism_type, param_name) >>
        (param_name) and give a warning. 
        As a note, mechanism_name refers to the .name variable of a Mechanism. mechanism_type refers to the .type variable of a Mechanism. 
        Either of these can be used as a mechanism_id. This allows for models to be constructed easily using default parameter values and 
        for parameters to be shared between different Mechanisms and/or Components.
        """

        #this is imported here because otherwise there are import loops
        from .mechanism import Mechanism

        found_entry = None

        if isinstance(mechanism, str):
            mech_name = mechanism
            mech_type = mechanism
        elif isinstance(mechanism, Mechanism):
            mech_name = mechanism.name
            mech_type = mechanism.mechanism_type
        elif mechanism is not None:
            raise ValueError(f"mechanism keyword must be or string or have name and mechanism_type attributes: recievied {mechanism}.")
        else:
            mech_name = None
            mech_type = None
        
        parameter_key_list = [
            ParameterKey(mechanism = mech_name, part_id = part_id, name = param_name), 
            ParameterKey(mechanism = mech_type, part_id = part_id, name = param_name),
            ParameterKey(mechanism = None, part_id = part_id, name = param_name),
            ParameterKey(mechanism = mech_name, part_id = None, name = param_name),
            ParameterKey(mechanism = mech_type, part_id = None, name = param_name),
            ParameterKey(mechanism = None, part_id = None, name = param_name),
            ]

        for key in parameter_key_list:
            if key in self.parameters and found_entry is None:
                found_entry = self.parameters[key]
                found_key = key
                break

        if found_entry is None:
            return None
        else:
            return_param = ModelParameter(found_entry.parameter_name,found_entry.value, (mech_name, part_id, param_name), found_key,
                parameter_key = found_entry.parameter_key, parameter_info = found_entry.parameter_info, unit = found_entry.unit)
            return return_param
