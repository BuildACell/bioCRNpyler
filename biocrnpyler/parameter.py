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

Please see the Parameter ipynb in examples for details on mechanism parameters.

#### Parameter File Overview:
Parameter files can be TSVs or CSVs. The first line of the file should contain column headings. 

The following headings are required (in any order): mechanism_id, part_id, param_name, param_val 
    (spaces can be substituted for underscores and headings are not case sensitive).

mechanism_id is the name of the Mechanism or the kind of mechanism that will use this parameter, 
    for example "transcription" or  "transcription_mm" for Mechalis-Menten transcription would go in this column. 
part_id refers to the name of the Component that will use this mechanism, 
    for example "ptet" for a tet repressed promoter. 
param_name refers to the name of the model parameter, 
    for example "ktx", "kb", or "ku". The value of these columns is case sensitive and underscores are different from spaces.

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

#### Multiple Parameter Files:
Components and Mixtures can both have one more multiple parameter files by passing in a list of filenames 
    instead of a single filename to the parameter_file keyword. 
    Components use parameters loaded from their file(s) before defaulting to the file(s) supplied to a Mixture. 
    The last file in any list will take precedent and overwrite parameter files which were written earlier.

#### Suppressing warnings
To suppress parameter warnings, use the keyword parameter_warnings = False inside a Mixture or Component constructor.

Below is an example csv with all the parameters for a tetR promoter undergoing Michalis Menten transcription and translation.

"""

import csv
from warnings import warn
from typing import List, Dict, Union


class Parameter(object):
    """Parameter value (reaction rates)"""

    def __init__(self, name, param_type, value, comment="", debug=False):
        if not isinstance(param_type, str):
            raise ValueError('parameter_type must be a string')

        self.name = name.strip()
        self.param_type = param_type.strip()
        self.comment = comment.strip()

        # Set the value of the parameter
        if debug:
            print("%s [%s] = %s" % (self.name, self.param_type, value))
        if param_type.strip() == 'Numeric':
            self.value = float(value)  # store as float
        elif param_type.strip() == 'Expression':
            self.value = value  # store as string
        else:
            raise ValueError("can't parse value of parameter %s" % name)

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

    @staticmethod
    def load_parameter_file(filename: str) -> Dict:
        """load the parameter configuration file into a dictionary
        :param filename: valid parameter file name
        :return: a dictionary of the parameters
        """
        assert isinstance(filename, str) and len(filename) > 0
        param_dict = {}
        # TODO implement search through possible parameter config file locations
        # Open up the CSV file for reaching
        with open(filename) as f:
            csvreader = csv.DictReader(f, delimiter='\t')

            accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id'],
                                    'param_name': ["parameter_name", "parameter", "param", "param_name"],
                                    'part_id': ['part_id', 'part'],
                                    'param_val': ["val", "value", "param_val", "parameter_value"]
                                    }

            field_names = Parameter._get_field_names(csvreader.fieldnames, accepted_field_names)

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

            for row in csvreader:
                # TODO what about integers? float might cause numerical drift in simulations, e.g. cooperativity=2.001
                param_value = float(row[field_names['param_val']])
                # TODO test all these cases!
                if row[field_names['param_name']] is None or len(row[field_names['param_name']]) == 0:
                    pass
                elif no_mechism_column and no_part_id_column:
                    param_name = row[field_names['param_name']]
                    param_dict[param_name] = param_value
                elif no_mechism_column and no_part_id_column is False:
                    if row[field_names['part_id']] is not None and len(row[field_names['part_id']]) > 0:
                        part_id = row[field_names['part_id']]
                        param_name = row[field_names['param_name']]
                        param_dict[(part_id, param_name)] = param_value
                    else:
                        param_name = row[field_names['param_name']]
                        param_dict[param_name] = param_value
                elif no_part_id_column and no_mechism_column is False:
                    if row[field_names['mechanism']] is not None and len(row[field_names['mechanism']]) > 0:
                        mech_name = row[field_names['mechanism']]
                        param_name = row[field_names['param_name']]
                        param_dict[(mech_name, param_name)] = param_value
                    else:
                        param_name = row[field_names['param_name']]
                        param_dict[param_name] = param_value
                else:
                    if row[field_names['part_id']] is not None and len(row[field_names['part_id']]) > 0:
                        if row[field_names['mechanism']] is not None and len(row[field_names['mechanism']]) > 0:
                            part_id = row[field_names['part_id']]
                            mech_name = row[field_names['mechanism']]
                            param_name = row[field_names['param_name']]
                            param_dict[(mech_name, part_id, param_name)] = param_value
                        else:
                            part_id = row[field_names['part_id']]
                            param_name = row[field_names['param_name']]
                            param_dict[(part_id, param_name)] = param_value
                    else:
                        if row[field_names['mechanism']] is not None and len(row[field_names['mechanism']]) > 0:
                            mech_name = row[field_names['mechanism']]
                            param_name = row[field_names['param_name']]
                            param_dict[(mech_name, param_name)] = param_value
                        else:
                            param_name = row[field_names['param_name']]
                            param_dict[param_name] = param_value
        return param_dict

    @staticmethod
    def create_parameter_dictionary(parameters: Union[None, Dict],
                                    parameter_file: Union[None, str, List[str]]) -> Union[None, Dict]:
        """
        Loads parameter config file(s) and merges the parameters with the existing parameters dictionary
        :param parameters: existing parameters dictionary or None
        :param parameter_file: valid parameter file(s)
        :return: updated parameters dictionary (empty dict if no parameter file was given)
        """

        # no parameter dictionary was given creating one
        if not isinstance(parameters, dict):
            parameters = {}

        # empty parameter_file no new parameters are loaded
        if parameter_file is None:
            return parameters

        assert isinstance(parameter_file, str) or isinstance(parameter_file, list)

        if isinstance(parameter_file, list):
            file_list = parameter_file
        else:
            file_list = [parameter_file]

        for file_name in file_list:
            new_parameters = Parameter.load_parameter_file(file_name)
            parameters.update(new_parameters)

        return parameters
