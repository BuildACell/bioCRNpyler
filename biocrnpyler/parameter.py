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
   strength) across multiple parameters.  Extract parameter values are
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
from warnings import warn
from typing import List, Dict, Union


class Parameter(object):
    """Parameter value (reaction rates)"""

    def __init__(self, name, param_type, value, comment="", debug=False):
        assert type(param_type) is str

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
        assert isinstance(field_names, list) and len(field_names) > 0
        assert isinstance(accepted_field_names, dict) and len(accepted_field_names) > 0

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
        :param parameters: existing parameters dictionary
        :param parameter_file: valid parameter file(s)
        :return: updated parameters dictionary
        """
        # empty call no parameters are loaded
        if parameters is None or parameter_file is None:
            return parameters

        assert isinstance(parameters, dict)
        assert isinstance(parameter_file, str) or isinstance(parameter_file, list)

        if isinstance(parameter_file, list):
            file_list = parameter_file
        else:
            file_list = [parameter_file]

        for file_name in file_list:
            new_parameters = Parameter.load_parameter_file(file_name)
            parameters.update(new_parameters)

        return parameters
