#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import patch, mock_open
from biocrnpyler import Parameter
import sys
from warnings import warn


class TestParameter(TestCase):

    def test_parameter_initialization(self):
        # test unknown parameter type
        with self.assertRaisesRegexp(ValueError,"can't parse value of parameter"):
            Parameter(name='test_param', param_type='Nothing', value=0)
        # test invalid parameter type
        with self.assertRaisesRegexp(ValueError,'parameter_type must be a string'):
            Parameter(name='test_param', param_type=[], value=0)

    def test__get_field_names(self):
        test_accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id']}
        valid_field_names = ['part_id']

        # test None as field_names
        with self.assertRaisesRegexp(ValueError, 'field_names must be a list of strings'):
            Parameter._get_field_names(field_names=None, accepted_field_names=test_accepted_field_names)
        # test invalid field_names type
        with self.assertRaisesRegexp(ValueError, 'field_names must be a list of strings'):
            Parameter._get_field_names(field_names={}, accepted_field_names=test_accepted_field_names)
        # test empty field_names list
        with self.assertRaisesRegexp(ValueError,'field_names cannot be empty list!'):
            Parameter._get_field_names(field_names=[], accepted_field_names=test_accepted_field_names)

        # test None as accepted_field_names
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names must be a dictionary'):
            Parameter._get_field_names(field_names=valid_field_names, accepted_field_names=None)
        # test invalid accepted_field_names type
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names must be a dictionary'):
            Parameter._get_field_names(field_names=valid_field_names, accepted_field_names=[])
        # test empty field_names list
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names cannot be empty dictionary'):
            Parameter._get_field_names(field_names=valid_field_names, accepted_field_names={})

        accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id'],
                                'param_name': ["parameter_name", "parameter", "param", "param_name"],
                                'part_id': ['part_id', 'part'],
                                'param_val': ["val", "value", "param_val", "parameter_value"]
                                }

        ret_dict = Parameter._get_field_names(field_names=[''], accepted_field_names=accepted_field_names)
        self.assertEqual(accepted_field_names.keys(), ret_dict.keys())

        with self.assertWarns(Warning):
            Parameter._get_field_names(valid_field_names, accepted_field_names)

        accepted_field_names = {'dummy': ['dumb', 'dumber'],
                                }

        with self.assertWarns(Warning):
            Parameter._get_field_names(valid_field_names, accepted_field_names)

    def test_load_parameter_file(self):

        with self.assertRaises(AssertionError):
            Parameter.load_parameter_file(filename=None)

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            # do NOT reformat this string below
            example_csv = """mechanism_id	part_id	param_name	param_val	comments\ntranscription_mm	ptet_tetR	kb	10.	extra columns are okay!\ntranscription_mm	ptet_tetR	ku	.1	These are the parameters for transcription"""

            with patch('builtins.open', mock_open(read_data=example_csv), create=True):
                rtn_dict = Parameter.load_parameter_file(filename='test_file')

                right_dict = {('transcription_mm', 'ptet_tetR', 'kb'): 10.0, ('transcription_mm', 'ptet_tetR', 'ku'): 0.1}

                self.assertEqual(rtn_dict, right_dict)
        else:
            warn('version below 3.6 was detected! This test was skipped')

    def test_create_parameter_dictionary(self):

        # test that no parameter dictionary or parameter files are given empty parameter dict is returned
        empty_dict = {}
        parameters = Parameter.create_parameter_dictionary(parameters=None, parameter_file=None)
        self.assertEqual(parameters, empty_dict)

        # test that no parameter file is given then the supplied parameter dictionary is returned
        param_dict = {'kb': 10.0}
        parameters = Parameter.create_parameter_dictionary(parameters=param_dict, parameter_file=None)
        self.assertEqual(parameters, param_dict)

        with self.assertRaises(FileNotFoundError):
            Parameter.create_parameter_dictionary(parameters={}, parameter_file='dummy_file')
