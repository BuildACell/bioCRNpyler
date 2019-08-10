#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import patch, mock_open


class TestParameter(TestCase):

    def test_parameter_initialization(self):
        from biocrnpyler import Parameter

        with self.assertRaises(ValueError):
            Parameter(name='test_param', param_type='Nothing', value=0)

        with self.assertRaises(AssertionError):
            Parameter(name='test_param', param_type=[], value=0)

    def test__get_field_names(self):
        from biocrnpyler import Parameter

        with self.assertRaises(AssertionError):
            Parameter._get_field_names(field_names=None, accepted_field_names=None)
        with self.assertRaises(AssertionError):
            Parameter._get_field_names(field_names={}, accepted_field_names={})

        accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id'],
                                'param_name': ["parameter_name", "parameter", "param", "param_name"],
                                'part_id': ['part_id', 'part'],
                                'param_val': ["val", "value", "param_val", "parameter_value"]
                                }

        ret_dict = Parameter._get_field_names(field_names=[''], accepted_field_names=accepted_field_names)
        self.assertEqual(accepted_field_names.keys(),ret_dict.keys())

        field_names = ['part_id']

        with self.assertWarns(Warning):
            Parameter._get_field_names(field_names, accepted_field_names)

        accepted_field_names = {'dummy': ['dumb', 'dumber'],
                                }

        with self.assertWarns(Warning):
            Parameter._get_field_names(field_names, accepted_field_names)

    def test_load_parameter_file(self):
        from biocrnpyler import Parameter

        with self.assertRaises(AssertionError):
            Parameter.load_parameter_file(filename=None)

        # do NOT reformat this string
        example_file = """mechanism_id	part_id	param_name	param_val	comments
transcription_mm	ptet_tetR	kb	10.	extra columns are okay!
transcription_mm	ptet_tetR	ku	.1	These are the parameters for transcription"""

        with patch('builtins.open', mock_open(read_data=example_file), create=True):
            rtn_dict = Parameter.load_parameter_file(filename='test_file')

            right_dict = {('transcription_mm', 'ptet_tetR', 'kb'): 10.0, ('transcription_mm', 'ptet_tetR', 'ku'): 0.1}

            self.assertEqual(rtn_dict,right_dict)

    def test_create_parameter_dictionary(self):
        from biocrnpyler import Parameter

        parameters = Parameter.create_parameter_dictionary(parameters=None, parameter_file=None)

        self.assertEqual(parameters, None)

        with self.assertRaises(FileNotFoundError):
            Parameter.create_parameter_dictionary(parameters={}, parameter_file='dummy_file')
