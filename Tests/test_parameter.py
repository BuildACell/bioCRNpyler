#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import patch, mock_open
from biocrnpyler import Parameter, ParameterEntry, ModelParameter, ParameterDatabase
import sys
from warnings import warn


class TestParameter(TestCase):

    def test_parameter(self):
        # test parameter name
        with self.assertRaisesRegexp(ValueError, f"parameter_name must be a string"):
            Parameter(parameter_name=None, parameter_value=1.0)
        # test parameter value
        with self.assertRaisesRegexp(ValueError, f"parameter_value must be a float or int"):
            Parameter(parameter_name="None", parameter_value=None)

        #Test string parameter values
        self.assertTrue(Parameter(parameter_name="None", parameter_value="1.0").value == 1.0)

    def test_parameter_entry(self):
        #Valid ParameterEntry Construction
        ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_keys = {"part_id":"id"}, parameter_info = {"comment":"comment"})

        #Invalid keys
        param_keys = (1, 2, 3)
        with self.assertRaisesRegexp(ValueError, "parameter_keys must be a None or a dictionary"):
            ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_keys = param_keys)

        #Invalid info
        param_info = "blah blah"
        with self.assertRaisesRegexp(ValueError, f"parameter_info must be None or a dictionary"):
            ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_info = param_info)

    def test_model_parameter(self):
        #valid ModelParameter Construction
        ModelParameter(parameter_name="None", parameter_value=1.0, search_key = ("that", "this", "k"), found_key = ("this", None, "k"))

        #Invalid keys
        with self.assertRaisesRegexp(ValueError,"search_key must be a tuple"):
            ModelParameter(parameter_name="None", parameter_value=1.0, search_key = "k", found_key = ("this", None, "k"))

        with self.assertRaisesRegexp(ValueError,"found_key must be a tuple"):
            ModelParameter(parameter_name="None", parameter_value=1.0, search_key = ("that", "this", "k"), found_key = "k")


    def test_get_field_names(self):
        PD = ParameterDatabase
        test_accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id']}
        valid_field_names = ['part_id']

        # test None as field_names
        with self.assertRaisesRegexp(ValueError, 'field_names must be a list of strings'):
            PD._get_field_names(field_names=None, accepted_field_names=test_accepted_field_names)
        # test invalid field_names type
        with self.assertRaisesRegexp(ValueError, 'field_names must be a list of strings'):
            PD._get_field_names(field_names={}, accepted_field_names=test_accepted_field_names)
        # test empty field_names list
        with self.assertRaisesRegexp(ValueError,'field_names cannot be empty list!'):
            PD._get_field_names(field_names=[], accepted_field_names=test_accepted_field_names)

        # test None as accepted_field_names
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names must be a dictionary'):
            PD._get_field_names(field_names=valid_field_names, accepted_field_names=None)
        # test invalid accepted_field_names type
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names must be a dictionary'):
            PD._get_field_names(field_names=valid_field_names, accepted_field_names=[])
        # test empty field_names list
        with self.assertRaisesRegexp(ValueError, 'accepted_field_names cannot be empty dictionary'):
            PD._get_field_names(field_names=valid_field_names, accepted_field_names={})

        accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id'],
                                'param_name': ["parameter_name", "parameter", "param", "param_name"],
                                'part_id': ['part_id', 'part'],
                                'param_val': ["val", "value", "param_val", "parameter_value"]
                                }

        ret_dict = PD._get_field_names(field_names=[''], accepted_field_names=accepted_field_names)
        self.assertEqual(accepted_field_names.keys(), ret_dict.keys())

        with self.assertWarns(Warning):
            PD._get_field_names(valid_field_names, accepted_field_names)

        accepted_field_names = {'dummy': ['dumb', 'dumber'],
                                }

        with self.assertWarns(Warning):
            PD._get_field_names(valid_field_names, accepted_field_names)

    def load_parameters_from_file(self):

        #Bad parameter file keyword
        with self.assertRaisesRegexp(ValueError, f"parameter_file must be a string representing a file name (and path)."):
            PD = ParameterDatabase(parameter_file = {})

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            # do NOT reformat this string below
            example_csv = """
            mechanism_id	part_id	param_name	param_val	comments\n
            transcription_mm	ptet_tetR	kb	10.	extra columns are okay!\n4
            transcription_mm	ptet_tetR	ku	.1	These are the parameters for transcription\n
                P kb    2.0 A parameter with no mechanism!\n
            mechanism   ku  100.    A parameter with no part_id!\n
                    kdefault    1.0   A default parameter!\n
                        3.3 A useless row!\n

            """

            with patch('builtins.open', mock_open(read_data=example_csv), create=True):
                PD = ParameterDatabase(parameter_file = 'test_file')

                right_dict = {
                ('transcription_mm', 'ptet_tetR', 'kb'): 10.0, 
                ('transcription_mm', 'ptet_tetR', 'ku'): 0.1,
                (None, 'P', 'kb'): 2.0,
                ("mechanism", None, "ku"): 100,
                (None, None, "kdefault"): 1.0
                }

                return_dict = {k:p.value for (k, p) in PD.parameters}
                self.assertEqual(return_dict, right_dict)
        else:
            warn('version below 3.6 was detected! This test was skipped')


        #Lets load an actual parameter files from the examples
        import os
        os.chdir('../')
        os.chdir('examples')
        PD = ParameterDatabase(parameter_file = "default_parameters.txt")
        self.assertTrue("ktx" in PD)
        self.assertTrue((None, "weak", "ktx") in PD)
        self.assertTrue(("simple_transcript", "weak", "ktx") in PD)


    def load_parameters_from_dictionary(self):

        #bad parameter_dictionary keyword
        with self.assertRaisesRegexp(ValueError, f"parameter_dictionary must be None or a dictionary!"):
            PD = ParameterDatabase(parameter_dictionary = 'test_file')

        #proper parameter dictionary
        parameter_dict = {
        "k":1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }
        PD = ParameterDatabase(parameter_dictionary = parameter_dict)
        return_dict = {k:p.value for (k, p) in PD.parameters}
        self.assertEqual(return_dict, right_dict)

        #improper parameter dictionary
        k = ("M", "k")
        parameter_dict = {k:2.0}
        with self.assertRaisesRegexp(ValueError, f"Invalid parameter key {k}. parameter_dictionary keys must be a 3-tuple or a single string: ('mechanism', 'part_id', 'param_name') OR 'param_name'. Note 'mechanism' and 'part_id' can be None."):
            PD = ParameterDatabase(parameter_dictionary = 'test_file')

        #duplicate parameter dictionary
        key = (None, None, "k")
        parameter_dict = {"k":1, (None, None, "k"):1}
        with self.assertRaisesRegexp(ValueError, f"Duplicate parameter detected. Parameter with key = {key} is already in the ParameterDatabase. To Overwrite existing parameters, use overwrite_parameters = True."):
            PD = ParameterDatabase(parameter_dictionary = 'test_file')

    def test_iterator(self):
        parameter_dict = {
        "k":1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }

        PD = ParameterDatabase(parameter_dictionary = parameter_dict)

        #All things in iterator are entries
        count = 0
        for entry in PD:
            count += 1
            self.assertTrue(isinstance(entry, ParameterEntry))

        #The correct number of entries
        self.assertEqual(count, len(parameter_dict))

    def test_contains(self):
        parameter_dict = {
        "k":1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }

        PD = ParameterDatabase(parameter_dictionary = parameter_dict)

        self.assertTrue("k" in PD)
        self.assertTrue(("M", "pid", "k") in PD)
        self.assertTrue(PD["k"] in PD)
        self.assertFalse("ktx" in PD)
        self.assertFalse(("M", "k") in PD)

    def test_indexing(self):
        parameter_dict = {
        "k":1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }

        PD = ParameterDatabase(parameter_dictionary = parameter_dict)

        #Test correct accessing
        self.assertTrue(PD["k"].value == 1)
        self.assertTrue(PD[("M", "pid", "k")].value == 2.0)
        self.assertTrue(PD[(None, "pid", "k")].value == 3.3)
        self.assertTrue(PD[("M", None, "k")].value == 4)

        #test incorrect accessing
        with self.assertRaisesRegexp(ValueError, f"Invalid parameter key"):
            PD[("M", "k")]

        #test accessing something not in the PD
        with self.assertRaisesRegexp(ValueError, f"Parameter"):
            PD["kb"]

        #test inserting values
        PD["kb"] = 100
        self.assertTrue(PD["kb"].value == 100)
        PD[(None, None, "ku")] = 200
        self.assertTrue(PD["ku"].value == 200)
        PD[("M", "pid", "ktx")] = 300
        self.assertTrue(PD[("M", "pid", "ktx")].value == 300)

        #Test inserting ParameterEntry
        PE = ParameterEntry("test", 1.0)
        PD["test"] = PE
        self.assertTrue(PE in PD)

        #Test Correct Overwriting
        PE = ParameterEntry("test", 2.0)
        PD["test"] = PE
        self.assertTrue(PD["test"].value == PE.value)

        #Test incorrect Overwriting
        PE = ParameterEntry("t", 1.0)
        with self.assertRaisesRegexp(ValueError, f"Parameter Key does not match"):
            PD["test"] = PE


        #Invalid parameter key
        with self.assertRaisesRegexp(ValueError, f"Invalid parameter"):
            PD[("M", "k")] = 10

        #Test overwriting
        PD["k"] = .1
        self.assertTrue(PD["k"].value == .1)
        PD[("M", "pid", "ktx")] = .333
        self.assertTrue(PD[("M", "pid", "ktx")].value == .333)


