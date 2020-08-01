#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import patch, mock_open
from biocrnpyler import Parameter, ParameterEntry, ModelParameter, ParameterDatabase, ParameterKey, Mechanism
import sys
from warnings import warn


class TestParameter(TestCase):

    def test_parameter(self):
        # test parameter name
        with self.assertRaisesRegex(ValueError, f"parameter_name must be a string"):
            Parameter(parameter_name=None, parameter_value=1.0)
        # test parameter value
        with self.assertRaisesRegex(ValueError, f"parameter_value must be a float or int"):
            Parameter(parameter_name="None", parameter_value=None)

        # test invalid value string
        with self.assertRaisesRegex(ValueError, f"No valid parameter value! Accepted format"):
            Parameter(parameter_name="None", parameter_value='2ba')

        # test string parameter values
        self.assertTrue(Parameter(parameter_name="None", parameter_value="1.0").value == 1.0)

        self.assertTrue(Parameter(parameter_name="None", parameter_value="1/2").value == 0.5)

        self.assertTrue(Parameter(parameter_name="None", parameter_value="1e2").value == 100)

        # testing invalid parameter name
        with self.assertRaisesRegex(ValueError, f"parameter_name should be at least one character and cannot start with a number!"):
            Parameter(parameter_name="2", parameter_value=2)

    def test_parameter_entry(self):
        #Valid ParameterEntry Construction
        ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_key = {"part_id":"id"}, parameter_info = {"comment":"comment"})

        #Invalid keys
        param_keys = Parameter(parameter_name="None", parameter_value=1.0)
        with self.assertRaisesRegex(ValueError, "parameter_key must be"):
            ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_key = param_keys)

        #Invalid info
        param_info = "blah blah"
        with self.assertRaisesRegex(ValueError, f"parameter_info must be None or a dictionary"):
            ParameterEntry(parameter_name="None", parameter_value=1.0, parameter_info = param_info)

    def test_model_parameter(self):
        #valid ModelParameter Construction
        ModelParameter(parameter_name="None", parameter_value=1.0, search_key = ("that", "this", "k"), found_key = ("this", None, "k"))

        #Invalid keys
        k = Parameter(parameter_name="None", parameter_value="1.0")
        with self.assertRaisesRegex(ValueError,"parameter_key must be None"):
            ModelParameter(parameter_name="None", parameter_value=1.0, search_key = k, found_key = ("this", None, "k"))

        with self.assertRaisesRegex(ValueError,"parameter_key must be None"):
            ModelParameter(parameter_name="None", parameter_value=1.0, search_key = ("that", "this", "k"), found_key = k)


    def test_get_field_names(self):
        PD = ParameterDatabase
        test_accepted_field_names = {'mechanism': ['mechanism', 'mechanism_id']}
        valid_field_names = ['part_id']

        # test None as field_names
        with self.assertRaisesRegex(ValueError, 'field_names must be a list of strings'):
            PD._get_field_names(field_names=None, accepted_field_names=test_accepted_field_names)
        # test invalid field_names type
        with self.assertRaisesRegex(ValueError, 'field_names must be a list of strings'):
            PD._get_field_names(field_names={}, accepted_field_names=test_accepted_field_names)
        # test empty field_names list
        with self.assertRaisesRegex(ValueError,'field_names cannot be empty list!'):
            PD._get_field_names(field_names=[], accepted_field_names=test_accepted_field_names)

        # test None as accepted_field_names
        with self.assertRaisesRegex(ValueError, 'accepted_field_names must be a dictionary'):
            PD._get_field_names(field_names=valid_field_names, accepted_field_names=None)
        # test invalid accepted_field_names type
        with self.assertRaisesRegex(ValueError, 'accepted_field_names must be a dictionary'):
            PD._get_field_names(field_names=valid_field_names, accepted_field_names=[])
        # test empty field_names list
        with self.assertRaisesRegex(ValueError, 'accepted_field_names cannot be empty dictionary'):
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

    def test_load_parameters_from_file(self):

        # Bad parameter file keyword
        with self.assertRaisesRegex(ValueError, f"parameter_file must be a string representing a file name and path."):
            ParameterDatabase(parameter_file={})

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            # !!! DO NOT reformat this string below !!!!
            example_csv = """mechanism_id	part_id	param_name	param_val	comments\ntranscription_mm	ptet_tetR	kb	10.	extra columns are okay!\ntranscription_mm	ptet_tetR	ku	.1	These are the parameters for transcription"""
            # !!! DO NOT reformat this string above !!!!

            with patch('builtins.open', mock_open(read_data=example_csv), create=True):
                PD = ParameterDatabase(parameter_file='test_file.tsv')

                
                right_dict = {
                ('transcription_mm', 'ptet_tetR', 'kb'): 10.0,
                ('transcription_mm', 'ptet_tetR', 'ku'): 0.1
                }
                returned_dict = {(k.mechanism, k.part_id, k.name):PD.parameters[k].value for k in PD.parameters}
                #raise ValueError(str(returned_dict))
                self.assertEqual(right_dict, returned_dict)
        else:
            warn('version below 3.6 was detected! This test was skipped')

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            example_csv = """mechanism_id"""

            with patch('builtins.open', mock_open(read_data=example_csv), create=True):
                with self.assertWarnsRegex(Warning, f"No param_name column was found, could not load parameter!"):
                    ParameterDatabase(parameter_file='test_file.tsv')

        else:
            warn('version below 3.6 was detected! This test was skipped')


    def test_load_parameters_from_dictionary(self):

        # bad parameter_dictionary keyword
        with self.assertRaisesRegex(ValueError, f"parameter_dictionary must be None or a dictionary!"):
            PD = ParameterDatabase(parameter_dictionary='test_file')

        # proper parameter dictionary
        parameter_dict = {
        (None, None, "k"):1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }
        PD = ParameterDatabase(parameter_dictionary=parameter_dict)
        return_dict = {(k.mechanism, k.part_id, k.name):PD.parameters[k].value for k in PD.parameters}
        self.assertEqual(return_dict, parameter_dict)

        #improper parameter dictionary
        k = ("M", "k")
        parameter_dict = {k:2.0}
        with self.assertRaisesRegex(ValueError, f"parameter_key must be"):
            PD = ParameterDatabase(parameter_dictionary = parameter_dict)

        #duplicate parameter dictionary
        key = (None, None, "k")
        parameter_dict = {"k":1, (None, None, "k"):1}
        with self.assertRaisesRegex(ValueError, f"Duplicate parameter detected"):
            PD = ParameterDatabase(parameter_dictionary = parameter_dict)

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

    def test_len(self):
        parameter_dict = {
        "k":1,
        ("M", "pid", "k"):2.0,
        (None, "pid", "k"):3.3,
        ("M", None, "k"): 4,
        }

        PD = ParameterDatabase(parameter_dictionary = parameter_dict)
        self.assertTrue(len(PD) == len(parameter_dict))

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
        with self.assertRaisesRegex(ValueError, f"parameter_key must be"):
            PD[("M", "k")]

        #test accessing something not in the PD
        with self.assertRaisesRegex(KeyError, f"ParameterKey"):
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
        with self.assertRaisesRegex(ValueError, f"Parameter Key does not match"):
            PD["test"] = PE


        #Invalid parameter key
        with self.assertRaisesRegex(ValueError, f"parameter_key must be"):
            PD[("M", "k")] = 10

        #Test overwriting
        PD["k"] = .1
        self.assertTrue(PD["k"].value == .1)
        PD[("M", "pid", "ktx")] = .333
        self.assertTrue(PD[("M", "pid", "ktx")].value == .333)

    def test_find_parameter(self):

        """Test the parameter defaulting heirarchy
        Parameter defaulting heirarchy:
        (mechanism_name, part_id, param_name) --> param_val. If that particular parameter key cannot be found, 
        the software will default to the following keys: 
        (mechanism_type, part_id, param_name) >> (part_id, param_name) >> 
        (mechanism_name, param_name) >> (mechanism_type, param_name) >>
        (param_name) and give a warning. """

        parameter_dict = {
            (None, None, "k"):1.0,
            ("M", None, "k"): 2.1,
            ("m", None, "k"): 2.2,
            (None, "pid", "k"):3,
            ("M", "pid", "k"):4.1,
            ("m", "pid", "k"):4.2,
        }

        M1 = Mechanism(name = "m", mechanism_type = "M")
        M2 = Mechanism(name = "m2", mechanism_type = "M")
        M3 = Mechanism(name = "m3", mechanism_type = "M2")

        PD = ParameterDatabase(parameter_dictionary = parameter_dict)

        self.assertEqual(PD.find_parameter(mechanism = M3, part_id="id", param_name = "k").value, 1.0)
        self.assertEqual(PD.find_parameter(mechanism = M2, part_id="id", param_name = "k").value, 2.1)
        self.assertEqual(PD.find_parameter(mechanism = M1, part_id="id", param_name = "k").value, 2.2)
        self.assertEqual(PD.find_parameter(mechanism = M3, part_id="pid", param_name = "k").value, 3.0)
        self.assertEqual(PD.find_parameter(mechanism = M2, part_id="pid", param_name = "k").value, 4.1)
        self.assertEqual(PD.find_parameter(mechanism = M1, part_id="pid", param_name = "k").value, 4.2)




