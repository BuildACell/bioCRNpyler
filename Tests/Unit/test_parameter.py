#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import sys
from unittest import TestCase
from unittest.mock import mock_open, patch
from warnings import warn

from biocrnpyler import (
    Mechanism,
    ModelParameter,
    Parameter,
    ParameterDatabase,
    ParameterEntry,
)


class TestParameter(TestCase):
    def test_parameter(self):
        # test parameter name
        with self.assertRaisesRegex(
            ValueError, 'parameter_name must be a string'
        ):
            Parameter(parameter_name=None, parameter_value=1.0)
        # test parameter value
        with self.assertRaisesRegex(
            ValueError, 'parameter_value must be a float or int'
        ):
            Parameter(parameter_name='None', parameter_value=None)

        # test invalid value string
        with self.assertRaisesRegex(
            ValueError, 'No valid parameter value! Accepted format'
        ):
            Parameter(parameter_name='None', parameter_value='2ba')

        # test invalid unit value
        with self.assertRaisesRegex(ValueError, 'All units must be strings'):
            Parameter(parameter_name='None', parameter_value=1.0, unit=0.0)

        # test string parameter values
        self.assertTrue(
            Parameter(parameter_name='None', parameter_value='1.0').value
            == 1.0
        )

        self.assertTrue(
            Parameter(parameter_name='None', parameter_value='1/2').value
            == 0.5
        )

        self.assertTrue(
            Parameter(parameter_name='None', parameter_value='1e2').value
            == 100
        )

        self.assertTrue(
            Parameter(parameter_name='None', parameter_value='1e-2').value
            == 0.01
        )
        self.assertTrue(
            Parameter(parameter_name='None', parameter_value='2e-2/2').value
            == 0.01
        )

        # Test unit setting
        self.assertTrue(
            Parameter(
                parameter_name='None', parameter_value='1.0', unit=None
            ).unit
            == ''
        )
        self.assertTrue(
            Parameter(
                parameter_name='None', parameter_value='1.0', unit='M'
            ).unit
            == 'M'
        )

        # testing invalid parameter name
        with self.assertRaisesRegex(
            ValueError,
            'parameter_name should be at least one character and cannot start with a number!',
        ):
            Parameter(parameter_name='2', parameter_value=2)

    def test_parameter_entry(self):
        # Valid ParameterEntry Construction
        # With unit
        P0 = ParameterEntry(
            parameter_name='None',
            parameter_value=1.0,
            parameter_key={'part_id': 'id'},
            parameter_info={'comment': 'comment', 'unit': 'M'},
        )
        P1 = ParameterEntry(
            parameter_name='None',
            parameter_value=1.0,
            parameter_key={'part_id': 'id'},
            parameter_info={'comment': 'comment'},
            unit='M',
        )

        # Without unit
        P2 = ParameterEntry(
            parameter_name='None',
            parameter_value=1.0,
            parameter_key={'part_id': 'id'},
            parameter_info={'comment': 'comment'},
        )

        # Assert unit is passed through to Parameter.unit from the parameter_info_dictionary
        self.assertTrue(P0.unit == 'M')
        self.assertTrue(P1.unit == 'M')
        self.assertTrue(P2.unit == '')

        # Test duplication of parameter information error
        with self.assertRaisesRegex(
            ValueError, 'Recieved multiple parameter units'
        ):
            P0 = ParameterEntry(
                parameter_name='None',
                parameter_value=1.0,
                parameter_key={'part_id': 'id'},
                parameter_info={'comment': 'comment', 'unit': 'M'},
                unit='uM',
            )

        # Invalid keys
        param_keys = Parameter(parameter_name='None', parameter_value=1.0)
        with self.assertRaisesRegex(ValueError, 'parameter_key must be'):
            ParameterEntry(
                parameter_name='None',
                parameter_value=1.0,
                parameter_key=param_keys,
            )

        # Invalid info
        param_info = 'blah blah'
        with self.assertRaisesRegex(
            ValueError, 'parameter_info must be None or a dictionary'
        ):
            ParameterEntry(
                parameter_name='None',
                parameter_value=1.0,
                parameter_info=param_info,
            )

    def test_model_parameter(self):
        # valid ModelParameter Construction
        ModelParameter(
            parameter_name='None',
            parameter_value=1.0,
            search_key=('that', 'this', 'k'),
            found_key=('this', None, 'k'),
        )

        # Invalid keys
        k = Parameter(parameter_name='None', parameter_value='1.0')
        with self.assertRaisesRegex(ValueError, 'parameter_key must be None'):
            ModelParameter(
                parameter_name='None',
                parameter_value=1.0,
                search_key=k,
                found_key=('this', None, 'k'),
            )

        with self.assertRaisesRegex(ValueError, 'parameter_key must be None'):
            ModelParameter(
                parameter_name='None',
                parameter_value=1.0,
                search_key=('that', 'this', 'k'),
                found_key=k,
            )

    def test_get_field_names(self):
        PD = ParameterDatabase
        test_accepted_field_names = {
            'mechanism': ['mechanism', 'mechanism_id']
        }
        valid_field_names = ['part_id']

        # test None as field_names
        with self.assertRaisesRegex(
            ValueError, 'field_names must be a list of strings'
        ):
            PD._get_field_names(
                field_names=None,
                accepted_field_names=test_accepted_field_names,
            )
        # test invalid field_names type
        with self.assertRaisesRegex(
            ValueError, 'field_names must be a list of strings'
        ):
            PD._get_field_names(
                field_names={}, accepted_field_names=test_accepted_field_names
            )
        # test empty field_names list
        with self.assertRaisesRegex(
            ValueError, 'field_names cannot be empty list!'
        ):
            PD._get_field_names(
                field_names=[], accepted_field_names=test_accepted_field_names
            )

        # test None as accepted_field_names
        with self.assertRaisesRegex(
            ValueError, 'accepted_field_names must be a dictionary'
        ):
            PD._get_field_names(
                field_names=valid_field_names, accepted_field_names=None
            )
        # test invalid accepted_field_names type
        with self.assertRaisesRegex(
            ValueError, 'accepted_field_names must be a dictionary'
        ):
            PD._get_field_names(
                field_names=valid_field_names, accepted_field_names=[]
            )
        # test empty field_names list
        with self.assertRaisesRegex(
            ValueError, 'accepted_field_names cannot be empty dictionary'
        ):
            PD._get_field_names(
                field_names=valid_field_names, accepted_field_names={}
            )

        accepted_field_names = {
            'mechanism': ['mechanism', 'mechanism_id'],
            'param_name': [
                'parameter_name',
                'parameter',
                'param',
                'param_name',
            ],
            'part_id': ['part_id', 'part'],
            'param_val': ['val', 'value', 'param_val', 'parameter_value'],
        }

        ret_dict = PD._get_field_names(
            field_names=[''], accepted_field_names=accepted_field_names
        )
        self.assertEqual(accepted_field_names.keys(), ret_dict.keys())

        with self.assertWarns(Warning):
            PD._get_field_names(valid_field_names, accepted_field_names)

        accepted_field_names = {
            'dummy': ['dumb', 'dumber'],
        }

        with self.assertWarns(Warning):
            PD._get_field_names(valid_field_names, accepted_field_names)

    def test_load_parameters_from_file(self):
        # Bad parameter file keyword
        with self.assertRaisesRegex(
            ValueError,
            'parameter_file must be a string representing a file name and path.',
        ):
            ParameterDatabase(parameter_file={})

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            # !!! DO NOT reformat this string below !!!!
            example_csv = """mechanism_id	part_id	param_name	param_val	comments\ntranscription_mm	ptet_tetR	kb	10.	extra columns are okay!\ntranscription_mm	ptet_tetR	ku	.1	These are the parameters for transcription"""
            # !!! DO NOT reformat this string above !!!!

            from pathlib import Path

            with (
                patch.object(Path, 'exists', return_value=True),
                patch.object(Path, 'is_file', return_value=True),
                patch(
                    'builtins.open',
                    mock_open(read_data=example_csv),
                    create=True,
                ),
            ):
                PD = ParameterDatabase(parameter_file='test_file.tsv')

                right_dict = {
                    ('transcription_mm', 'ptet_tetR', 'kb'): 10.0,
                    ('transcription_mm', 'ptet_tetR', 'ku'): 0.1,
                }
                returned_dict = {
                    (k.mechanism, k.part_id, k.name): PD.parameters[k].value
                    for k in PD.parameters
                }
                # raise ValueError(str(returned_dict))
                self.assertEqual(right_dict, returned_dict)
        else:
            warn('version below 3.6 was detected! This test was skipped')

        # TODO track down why this test fails in python 3.6!
        if sys.version_info[1] >= 7:
            example_csv = """mechanism_id"""

            with (
                patch.object(Path, 'exists', return_value=True),
                patch.object(Path, 'is_file', return_value=True),
                patch(
                    'builtins.open',
                    mock_open(read_data=example_csv),
                    create=True,
                ),
            ):
                with self.assertWarnsRegex(
                    Warning,
                    'No param_name column was found, could not load parameter!',
                ):
                    ParameterDatabase(parameter_file='test_file.tsv')

        else:
            warn('version below 3.6 was detected! This test was skipped')

    def test_load_parameters_from_file_preserves_units_column(self):
        if sys.version_info[1] >= 7:
            example_tsv = (
                "mechanism_id\tparam_name\tparam_val\tunits\n"
                "transcription\tktx\t0.05\t/sec\n"
                "binding\tku\t10\tuM/sec"
            )

            from pathlib import Path

            with (
                patch.object(Path, 'exists', return_value=True),
                patch.object(Path, 'is_file', return_value=True),
                patch(
                    'builtins.open',
                    mock_open(read_data=example_tsv),
                    create=True,
                ),
            ):
                pd = ParameterDatabase(parameter_file='test_file.tsv')

            units_by_key = {
                (k.mechanism, k.part_id, k.name): pd.parameters[k].unit
                for k in pd.parameters
            }
            self.assertEqual(
                units_by_key[('transcription', None, 'ktx')], 'per_sec'
            )
            self.assertEqual(
                units_by_key[('binding', None, 'ku')], 'uM_per_sec'
            )
        else:
            warn('version below 3.6 was detected! This test was skipped')

    def test_load_parameters_from_dictionary(self):
        # bad parameter_dictionary keyword
        with self.assertRaisesRegex(
            ValueError, 'parameter_dictionary must be None or a dictionary!'
        ):
            PD = ParameterDatabase(parameter_dictionary='test_file')

        # proper parameter dictionary
        parameter_dict = {
            (None, None, 'k'): 1,
            ('M', 'pid', 'k'): 2.0,
            (None, 'pid', 'k'): 3.3,
            ('M', None, 'k'): 4,
        }
        PD = ParameterDatabase(parameter_dictionary=parameter_dict)
        return_dict = {
            (k.mechanism, k.part_id, k.name): PD.parameters[k].value
            for k in PD.parameters
        }
        self.assertEqual(return_dict, parameter_dict)

        # improper parameter dictionary
        k = ('M', 'k')
        parameter_dict = {k: 2.0}
        with self.assertRaisesRegex(ValueError, 'parameter_key must be'):
            PD = ParameterDatabase(parameter_dictionary=parameter_dict)

        # duplicate parameter dictionary
        key = (None, None, 'k')  # noqa: F841
        parameter_dict = {'k': 1, (None, None, 'k'): 1}
        with self.assertRaisesRegex(
            ValueError, 'Duplicate parameter detected'
        ):
            PD = ParameterDatabase(parameter_dictionary=parameter_dict)

    def test_iterator(self):
        parameter_dict = {
            'k': 1,
            ('M', 'pid', 'k'): 2.0,
            (None, 'pid', 'k'): 3.3,
            ('M', None, 'k'): 4,
        }

        PD = ParameterDatabase(parameter_dictionary=parameter_dict)

        # All things in iterator are entries
        count = 0
        for entry in PD:
            count += 1
            self.assertTrue(isinstance(entry, ParameterEntry))

        # The correct number of entries
        self.assertEqual(count, len(parameter_dict))

    def test_len(self):
        parameter_dict = {
            'k': 1,
            ('M', 'pid', 'k'): 2.0,
            (None, 'pid', 'k'): 3.3,
            ('M', None, 'k'): 4,
        }

        PD = ParameterDatabase(parameter_dictionary=parameter_dict)
        self.assertTrue(len(PD) == len(parameter_dict))

    def test_contains(self):
        parameter_dict = {
            'k': 1,
            ('M', 'pid', 'k'): 2.0,
            (None, 'pid', 'k'): 3.3,
            ('M', None, 'k'): 4,
        }

        PD = ParameterDatabase(parameter_dictionary=parameter_dict)

        self.assertTrue('k' in PD)
        self.assertTrue(('M', 'pid', 'k') in PD)
        self.assertTrue(PD['k'] in PD)
        self.assertFalse('ktx' in PD)
        self.assertFalse(('M', 'k') in PD)

    def test_indexing(self):
        parameter_dict = {
            'k': 1,
            ('M', 'pid', 'k'): 2.0,
            (None, 'pid', 'k'): 3.3,
            ('M', None, 'k'): 4,
        }

        PD = ParameterDatabase(parameter_dictionary=parameter_dict)

        # Test correct accessing
        self.assertTrue(PD['k'].value == 1)
        self.assertTrue(PD[('M', 'pid', 'k')].value == 2.0)
        self.assertTrue(PD[(None, 'pid', 'k')].value == 3.3)
        self.assertTrue(PD[('M', None, 'k')].value == 4)

        # test incorrect accessing
        with self.assertRaisesRegex(ValueError, 'parameter_key must be'):
            PD[('M', 'k')]

        # test accessing something not in the PD
        with self.assertRaisesRegex(KeyError, 'ParameterKey'):
            PD['kb']

        # test inserting values
        PD['kb'] = 100
        self.assertTrue(PD['kb'].value == 100)
        PD[(None, None, 'ku')] = 200
        self.assertTrue(PD['ku'].value == 200)
        PD[('M', 'pid', 'ktx')] = 300
        self.assertTrue(PD[('M', 'pid', 'ktx')].value == 300)

        # Test inserting ParameterEntry
        PE = ParameterEntry('test', 1.0)
        PD['test'] = PE
        self.assertTrue(PE in PD)

        # Test Correct Overwriting
        PE = ParameterEntry('test', 2.0)
        PD['test'] = PE
        self.assertTrue(PD['test'].value == PE.value)

        # Test incorrect Overwriting
        PE = ParameterEntry('t', 1.0)
        with self.assertRaisesRegex(
            ValueError, 'Parameter Key does not match'
        ):
            PD['test'] = PE

        # Invalid parameter key
        with self.assertRaisesRegex(ValueError, 'parameter_key must be'):
            PD[('M', 'k')] = 10

        # Test overwriting
        PD['k'] = 0.1
        self.assertTrue(PD['k'].value == 0.1)
        PD[('M', 'pid', 'ktx')] = 0.333
        self.assertTrue(PD[('M', 'pid', 'ktx')].value == 0.333)

    def test_find_parameter(self):
        """Test the parameter defaulting heirarchy
        Parameter defaulting heirarchy:
        (mechanism_name, part_id, param_name) --> param_val. If that particular parameter key cannot be found,
        the software will default to the following keys:
        (mechanism_type, part_id, param_name) >> (part_id, param_name) >>
        (mechanism_name, param_name) >> (mechanism_type, param_name) >>
        (param_name) and give a warning."""

        parameter_dict = {
            (None, None, 'k'): 1.0,
            ('M', None, 'k'): 2.1,
            ('m', None, 'k'): 2.2,
            (None, 'pid', 'k'): 3,
            ('M', 'pid', 'k'): 4.1,
            ('m', 'pid', 'k'): 4.2,
        }

        M1 = Mechanism(name='m', mechanism_type='M')
        M2 = Mechanism(name='m2', mechanism_type='M')
        M3 = Mechanism(name='m3', mechanism_type='M2')

        PD = ParameterDatabase(parameter_dictionary=parameter_dict)

        self.assertEqual(
            PD.find_parameter(
                mechanism=M3, part_id='id', param_name='k'
            ).value,
            1.0,
        )
        self.assertEqual(
            PD.find_parameter(
                mechanism=M2, part_id='id', param_name='k'
            ).value,
            2.1,
        )
        self.assertEqual(
            PD.find_parameter(
                mechanism=M1, part_id='id', param_name='k'
            ).value,
            2.2,
        )
        self.assertEqual(
            PD.find_parameter(
                mechanism=M3, part_id='pid', param_name='k'
            ).value,
            3.0,
        )
        self.assertEqual(
            PD.find_parameter(
                mechanism=M2, part_id='pid', param_name='k'
            ).value,
            4.1,
        )
        self.assertEqual(
            PD.find_parameter(
                mechanism=M1, part_id='pid', param_name='k'
            ).value,
            4.2,
        )


def test_findpath(tmp_path):
    import os
    import platform

    import biocrnpyler as bcp

    # Make sure files that don't exist return None
    assert bcp.find_file_in_bcp_path('does_not_exist') is None

    # Make sure that files in package can be found
    assert os.path.exists(
        bcp.find_file_in_bcp_path('mechanisms/txtl_parameters.tsv')
    )

    prev_cwd = os.getcwd()
    prev_bcp_path = os.environ.get('BCP_PATH')
    child_dir = tmp_path / 'child'
    child_dir.mkdir()
    os.chdir(child_dir)
    try:
        # Make sure we can find files in current directory
        assert bcp.find_file_in_bcp_path('__testfile__.tsv') is None
        (child_dir / '__testfile__.tsv').write_text('')
        assert os.path.exists(bcp.find_file_in_bcp_path('__testfile__.tsv'))
        os.remove(child_dir / '__testfile__.tsv')

        # Make sure we can find files in the path
        (tmp_path / '__testfile__.tsv').write_text('')
        assert bcp.find_file_in_bcp_path('__testfile__.tsv') is None
        if platform.system() == 'Windows':
            os.environ['BCP_PATH'] = f'/tmp;.;{tmp_path}'
        else:
            os.environ['BCP_PATH'] = f'/tmp:.:{tmp_path}'
        assert os.path.exists(bcp.find_file_in_bcp_path('__testfile__.tsv'))
    finally:
        os.chdir(prev_cwd)
        if prev_bcp_path is None:
            os.environ.pop('BCP_PATH', None)
        else:
            os.environ['BCP_PATH'] = prev_bcp_path


def test_parameter_ordering():
    import biocrnpyler as bcp

    # Create the components for testing
    TetR = bcp.Protein('TetR')
    aTc = bcp.Species('aTc')
    TetR_inactive=bcp.ChemicalComplex([TetR.species, aTc])
    ptet = bcp.RegulatedPromoter('ptet', TetR)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP')

    # Default case: make sure we can compile with no additional paramters
    default_mixture = bcp.BasicPURE(
        'default', components=[dna_GFP, TetR_inactive])
    default_crn = default_mixture.compile_crn()

    # Find the binding reactions of TetR to DNA and aTc
    aTc_TetR_index = dna_TetR_index = -1
    for i, reaction in enumerate(default_crn.reactions):
        if not isinstance(reaction.propensity_type, bcp.MassAction):
            continue

        inputs = [input.species for input in reaction.inputs]

        # Find the binding of aTc to TetR
        if aTc in inputs and TetR.species in inputs:
            aTc_TetR_index = i

        # Find the binding of TetR to DNA
        if dna_GFP.species in inputs and TetR.species in inputs:
            dna_TetR_index = i

    assert aTc_TetR_index != -1 and dna_TetR_index != -1
    assert default_crn.reactions[
        dna_TetR_index].propensity_type.k_forward != 1
    assert default_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward != 1

    # Rebuild using given parameters
    baseline_parameters = {
        ('binding', None, 'kb'): 1,
        ('binding', None, 'ku'): 0.1,
        ('binding', None, 'cooperativity'): 2,
    }
    TetR_inactive=bcp.ChemicalComplex(
        [TetR.species, aTc], parameters=baseline_parameters)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP',
        parameters=baseline_parameters)
    baseline_mixture = bcp.BasicPURE(
        'baseline', components=[dna_GFP, TetR_inactive])
    baseline_crn = baseline_mixture.compile_crn()
    assert baseline_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 1
    assert baseline_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 1

    # Override using mechanism subtypes
    mechtype_parameters = baseline_parameters | {
        ('binding', 'dna_protein', 'kb'): 2,
        ('binding', 'chemical_complex', 'kb'): 3,
    }
    TetR_inactive=bcp.ChemicalComplex(
        [TetR.species, aTc], parameters=mechtype_parameters)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP',
        parameters=mechtype_parameters)
    mechtype_mixture = bcp.BasicPURE(
        'mechtype', components=[dna_GFP, TetR_inactive])
    mechtype_crn = mechtype_mixture.compile_crn()
    assert mechtype_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 2
    assert mechtype_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 3

    # Override using part IDs
    partid_parameters = mechtype_parameters | {
        ('binding', 'ptet_TetR', 'kb'): 4,
        ('binding', 'aTc_protein_TetR', 'kb'): 5,
    }
    TetR_inactive=bcp.ChemicalComplex(
        [TetR.species, aTc], parameters=partid_parameters)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP',
        parameters=partid_parameters)
    partid_mixture = bcp.BasicPURE(
        'partid', components=[dna_GFP, TetR_inactive])
    partid_crn = partid_mixture.compile_crn()
    assert partid_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 4
    assert partid_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 5

    # Override using mechanism name
    mechname_parameters = baseline_parameters | {
        ('one_step_cooperative_binding', None, 'kb'): 6,
        ('binding', 'chemical_complex', 'kb'): 7,
    }
    TetR_inactive=bcp.ChemicalComplex(
        [TetR.species, aTc], parameters=mechname_parameters)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP',
        parameters=mechname_parameters)
    mechname_mixture = bcp.BasicPURE(
        'mechname', components=[dna_GFP, TetR_inactive])
    mechname_crn = mechname_mixture.compile_crn()
    assert mechname_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 6
    assert mechname_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 7

    # Override in mixture instead of components
    mixture_parameters = baseline_parameters | {
        ('one_step_cooperative_binding', None, 'kb'): 8,
        ('binding', 'chemical_complex', 'kb'): 9,
    }
    TetR_inactive=bcp.ChemicalComplex([TetR.species, aTc])
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP')
    mixture_mixture = bcp.BasicPURE(
        'mixture', components=[dna_GFP, TetR_inactive],
        parameters=mixture_parameters)
    mixture_crn = mixture_mixture.compile_crn()
    assert mixture_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 8
    assert mixture_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 9

    # Override in components and not mixture
    component_parameters = baseline_parameters | {
        ('one_step_cooperative_binding', None, 'kb'): 10,
        ('binding', 'chemical_complex', 'kb'): 11,
    }
    TetR_inactive=bcp.ChemicalComplex(
        [TetR.species, aTc], parameters=component_parameters)
    dna_GFP = bcp.DNAassembly(
        'GFP', promoter=ptet, rbs='RBS', protein='GFP',
        parameters=component_parameters)
    mixture_mixture = bcp.BasicPURE(
        'mixture', components=[dna_GFP, TetR_inactive],
        parameters=mixture_parameters)
    mixture_crn = mixture_mixture.compile_crn()
    assert mixture_crn.reactions[
        dna_TetR_index].propensity_type.k_forward == 10
    assert mixture_crn.reactions[
        aTc_TetR_index].propensity_type.k_forward == 11
