# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from unittest import TestCase

import libsbml

from biocrnpyler.utils.sbmlutil import create_sbml_model
from biocrnpyler.utils.units import *


def check(value, message):
    """Utility function for printing SBML error messages.

    If `value` is None, prints an error message constructed using 'message'
    and then exits with status code 1 (for libsbml). If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.

    """
    if value == None:
        raise SystemExit(
            'LibSBML returned a null value trying to ' + message + '.'
        )
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = (
                'Error encountered trying to '
                + message
                + '.'
                + 'LibSBML returned error code '
                + str(value)
                + ': "'
                + libsbml.OperationReturnValue_toString(value).strip()
                + '"'
            )
            raise SystemExit(err_msg)
    else:
        return


class TestUnits(TestCase):
    def test_units_initialization(self):
        document, sbml_model = create_sbml_model()
        supported_units = biocrnpyler_supported_units().keys()
        for unit in supported_units:
            unit_definition = create_new_unit_definition(sbml_model, unit)
            check(unit_definition, 'create new unit definition')
        with self.assertRaisesRegex(
            ValueError,
            'Arguments are not of expected type. unit_id must be a string.',
        ):
            create_new_unit_definition(sbml_model, 24)


def test_supported_units():
    supported_units = biocrnpyler_supported_units()
    for key in [
        'nM',
        'uM',
        'mM',
        'M',
        'nL',
        'uL',
        'mL',
        'L',
        'sec',
        'min',
        'hrs',
    ]:
        assert key in supported_units

    for base, aliases in [
        ('second', ['s', 'sec', 'secs', 'seconds']),
        ('minute', ['min', 'mins', 'minutes']),
        ('hour', ['hr', 'hrs', 'hours']),
    ]:
        for alias in aliases:
            assert supported_units[alias] == supported_units[base]


def test_unit_conversion_variables():
    from biocrnpyler.utils.units import (
        nM,
        uM,
        mM,
        M,
        nL,
        uL,
        mL,
        L,
        sec,
        min,
        hrs,
    )

    assert nM == uM * 1e-3
    assert mM == uM * 1e3
    assert M == uM * 1e6

    assert nL == uL * 1e-3
    assert mL == uL * 1e3
    assert L == uL * 1e6

    assert min == sec * 60
    assert hrs == min * 60
