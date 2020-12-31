
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler.sbmlutil import create_sbml_model
import biocrnpyler.units as units
import pytest
import libsbml
from biocrnpyler.units import *

def check(value, message):
        """If 'value' is None, prints an error message constructed using
        'message' and then exits with status code 1 (for libsbml). If 'value' is an integer,
        it assumes it is a libSBML return status code.  If the code value is
        LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
        prints an error message constructed using 'message' along with text from
        libSBML explaining the meaning of the code, and exits with status code 1.
        """
        if value == None:
            raise SystemExit(
                'LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value == libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                    + 'LibSBML returned error code ' + str(value) + ': "' \
                    + libsbml.OperationReturnValue_toString(value).strip() + '"'
                raise SystemExit(err_msg)
        else:
            return

class TestUnits(TestCase):
    def test_units_initialization(self):
        document, sbml_model = create_sbml_model()
        supported_units = biocrnpyler_supported_units().keys()
        for unit in supported_units:
            unit_definition = create_new_unit_definition(sbml_model, unit)
            check(unit_definition, "create new unit definition")
        with self.assertRaisesRegex(ValueError, 'The arguments are not of expected type. unit_id must be a string.'):
            create_new_unit_definition(sbml_model, 24)