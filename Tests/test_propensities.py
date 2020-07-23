#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import MassAction
from biocrnpyler import ParameterEntry
import pytest


def test_massaction_forward_rate():

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=0)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=-1)


def test_massaction_reserve_rate():
    with pytest.raises(TypeError, match=r"missing 1 required positional argument: 'k_forward'"):
        MassAction(k_reverse=0.1)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=1, k_reverse=0)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=1, k_reverse=-1)


def test_massaction_is_reverable():
    mak = MassAction(k_forward=1, k_reverse=0.1)
    assert mak.is_reversible

def test_massaction_from_parameters():
    mak = MassAction(k_forward = ParameterEntry("k", .1))
    assert(mak.k_forward == .1)

    mak = MassAction(k_forward = ParameterEntry("k", .1), k_reverse = ParameterEntry("k", .01))
    assert(mak.k_reverse == .01)




