#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species
from biocrnpyler.utils import process_initial_concentration_dict


def test_process_initial_concentration_dict():
    A = Species('A')
    B = Species('B')
    C = Species('C')
    initial_concentration = {A: 10, B: 25.34, str(C): 77.24}
    processed = process_initial_concentration_dict(initial_concentration)
    for (key_old, key_new) in zip(initial_concentration.keys(), processed.keys()):
        assert str(key_old) == key_new
