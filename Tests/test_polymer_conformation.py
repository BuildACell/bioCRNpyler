#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import *

from biocrnpyler import OrderedPolymerSpecies, PolymerConformation, ComplexSpecies, Species

def test_single_polymer_instantiation():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p[0], p[1]])


