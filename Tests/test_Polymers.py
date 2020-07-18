#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies

class TestOrderedMonomer(TestCase):
    def test_ordered_monomer_initialization(self):
        x = OrderedMonomer(direction="forward")
        self.assertEqual(x.direction,"forward")