#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies

class TestOrderedMonomer(TestCase):

    def test_ordered_monomer_initialization(self):

        #Correct instantiation
        x = OrderedMonomer()
        x = OrderedMonomer(direction="forward")
        self.assertEqual(x.direction,"forward")
        self.assertEqual(x.position,None)
        self.assertEqual(x.parent,None)
        self.assertEqual(OrderedMonomer(direction="reverse"),x.set_dir("reverse"))

        #Position with no parent
        with self.assertRaisesRegexp(ValueError, f"OrderedMonomer"):
            m = OrderedMonomer(position = 1)

        #Bad parent
        with self.assertRaisesRegexp(ValueError, f"parent must be an OrderedPolymer"):
            m = OrderedMonomer(parent = 1)
        
        p = OrderedPolymer(parts = [])
        #Parent with no position
        with self.assertRaisesRegexp(ValueError, f"OrderedMonomer"):
            m = OrderedMonomer(parent = p)

        #Correct instantiation with parent and position
        m = OrderedMonomer(parent = p, position = 0)

    def test_properties(self):
        x = OrderedMonomer()

        #Test None initialization of properties
        self.assertTrue(x.parent is None)
        self.assertTrue(x.direction is None)
        self.assertTrue(x.position is None)

        p = OrderedPolymer(parts = [])
        x = OrderedMonomer(parent = p, position = 0, direction = "r")
        self.assertTrue(x.parent == p)
        self.assertTrue(x.position == 0)
        self.assertTrue(x.direction == "r")


class TestOrderedPolymer(TestCase):

    def test_ordered_polymer_intialization(self):

        x = OrderedMonomer(direction="forward")
        y = OrderedPolymer([x,[OrderedMonomer(),"reverse"],OrderedMonomer()])
        self.assertEqual(y[0],x)
        self.assertEqual(y[1].position,1)
        self.assertEqual(y[1].direction,"reverse")
        self.assertEqual(y[1].parent,y)
        self.assertEqual(y[2].position,2)
        self.assertEqual(y[2].direction,None)
        self.assertEqual(y[2].parent,y)

    def test_ordered_polymer_manipulations(self):
        z = OrderedMonomer(direction="forward")
        x = OrderedMonomer(direction="forward")
        y = OrderedPolymer([x,[OrderedMonomer(),"reverse"],OrderedMonomer()])
        y.replace(1,z)
        truthvalue = OrderedPolymer([OrderedMonomer(direction="forward"),OrderedMonomer(direction="forward"),OrderedMonomer(direction=None)])
        self.assertEqual(y,truthvalue)