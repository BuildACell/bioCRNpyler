#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies, Species

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

        self.assertFalse(y[0]==x)#Parts should be copied, not inserted!

        #Test faithful copying of a part
        self.assertEqual(y[0].direction, x.direction)
        self.assertEqual(y[0].parent, y)
        self.assertEqual(y[0].position, 0)

        self.assertEqual(y[1].position,1)
        self.assertEqual(y[1].direction,"reverse")
        self.assertEqual(y[1].parent,y)
        self.assertEqual(y[2].position,2)
        self.assertEqual(y[2].direction,None)
        self.assertEqual(y[2].parent,y)

    def test_ordered_polymer_manipulations(self):
        z = OrderedMonomer(direction="forward")
        x = OrderedMonomer(direction="reverse")
        y = OrderedPolymer([x,[OrderedMonomer(),"reverse"],OrderedMonomer()])
        #replace
        y.replace(1,z)
        truthvalue = OrderedPolymer([OrderedMonomer(direction="reverse"),OrderedMonomer(direction="forward"),OrderedMonomer(direction=None)])
        self.assertEqual(y,truthvalue)
        #insert
        y = OrderedPolymer([x,OrderedMonomer()])
        ins = OrderedMonomer(direction="forward")
        y.insert(1,ins)
        self.assertEqual(y,truthvalue)
        #append
        y = OrderedPolymer([x,OrderedMonomer(direction="forward")])
        y.append(OrderedMonomer())
        self.assertEqual(y,truthvalue)
        #reverse
        y = OrderedPolymer([OrderedMonomer(),OrderedMonomer(direction="reverse"),\
                                    OrderedMonomer(direction="forward")])
        y.reverse()
        self.assertEqual(y,truthvalue)
        #delpart
        y = OrderedPolymer([OrderedMonomer(direction="reverse"),\
                            OrderedMonomer(direction="forward"),\
                            OrderedMonomer(direction="reverse"),\
                            OrderedMonomer()])
        z = y[2]
        self.assertEqual(z.parent,y)
        y.delpart(2)
        self.assertEqual(z.parent,None)
        self.assertEqual(y,truthvalue)
class TestOrderedPolymerSpecies(TestCase):
    def test_ordered_polymer_species_initialization(self):
        a = Species("A")
        x = OrderedPolymerSpecies([Species("A"),[Species("B"),"forward"],\
                            Species("C").set_dir("reverse")],attributes=["ooga"])
        #repr
        self.assertEqual(str(x),"ordered_polymer_A_B_f_C_r_ooga")
        #pretty_print
        self.assertEqual(x.pretty_print(),"ordered_polymer[A:B-forward:C-reverse(ooga)]")
        #make sure we're copying
        self.assertEqual(a.parent,None)
        #make sure we're setting the properties of the components
        self.assertEqual(x[0].parent,x)
        self.assertEqual(x[1].parent,x)
        self.assertEqual(x[2].parent,x)
        self.assertEqual(x[1].direction,"forward")
        self.assertEqual(x[2].direction,"reverse")

        self.assertEqual(x[0].position,0)
        self.assertEqual(x[1].position,1)
        self.assertEqual(x[2].position,2)

        #circularity
        x.circular = True
        self.assertEqual(x.circular, True)
        
    def test_ordered_polymer_species_manipulations(self):
        a = Species("A")
        b = Species("B").set_dir("forward")
        c = Species("C").set_dir("reverse")
        truth = OrderedPolymerSpecies([Species("A"),[Species("B"),"forward"],\
                            Species("C").set_dir("reverse")],attributes=["ooga"])
        reversd = OrderedPolymerSpecies([Species("C").set_dir("forward"),[Species("B"),"reverse"],\
                            Species("A")],attributes=["ooga"])
        unreplaced = OrderedPolymerSpecies([Species("z"),[Species("B"),"forward"],\
                            Species("C").set_dir("reverse")],attributes=["ooga"])
        uninserted = OrderedPolymerSpecies([Species("A"),\
                            Species("C").set_dir("reverse")],attributes=["ooga"])
        unappended = OrderedPolymerSpecies([Species("A"),[Species("B"),"forward"],\
                            ],attributes=["ooga"])
        #replace
        unreplaced.replace(0,a)
        self.assertEqual(unreplaced,truth)
        #insert
        uninserted.insert(1,b)
        self.assertEqual(uninserted,truth)
        #append
        unappended.append(c)
        self.assertEqual(unappended,truth)
        #reverse
        reversd.reverse()
        self.assertEqual(reversd,truth)