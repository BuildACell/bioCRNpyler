#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies, Species, Complex, MonomerCollection

class TestMonomerCollection(TestCase):

    def test_initialization(self):
        x = OrderedMonomer()
        y = OrderedMonomer()
        c = MonomerCollection([x, y])

        #test setter
        assert len(c.monomers) == 2
        assert type(c.monomers) == tuple
        #test copying
        assert x not in c.monomers and y not in c.monomers
        assert x.parent is None and y.parent is None
        assert [m.parent == c for m in c.monomers]

    def test_monomers_setter_and_getter(self):
        x = OrderedMonomer()
        y = OrderedMonomer()
        c = MonomerCollection([x, y])

        xr = OrderedMonomer(direction = "reverse")

        c.monomers = [x, xr, y]
        assert len(c.monomers) == 3
        assert c.monomers[1].direction == "reverse"
        assert all([m.parent is None for m in [x, xr, y]])
        assert all([m.parent is c for m in c.monomers])

class TestOrderedMonomer(TestCase):

    def test_ordered_monomer_initialization(self):

        #Correct instantiation
        x = OrderedMonomer()
        x = OrderedMonomer(direction="forward")
        self.assertEqual(x.direction,"forward")
        self.assertEqual(x.position,None)
        self.assertEqual(x.parent,None)
        self.assertEqual(OrderedMonomer(direction="reverse"),x.set_dir("reverse"))

        #Bad parent
        with self.assertRaisesRegex(ValueError, "parent must be an MonomerCollection"):
            m = OrderedMonomer(parent = 1)
        
        p = OrderedPolymer(parts = [])

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

    def test_naming_convention(self):
        A = Species("A", material_type = "a")
        B = Species('B', attributes = "b")
        C = Complex([Species("S"), Species("S")])
        p = OrderedPolymerSpecies([A, B, C], attributes = ["a"])
        print(str(p))
        print(f"{p.material_type}_{str(A)}_{str(B)}_{str(C)}_{p.attributes[0]}_")
        self.assertTrue(str(p) == f"{p.material_type}_{str(A)}_{str(B)}_{str(C)}_{p.attributes[0]}_")

    def test_ordered_polymer_species_initialization(self):
        a = Species("A")
        x = OrderedPolymerSpecies([Species("A"),[Species("B"),"forward"],\
                            Species("C").set_dir("reverse")],attributes=["ooga"])
        #repr
        revtest = OrderedPolymerSpecies([Species("A",direction="forward"),Species("A")])
        r1 = str(revtest)
        revtest.reverse()
        r2 = str(revtest)
        self.assertEqual(type(r1),str)
        self.assertNotEqual(r1,r2)
        #pretty_print
        self.assertEqual(type(x.pretty_print()),str)
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
        self.assertIn("circular",x.attributes)
        x.circular = False
        self.assertEqual(x.circular, False)
        self.assertNotIn("circular",x.attributes)
        
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
        print(reversd)
        print(truth)
        self.assertEqual(reversd,truth)

    def test_ordered_polymer_species_contains(self):
        a = Species("A")
        b = Species("B")
        bf = Species("B").set_dir("forward")
        c = Complex([a, a])
        p = OrderedPolymerSpecies([bf, b, c])

        #In these cases, parent doesn't matter
        self.assertTrue(a in p[2])
        self.assertTrue(a in c)
        self.assertTrue(b in p)
        self.assertTrue(c in p)

        p2 = OrderedPolymerSpecies([bf, a, c])
        #In these cases parents matter
        self.assertFalse(p[0] in p2)
        self.assertFalse(p2[0] in p)

        #In this case, direciton matters
        self.assertFalse(b in p2)

    def test_ordered_polymer_monomer_eq(self):
        a = Species("A")
        b = Species("B")
        bf = Species("B").set_dir("forward")
        br = Species("B").set_dir("reverse")
        c = Complex([a, a])
        p = OrderedPolymerSpecies([bf, b, c, br])

        self.assertFalse(b == p[1]) #Equality depends on the position and parent

        #monomer_equality does not depend on parent or position (direction still matters if they are different)
        self.assertTrue(b.monomer_eq(p[0])) 
        self.assertTrue(b.monomer_eq(p[1]))
        self.assertTrue(bf.monomer_eq(p[0])) 
        self.assertTrue(bf.monomer_eq(p[1]))

        #Also works for things in the same Parent
        self.assertTrue(p[0].monomer_eq(p[1]))

        #Also works for Complexes
        self.assertTrue(c.monomer_eq(p[2]))
        
        #In this case, direction will matter because it is not None in both Monomers
        self.assertFalse(p[0].monomer_eq(p[3]))
        