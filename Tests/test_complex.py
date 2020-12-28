#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies, Species, ComplexSpecies, Complex, OrderedComplexSpecies

class TestComplex(TestCase):
    def test_complex(self):
        a = Species('A')
        b = Species('B')
        c = Complex([a,b])
        d = OrderedPolymerSpecies([a,b,a])
        #normal complex that makes complex species
        self.assertEqual(set(c.species),set([a,b]))
        self.assertTrue(type(c)==ComplexSpecies)
        #Polymer Complex
        testpoly = Complex([d[1],a]).parent
        truth = OrderedPolymerSpecies([a,Complex([b,a]),a])
        self.assertEqual(testpoly,truth)
        #Ordered Complex
        truth = OrderedComplexSpecies([a,Complex([b,a]),a])
        testcomplx = Complex([a,Complex([b,a]),a],ordered=True)
        self.assertEqual(truth,testcomplx)
