#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from unittest import TestCase

#This file tests ComplexSpecies, OrderedComplexSpecies, and Multimers which are all subclasses of species.

class TestComplexSpecies(TestCase):

    def test_species_initialization(self):
        from biocrnpyler import Species
        from biocrnpyler import ComplexSpecies
        from biocrnpyler import OrderedComplexSpecies
        from biocrnpyler import Multimer
        
        
        s1 = Species(name='s1', material_type = "m1")
        s2 = Species(name = 's2', material_type = "m2")
        
        
        
        #Check invalidity of ComplexSpecies with fewer than 1 component
        with self.assertRaises(ValueError):
            c = ComplexSpecies([s1])
        
        #Check invalidity of OrderedComplexSpecies with fewer than 1 component
        with self.assertRaises(ValueError):
            oc = OrderedComplexSpecies([s1])
        
        #Check invalidity of multimers with fewer than 1 component
        with self.assertRaises(ValueError):
            m = Multimer(s1, 1)
            
        #Check invalidity of ComplexSpecies with non-species
        with self.assertRaises(ValueError):
            c = ComplexSpecies(["s1", "s2"])
        
        #Check invalidity of OrderedComplexSpecies with non-species
        with self.assertRaises(ValueError):
            oc = OrderedComplexSpecies(["s1", "s2"])
        
        #Check invalidity of multimers with non-species
        with self.assertRaises(ValueError):
            m = Multimer("s1", 2)
        
        #Check the naming conventions
        oc1 = OrderedComplexSpecies([s2, s1])
        c1 = ComplexSpecies([s2, s1])
        m1 = Multimer(s1, 2)
     
        self.assertEqual(repr(c1), "complex_"+repr(s1)+"_"+repr(s2))
        self.assertEqual(repr(oc1), "ordered_complex_"+repr(s2)+"_"+repr(s1))
        self.assertEqual(repr(m1), "multimer_2x_"+repr(s1))
        
        
        
        
    def test_species_equality(self):
        from biocrnpyler import Species
        from biocrnpyler import ComplexSpecies
        from biocrnpyler import OrderedComplexSpecies
        
        s1 = Species(name='s1', material_type = "m1")
        s2 = Species(name = 's2', material_type = "m2")
        
        c1 = ComplexSpecies([s1, s2])
        c2 = ComplexSpecies([s2, s1])
        #check equality of differently ordered complexes
        self.assertEqual(c1, c2)
                         
        oc1 = OrderedComplexSpecies([s2, s1])
        oc2 = OrderedComplexSpecies([s1, s2])
        self.assertFalse(oc1==oc2)


        