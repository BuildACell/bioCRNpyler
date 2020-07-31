#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from unittest import TestCase
from biocrnpyler import Species
from biocrnpyler import ComplexSpecies
from biocrnpyler import OrderedComplexSpecies
from biocrnpyler import Multimer

"""This file tests ComplexSpecies, OrderedComplexSpecies, and Multimers which are all subclasses of species."""


class TestComplexSpecies(TestCase):

    def test_species_initialization(self):

        s1 = Species(name='s1')
        s2 = Species(name='s2', material_type="m2")

        # Check invalidity of ComplexSpecies with fewer than 2 component
        with self.assertRaisesRegexp(ValueError,'chemical_reaction_network.complex requires '
                                                '2 or more species in its constructor'):
            ComplexSpecies([s1], called_from_complex = True)
        
        # Check invalidity of OrderedComplexSpecies with fewer than 2 component
        with self.assertRaisesRegexp(ValueError, 'chemical_reaction_network.complex requires 2 '
                                                 'or more species in its constructor.'):
            OrderedComplexSpecies([s1], called_from_complex = True)
        
        # Check invalidity of multimers with fewer than 2 component
        with self.assertRaisesRegexp(ValueError, 'chemical_reaction_network.complex requires 2 '
                                                 'or more species in its constructor'):
            Multimer(s1, 1, called_from_complex = True)

        # Check the naming conventions
        oc1 = OrderedComplexSpecies([s2, s1], called_from_complex = True)
        c1 = ComplexSpecies([s2, s1], called_from_complex = True)
        m1 = Multimer(s1, 2, called_from_complex = True)
        c3 = ComplexSpecies([s1, s1], called_from_complex = True)

        # Check invalidity of ComplexSpecies, Multimers and OrderedComplexSpecies with strings instead of species
        with self.assertRaisesRegexp(TypeError, 'recieved a non-species as a member of the list species'):
            self.assertEqual(OrderedComplexSpecies([s2, "s1"]), oc1)
            self.assertEqual(ComplexSpecies([s2, "s1"]), c1)
            self.assertEqual(Multimer("s1", 2), m1)
            
        
        # ComplexSpecies should sort the species added alphabetically by representation
        #an extra underscore is added at the end of the representation to deal with embedded complexes.
        self.assertEqual(repr(c1), "complex_"+repr(s2)+"_"+repr(s1)+"_")
        
        # OrderedComplexSpecies do not sort their internal species
        self.assertEqual(repr(oc1), "ordered_complex_"+repr(s2)+"_"+repr(s1)+"_")
        
        # Multimers are just complexes with multiplicity
        self.assertEqual(repr(m1), "complex_"+repr(s1)+"_2x_")
        self.assertEqual(repr(c3), repr(m1))

        # Nested list creation of ComplexSpecies
        c1 = ComplexSpecies([s1, [s2, s1]], called_from_complex = True)
        c2 = ComplexSpecies([s1, s2, s1], called_from_complex = True)
        self.assertEqual(c1, c2)

        c1 = OrderedComplexSpecies([s1, [s2, s1]], called_from_complex = True)
        c2 = OrderedComplexSpecies([s1, s2, s1], called_from_complex = True)
        c3 = OrderedComplexSpecies([s1, [s1, s2]], called_from_complex = True)
        self.assertEqual(c1, c2)
        self.assertFalse(c1==c3)

    def test_species_equality(self):

        s1 = Species(name='s1', material_type="m1")
        s2 = Species(name='s2', material_type="m2")

        c1 = ComplexSpecies([s1, s2], called_from_complex = True)
        c2 = ComplexSpecies([s2, s1], called_from_complex = True)
        # Check equality of differently ordered complexes
        self.assertEqual(c1, c2)
                         
        oc1 = OrderedComplexSpecies([s2, s1], called_from_complex = True)
        oc2 = OrderedComplexSpecies([s1, s2], called_from_complex = True)
        self.assertFalse(oc1 == oc2)
        #check of __contains__
        self.assertTrue(s1 in c1)
