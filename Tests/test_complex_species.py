#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from unittest import TestCase
from biocrnpyler import Species
from biocrnpyler import ComplexSpecies
from biocrnpyler import OrderedComplexSpecies, OrderedPolymerSpecies

"""This file tests ComplexSpecies, OrderedComplexSpecies which are all subclasses of species."""


class TestComplexSpecies(TestCase):

    def test_species_initialization(self):

        s1 = Species(name='s1')
        s2 = Species(name='s2', material_type="m2")

        # Check invalidity of ComplexSpecies with fewer than 2 component
        with self.assertRaisesRegex(ValueError,'chemical_reaction_network.complex requires '
                                                '2 or more species in its constructor'):
            ComplexSpecies([s1], called_from_complex = True)
        
        # Check invalidity of OrderedComplexSpecies with fewer than 2 component
        with self.assertRaisesRegex(ValueError, 'chemical_reaction_network.complex requires 2 '
                                                 'or more species in its constructor.'):
            OrderedComplexSpecies([s1], called_from_complex = True)
        

        # Check the naming conventions
        oc1 = OrderedComplexSpecies([s2, s1], called_from_complex = True)
        c1 = ComplexSpecies([s2, s1], called_from_complex = True)
        c3 = ComplexSpecies([s1, s1], called_from_complex = True)

        # Check invalidity of ComplexSpecies and OrderedComplexSpecies with strings instead of species
        with self.assertRaisesRegex(TypeError, 'recieved a non-species as a member of the list species'):
            self.assertEqual(OrderedComplexSpecies([s2, "s1"]), oc1)
        with self.assertRaisesRegex(TypeError, 'recieved a non-species as a member of the list species'):
            self.assertEqual(ComplexSpecies([s2, "s1"]), c1)
            
        
        # ComplexSpecies should sort the species added alphabetically by representation
        #an extra underscore is added at the end of the representation to deal with embedded complexes.
        self.assertEqual(repr(c1), "complex_"+repr(s2)+"_"+repr(s1)+"_")
        
        # OrderedComplexSpecies do not sort their internal species
        self.assertEqual(repr(oc1), "ordered_complex_"+repr(s2)+"_"+repr(s1)+"_")

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


    def test_contains(self):
        s1 = Species(name='s1', material_type="m1")
        s2 = Species(name='s2', material_type="m2")
        s3 = Species("s3")

        c1 = ComplexSpecies([s1, s2, s1], called_from_complex = True)
        c2 = ComplexSpecies([c1, s3], called_from_complex = True)

        self.assertTrue(s1 in c1)
        self.assertTrue(s2 in c1)
        self.assertFalse(c1 in c1)
        self.assertTrue(s1 in c2)
        self.assertTrue(s2 in c2)
        self.assertTrue(c1 in c2)
        self.assertTrue(s3 in c2)

        #contains checks the direction, position, and parent as well.
        p = OrderedPolymerSpecies([s1, [s2, "reverse"]])
        p2 = OrderedPolymerSpecies([s2, [s2, "reverse"]])
        c3 = ComplexSpecies([p[0], p[1]], called_from_complex = True)
        assert not s1 in c3
        assert not s2 in c3
        assert p[0] in c3
        assert p[1] in c3
        assert not p2[1] in c3
        assert p[0] not in c1

    def test_contains_species_monomer(self):
        """Checks if the ComplexSpecies has a monomer (Species) inside of it, 
        but without checking Species.parent or Species.position. In effect, a
        less stringent version of __contains__. This is useful for checking
        complexes containing monomers from Polymers."""

        s1 = Species(name='s1', material_type="m1")
        s2 = Species(name='s2', material_type="m2")
        p = OrderedPolymerSpecies([s1, [s2, "reverse"]])
        p2 = OrderedPolymerSpecies([[s2, "reverse"], s1])
        c2 = ComplexSpecies([p[0], p[1]], called_from_complex = True)

        #Contains parent doesn't matter
        assert c2.contains_species_monomer(s1)
        assert c2.contains_species_monomer(p[0])
        assert c2.contains_species_monomer(p2[0])

        #position doesn't matter
        assert c2.contains_species_monomer(p2[1])

        #Direction does matter (by default)
        assert not c2.contains_species_monomer(s2)

        #Direction does matter (with keyword)
        print("*****")
        print(s2 is p2[0], s2, p2[0])        
