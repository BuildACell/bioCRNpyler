# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Species, Compartment
import pytest


class TestSpecies(TestCase):
    def test_compartment_initialization(self):
        # tests naming convention repr without species type or attributes
        species = Species(name='test_species')
        compartment = Compartment(name='test_compartment')
        species.compartment = compartment
        self.assertEqual(compartment.name, 'test_compartment')
        self.assertEqual(species.compartment.name, compartment.name)

        # Test non-string named compartment
        with self.assertRaisesRegex(ValueError, 'Compartment name must be a string.'):
            compartment = Compartment(name = 24)

        with self.assertRaisesRegex(TypeError, 'Compartment name must be a string.'):
            compartment = Compartment(name = "test_compartment")
            compartment.name = None

        # tests naming convention for species with name and compartment
        compartment = Compartment(
            name='test_compartment', spatial_dimensions=1, size=1e-4)
        self.assertEqual(compartment.spatial_dimensions, 1)
        self.assertEqual(compartment.size, 1e-4)
        with self.assertRaises(ValueError):
            Compartment(name='2test_compartment')

        with self.assertRaisesRegex(ValueError, 'Compartment name must be a string'):
            compartment = Compartment(name="test_compartment")
            compartment.name = 24

        with self.assertRaisesRegex(ValueError, 'Compartment size must be a float or int.'):
            Compartment(name="test_compartment", size='2.5')

        with self.assertRaisesRegex(ValueError, 'Compartment size must be non-negative.'):
            Compartment(name="test_compartment", size=-2)

        with self.assertRaisesRegex(ValueError, 'Compartment spatial dimension must be an integer.'):
            Compartment(name="test_compartment", spatial_dimensions=2.5)

        with self.assertRaisesRegex(ValueError, 'Compartment spatial dimension must be non-negative.'):
            Compartment(name="test_compartment", spatial_dimensions=-2)

        with self.assertRaisesRegex(ValueError, 'Compartments with same names must have the same size and spatial dimensions.'):
            compartment1 = Compartment(
                name="test_compartment", spatial_dimensions=2)
            compartment2 = Compartment(
                name="test_compartment", spatial_dimensions=3)
            self.assertEqual(compartment1, compartment2)
        compartment = Compartment(name="test_compartment", unit="litre")
        self.assertEqual(compartment.unit, "litre")
        compartment.unit = "nanolitre"
        self.assertEqual(compartment.unit, "nanolitre")
        with self.assertRaisesRegex(ValueError, 'Unit of compartment must be a string representing compartment size.'):
            compartment = Compartment(name="test_compartment")
            compartment.unit = 24

        compartment = Compartment(name="test_compartment")
        self.assertEqual(compartment.unit, None)
        
    def test_add_compartment(self):
        c1 = Compartment("c1")
        c2 = Compartment("c2")
        c1.add_compartment("internal", c2)
        
        self.assertTrue(c2 in c1.compartment_dict["internal"])
        
        # testing errors
       
        with self.assertRaisesRegex(ValueError,'Relationship must be a string!'):
            c1.add_compartment(None, c2)
            
        with self.assertRaisesRegex(ValueError,"You did not input a valid Compartment or list of Compartment objects"):
            c1.add_compartment("random", "I am not a compartment")
            
        with self.assertRaisesRegex(ValueError,"You did not input a valid Compartment or list of Compartment objects"):
            c1.add_compartment("internal", "I am not a compartment")
        
    def test_set_compartment(self):
        c1 = Compartment("c1")
        c2 = Compartment("c2")
        c3 = Compartment("c3")
        c1.add_compartment("internal", c2)
        c1.add_compartment("internal", c3)
        
        c1.set_compartment("internal", c3)
      
        self.assertTrue(c3 in c1.get_compartment("internal"))
        self.assertTrue(1 == len( c1.get_compartment("internal")))
        
        with self.assertRaisesRegex(ValueError,"You did not input a valid Compartment or list of Compartment objects"):
            c1.set_compartment("internal", "external")
        with self.assertRaisesRegex(ValueError,'Relationship must be a string!'):
            c1.set_compartment(None, c2)
    def test_get_compartment(self):
        c1 = Compartment("c1")
        c2 = Compartment("c2")
        c3 = Compartment("c3")
        c1.add_compartment("internal", c2)
        c1.add_compartment("internal", c3)
        
        result = c1.get_compartment("internal")
        self.assertTrue(c2 in c1.get_compartment("internal") and c3 in c1.get_compartment("internal") and len(c1.get_compartment("internal")) == 2)
        
        with self.assertRaisesRegex(KeyError,"No compartment exists by that name!"):
            c1.get_compartment("external")

    def test_get_compartment_dict(self):
        c1_new = Compartment("c1")
        c2 = Compartment("c2")
        c1_new.add_compartment("internal", c2)
    
        self.assertTrue(len((c1_new.get_compartment_dict()).keys()) is 1)
        self.assertTrue("internal" in c1_new.get_compartment_dict())
        self.assertTrue(c2 in( c1_new.get_compartment_dict())["internal"] )
