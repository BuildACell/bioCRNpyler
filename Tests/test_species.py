#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Species


class TestSpecies(TestCase):

    def test_species_initialization(self):
        # species should not have type 'complex' should use ComplexSpecies class
        with self.assertWarnsRegex(Warning, f'species which are formed of two species or more should be called using'):
            Species(name='test_species', material_type='complex')

        # tests naming convention repr without species type or attributes
        species = Species(name='test_species')
        self.assertEqual(repr(species), species.name)

        # tests material type
        species = Species(name='test_species', material_type="dna")
        self.assertTrue(species.material_type == "dna")

        # tests emtpy attributes
        self.assertTrue(isinstance(species.attributes, list))

        # tests naming convention via repr without attributes
        self.assertEqual(repr(species), species.material_type + "_" + species.name)

        # tests adding attributes
        attr_list = ['atr1', 'atr2']
        species = Species(name='test_species', attributes=attr_list)
        self.assertEqual(attr_list, species.attributes)

        # tests naming convention with attributes and no material
        correct_name = species.name
        for attribute in species.attributes:
            correct_name += "_"+attribute
        self.assertEqual(repr(species), correct_name)

        # tests initial condition by default should be 0
        species = Species(name='test_species')
        self.assertTrue(species.initial_concentration == 0)

        # tests setting correct initial concentration
        initial_concentration = 10
        species = Species(name='test_species', initial_concentration=initial_concentration)
        self.assertEqual(species.initial_concentration, initial_concentration)

    def test_add_attribute(self):
        species = Species(name='test_species')
        # an attribute must be a string
        with self.assertRaisesRegexp(AssertionError,f'must be an alpha-numeric string'):
            species.add_attribute({'k': 'v'})

        species.add_attribute('attribute')
        # testing whether a valid attribute has been added to the attribute list
        self.assertTrue('attribute' in species.attributes)

    def test_species_equality(self):
        # testing species equality
        s1 = Species(name='a', material_type='mat1', attributes=['red'])
        s2 = Species(name='a', material_type='mat1', attributes=['red'])
        # if two species have the same name, material_type and attributes, then they are the same
        self.assertTrue(s1 == s2)

        s3 = Species(name='b', material_type='mat1', attributes=['red'])
        # different species name: not the same species
        self.assertFalse(s1 == s3)

        s4 = Species(name='a', material_type='mat2', attributes=['red'])
        # different material type: not the same species
        self.assertFalse(s1 == s4)

        s5 = Species(name='a', material_type='mat1', attributes=['red', 'large'])
        # different attributes: not the same species
        self.assertFalse(s1 == s5)
