#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase


class TestSpecies(TestCase):

    def test_species_initialization(self):
        from biocrnpyler import Species

        with self.assertWarns(Warning):
            Species(name='test_species', material_type='complex')

        species = Species(name='test_species')
        self.assertTrue(isinstance(species.attributes, list))

        attr_list = ['atr1', 'atr2']
        species = Species(name='test_species', attributes=attr_list)
        self.assertEqual(attr_list, species.attributes)

    def test_add_attribute(self):
        from biocrnpyler import Species

        species = Species(name='test_species')

        with self.assertRaises(AssertionError):
            species.add_attribute({'k': 'v'})

        species.add_attribute('attribute')

    def test_species_equality(self):
        from biocrnpyler import Species

        s1 = Species(name='a', material_type='mat1', attributes=['red'])
        s2 = Species(name='a', material_type='mat1', attributes=['red'])

        self.assertTrue(s1 == s2)

        s3 = Species(name='b', material_type='mat1', attributes=['red'])

        self.assertFalse(s1 == s3)

        s4 = Species(name='a', material_type='mat2', attributes=['red'])

        self.assertFalse(s1 == s4)

        s5 = Species(name='a', material_type='mat1', attributes=['red', 'large'])

        self.assertFalse(s1 == s5)
