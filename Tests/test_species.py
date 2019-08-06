from unittest import TestCase


class TestSpecies(TestCase):

    def test_species_initialization(self):
        from biocrnpyler import Species

        with self.assertWarns(Warning):
            Species(name='test_species', material_type='complex')

        species = Species(name='test_species')
        self.assertTrue(isinstance(species.attributes, list))

        species = Species(name='test_species', attributes=[None, None, None])

        self.assertTrue(None not in species.attributes)

    def test_add_attribute(self):
        from biocrnpyler import Species

        species = Species(name='test_species')

        with self.assertRaises(ValueError):
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
