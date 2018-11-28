import unittest

from biocrnpyler import *


class DNATest(unittest.TestCase):
    def test_DNA(self):
        name = "test_DNA"
        dna = DNA(name=name)

        self.assertTrue(dna.name == name)
        self.assertTrue(dna.length == 0)

        with self.assertRaises(ValueError) as context:
            dna.length = -1

            self.assertTrue('non-negative length enforced' in context.exception)

        self.assertTrue(isinstance(dna.mechanisms, dict) and len(dna.mechanisms) == 0)
        self.assertTrue(isinstance(dna.parameters, dict) and len(dna.parameters) == 0)
        self.assertTrue(isinstance(dna.attributes, list) and len(dna.attributes) == 0)
