import unittest

from biocrnpyler import *


class DNAassemblyTest(unittest.TestCase):
    def test_empty_assembly(self):
        name = "test_assembly"
        construct = DNAassembly(name=name)

        self.assertTrue(construct.name == name)
        self.assertTrue(isinstance(construct.dna, Specie))
        self.assertTrue(construct.rbs is None)
        self.assertTrue(construct.promoter is None)
        self.assertTrue(isinstance(construct.transcript, Specie))
        self.assertTrue(isinstance(construct.protein, Specie))


