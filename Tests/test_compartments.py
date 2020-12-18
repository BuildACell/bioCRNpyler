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

        # tests naming convention for species with name and compartment
        compartment = Compartment(name='test_compartment', spatial_dimensions = 1, volume = 1e-4)
        self.assertEqual(compartment.spatial_dimensions, 1)
        self.assertEqual(compartment.volume, 1e-4)
        with pytest.raises(ValueError):
            compartment = Compartment(name='2test_compartment')