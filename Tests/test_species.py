#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Species, WeightedSpecies
import pytest


class TestSpecies(TestCase):

    def test_species_initialization(self):
        # tests naming convention repr without species type or attributes
        species = Species(name='test_species')
        self.assertEqual(repr(species), species.name)

        # tests naming convention for species with name and compartment
        species = Species(name='test_species', compartment = 'test_compartment')
        self.assertEqual(repr(species), species.name)
        self.assertEqual(species.compartment, 'test_compartment')

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

        #test OrderedMonomer subclass
        self.assertTrue(species.parent is None)
        self.assertTrue(species.position is None)
        self.assertTrue(species.direction is None)

    def test_add_attribute(self):
        species = Species(name='test_species')
        # an attribute must be a string
        with self.assertRaisesRegex(AssertionError,f'must be an alpha-numeric string'):
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

        s6 = Species(name='a', material_type='mat1', attributes=['red'], compartment = 'test_compartment')
        # same species name in different compartments: not the same species
        self.assertFalse(s1 == s6)


def test_weighted_species_init():
    s1 = Species(name='a')
    # normal operations
    ws1 = WeightedSpecies(species=s1)
    assert ws1.species == s1
    assert ws1.stoichiometry == 1

    with pytest.raises(ValueError, match='Stoichiometry must be positive integer!'):
        WeightedSpecies(species=s1, stoichiometry=0)

    with pytest.raises(ValueError, match='Stoichiometry must be positive integer!'):
        WeightedSpecies(species=s1, stoichiometry=-1)

    ws2 = WeightedSpecies(species=s1, stoichiometry=1.34)
    assert ws2.stoichiometry == 1


def test_merging_weighted_species():
    s1 = Species(name='a')

    ws1 = WeightedSpecies(species=s1, stoichiometry=2)
    ws2 = WeightedSpecies(species=s1, stoichiometry=5)
    ws_list = [ws1, ws2]

    freq_dict = WeightedSpecies._count_weighted_species(ws_list)
    assert len(freq_dict) == 1
    ws_merged = list(freq_dict.values())
    assert ws_merged[0] == 7

    s2 = Species(name='b')

    ws1 = WeightedSpecies(species=s1, stoichiometry=2)
    ws2 = WeightedSpecies(species=s2, stoichiometry=5)
    ws_list = [ws1, ws2]
    freq_dict = WeightedSpecies._count_weighted_species(ws_list)
    assert len(freq_dict) == 2

    ws_merged = list(freq_dict.values())
    assert ws_merged[0] == 2
    assert ws_merged[1] == 5


def test_weighted_species_equality():
    s1 = Species(name='a')
    ws1 = WeightedSpecies(species=s1, stoichiometry=1)

    s2 = Species(name='a')
    ws2 = WeightedSpecies(species=s2, stoichiometry=3)
    assert ws1 != ws2

    s2 = Species(name='b')
    ws3 = WeightedSpecies(species=s2, stoichiometry=1)

    assert ws1 != ws3













