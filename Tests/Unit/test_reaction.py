
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Reaction, Species, MassAction, WeightedSpecies
import pytest


class TestReaction(TestCase):

    def test_reaction_initialization(self):
        # warns if both input and output species are empty
        mak = MassAction(k_forward=0.1)
        with self.assertWarns(Warning):
            Reaction(inputs=[], outputs=[], propensity_type=mak)

        # test for invalid propensity type
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], propensity_type=Species)

        # input must be a valid species object
        with self.assertRaises(TypeError):
            Reaction(inputs=['a'], outputs=[], propensity_type=mak)
        # output must be a valid species object
        with self.assertRaises(TypeError):
            Reaction(inputs=[], outputs=['b'], propensity_type=mak)

        rxn = Reaction.from_massaction(inputs=[], outputs=[], k_forward=0.1, k_reverse=1)
        # test whether the reaction is registered as reversible
        self.assertTrue(rxn.is_reversible)
        # test whether the reaction is registered as massaction
        self.assertTrue(isinstance(rxn.propensity_type, MassAction))

        # test WeightedSpecies inputs
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')
        chem_com_sp1 = WeightedSpecies(species=sp1, stoichiometry=2)
        chem_com_sp2 = WeightedSpecies(species=sp2, stoichiometry=1)
        Reaction(inputs=[chem_com_sp1], outputs=[chem_com_sp2], propensity_type=MassAction(k_forward=1))

        # test different input and output lists
        Reaction(inputs=[chem_com_sp1], outputs=[sp2], propensity_type=MassAction(k_forward=1))

        # mixing WeightedSpecies and Species is not allowed
        with self.assertRaises(TypeError):
            Reaction(inputs=[chem_com_sp1, sp2], outputs=[sp1], propensity_type=MassAction(k_forward=1))


def test_species_merging():
    sp1 = Species(name='test_species_a')
    chem_complexes = [WeightedSpecies(species=sp1, stoichiometry=2),
                      WeightedSpecies(species=sp1, stoichiometry=1)]
    rxn = Reaction(inputs=chem_complexes, outputs=[], propensity_type=MassAction(k_forward=1))
    # same species with different stoichiometry gets merged into one species
    assert len(rxn.inputs) == 1

    sp1 = Species(name='test_species_a')
    chem_complexes = [WeightedSpecies(species=sp1, stoichiometry=2),
                      WeightedSpecies(species=sp1, stoichiometry=1)]
    rxn = Reaction(inputs=[], outputs=chem_complexes, propensity_type=MassAction(k_forward=1))

    # same species with different stoichiometry gets merged into one species
    assert len(rxn.outputs) == 1


def test_reaction_equality():
    """test for the_equality operator"""
    sp1 = Species(name='test_species_a')
    sp2 = Species(name='test_species_b')
    rxn1 = Reaction(inputs=[sp1, sp2], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=1))
    rxn2 = Reaction(inputs=[sp2, sp1], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=1))
    rxn3 = Reaction(inputs=[sp2, sp1], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=10))
    rxn4 = Reaction(inputs=[sp2, sp1], outputs=[sp2], propensity_type=MassAction(k_forward=1))
    assert rxn1 == rxn2
    assert rxn1 != rxn3
    assert rxn1 != rxn4


def test_reaction_list_flattening():
    sp1 = Species(name='test_species_a')
    sp2 = Species(name='test_species_b')
    k_f = 1
    mak = MassAction(k_forward=k_f)
    rxn1 = Reaction(inputs=[sp1, [sp1, sp2]], outputs=[[sp2, sp2], sp1], propensity_type=mak)
    rxn2 = Reaction(inputs=[sp1, sp1, sp2], outputs=[sp1, sp2, sp2], propensity_type=mak)
    assert rxn1 == rxn2
