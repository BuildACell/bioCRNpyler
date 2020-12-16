#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from unittest import TestCase
from biocrnpyler import Species, Mixture, Component, ParameterDatabase, ChemicalReactionNetwork

class TestOrderedMonomer(TestCase):

    def test_species_initial_condition_in_crn(self):
        S = Species("S", initial_concentration = 1.1)
        CRN = ChemicalReactionNetwork(species = [S], reactions = [])
        self.assertTrue(CRN.species[0].initial_concentration == 1.1)



    def test_species_initial_condition_in_mixture(self):
        S = Species("S", initial_concentration = 1.0)
        M = Mixture(species = [S])
        CRN = M.compile_crn()
        print("S", CRN.species[0])
        print("initial concentration", CRN.species[0].initial_concentration)
        self.assertTrue(CRN.species[0].initial_concentration == 1.1)



	

