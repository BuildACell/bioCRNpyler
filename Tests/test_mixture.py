#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Mixture, Species, Component, DNA, Reaction, ChemicalReactionNetwork, Component


class TestMixture(TestCase):

    def test_add_species(self):
        species = Species('test_species')
        mixture = Mixture()
        mixture.add_species(species)

        # test that the new mixture only contain one species
        self.assertEqual([species], mixture.added_species)

        # adding invalid species
        with self.assertRaisesRegexp(AssertionError,'only Species type is accepted!'):
            mixture.add_species(['ok', 'ok'])

    def test_add_components(self):

        mixture = Mixture()
        # test that new mixture has no components
        self.assertTrue(len(mixture.components) == 0)
        component = Component('test_comp')
        mixture.add_components(component)
        # test that the new components was added to the mixture
        self.assertTrue(component in mixture.components)

        species = Species('test_species')
        # species are invalid components
        with self.assertRaisesRegexp(AssertionError,f'the object: {species} passed into mixture as component must be of the class Component'):
            mixture.add_components(species)

    def test_update_species(self):
        species = Species(name='H2O')
        mixture = Mixture(species=[species])

        self.assertTrue(species in mixture.update_species())

        dna = DNA(name='test_DNA')
        mixture.add_components(dna)

        crn_list = mixture.update_species()

        for s_dna in dna.update_species():
            self.assertTrue(s_dna in crn_list)

        # Currently, there is no global mechanism that creates new species

        # dilution_mechanism = Dilution()
        # global_mechanisms = {"dilution": dilution_mechanism}
        #
        # mixture = Mixture(global_mechanisms=global_mechanisms)
        # mixture.update_species()

    def test_update_reactions(self):
        mixture = Mixture()
        with self.assertRaisesRegexp(AttributeError, 'Mixture.crn_species not defined.'):
            mixture.update_reactions()

        component = Component(name='test_component')

        # creating a mock update function to decouple the update process from the rest of the code
        def mock_update_reactions():
            rxn = Reaction(inputs=[], outputs=[], k=0.1)
            return [rxn]

        component.update_reactions = mock_update_reactions

        mixture.add_components(component)
        mixture.update_species()
        crn_rxn = mixture.update_reactions()
        crn_rxn_mock = mock_update_reactions()
        self.assertEqual(crn_rxn, crn_rxn_mock)

        # TODO add test for reactions added by global mechanisms

    def test_compile_crn(self):
        a = Species(name='a')
        b = Species(name='b')

        species_list = [a, b]

        

        rxn = Reaction(inputs=[a], outputs=[b], k=0.1)

        CRN = ChemicalReactionNetwork(species_list, [rxn])

        #create a component
        component = Component("comp")

        # creating a mock update function to decouple the update process from the rest of the code
        def mock_update_reactions():
            rxn = Reaction(inputs=[a], outputs=[b], k=0.1)
            return [rxn]
        def mock_update_species():
            return [a, b]

        component.update_species = mock_update_species
        component.update_reactions = mock_update_reactions
        
        mixture = Mixture(components = [component])


        crn_from_mixture = mixture.compile_crn()
        # test that the mixture has the same species as the manually build CRN object
        self.assertEqual(CRN.species, crn_from_mixture.species)
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture.reactions)
