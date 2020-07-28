#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Mixture, Species, DNA, Reaction, ChemicalReactionNetwork, Component, SimpleTranscription, SimpleTranslation


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
        mixture.add_component(component)

        # test that the new components was added to the mixture
        self.assertTrue(mixture.get_component(component = component) is not None)
        #test that it was added by copying
        self.assertTrue(component not in mixture._components)
        #test that the same component cannot be added again
        with self.assertRaisesRegexp(ValueError,f'{component} of the same type and name already in Mixture!'):
            mixture.add_component(component)


        #use the constructor the other way
        mixture = Mixture(components = [component])
        # test that the new components was added to the mixture
        self.assertTrue(mixture.get_component(component = component) is not None)
        #test that it was added by copying
        self.assertTrue(component not in mixture._components)


        species = Species('test_species')
        # species are invalid components
        with self.assertRaisesRegexp(ValueError,f'add_components expected a list of Components.'):
            mixture.add_components(species)

    def test_get_component(self):
        c1 = Component('c1')
        c2 = Component('c2')
        c3 = DNA('c1')
        c4 = DNA("c2")
        mixture = Mixture()
        mixture.add_components([c1, c2, c3])

        #Test with name
        #nothing named c3 is in the mixture
        self.assertTrue(mixture.get_component(name = "c3") is None)
        #two things named c2 are in the mixture
        self.assertTrue(type(mixture.get_component(name = "c1")) is list)
        #one thing named c1 is in the mixture
        self.assertTrue(type(mixture.get_component(name = "c2")) is Component)

        #test with component
        self.assertTrue(type(mixture.get_component(component = c1)) is Component)
        self.assertTrue(type(mixture.get_component(component = c2)) is Component)
        self.assertTrue(type(mixture.get_component(component = c3)) is DNA)
        self.assertTrue(mixture.get_component(component = c4) is None)

        #Test with index
        self.assertTrue(type(mixture.get_component(index = 1)) is Component)
        self.assertTrue(type(mixture.get_component(index = 2)) is DNA)
        with self.assertRaisesRegexp(IndexError,""):
            mixture.get_component(index = 3)

        with self.assertRaisesRegexp(ValueError,""):
            mixture.get_component(index = 3, component = c1)

        with self.assertRaisesRegexp(ValueError,""):
            mixture.get_component(component = "c1")

        with self.assertRaisesRegexp(ValueError,""):
            mixture.get_component(name = c1)

        with self.assertRaisesRegexp(ValueError,""):
            mixture.get_component(index = c1)


    def test_update_mechanisms(self):

        mixture = Mixture()

        tx = SimpleTranscription()
        tl = SimpleTranslation()

        test_mech = {tx.mechanism_type: tx, tl.mechanism_type: tl}

        # test that mixture has no mechanism
        self.assertTrue(isinstance(mixture.mechanisms, dict) and len(mixture.mechanisms) == 0)

        #test mechanism setter
        mixture.mechanisms = test_mech
        self.assertEqual(mixture.mechanisms.keys(), test_mech.keys())

        #test mechanisms are copied
        self.assertEqual(type(mixture.mechanisms[tx.mechanism_type]),  type(test_mech[tx.mechanism_type]))
        self.assertFalse(mixture.mechanisms[tx.mechanism_type] == test_mech[tx.mechanism_type])

        #remove all mechanisms
        mixture.mechanisms = {}
        self.assertEqual(mixture.mechanisms, {})

        #test add_mechanism
        mixture.add_mechanism(tx, tx.mechanism_type)
        mixture.add_mechanism(tl)
        self.assertEqual(mixture.mechanisms.keys(), test_mech.keys())

        #test add_mechanisms with list
        mixture.mechanisms = {}
        test_mech_list = list(test_mech.values())
        mixture.add_mechanisms(test_mech_list)
        self.assertEqual(mixture.mechanisms.keys(), test_mech.keys())


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
            rxn = Reaction.from_massaction(inputs=[], outputs=[], k_forward=0.1)
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

        rxn = Reaction.from_massaction(inputs=[a], outputs=[b], k_forward=0.1)

        CRN = ChemicalReactionNetwork(species_list, [rxn])

        # create a component
        component = Component("comp")

        # creating a mock update function to decouple the update process from the rest of the code
        def mock_update_reactions():
            rxn = Reaction.from_massaction(inputs=[a], outputs=[b], k_forward=0.1)
            return [rxn]

        def mock_update_species():
            return [a, b]

        component.update_species = mock_update_species
        component.update_reactions = mock_update_reactions

        mixture = Mixture(components=[component])

        crn_from_mixture = mixture.compile_crn()
        # test that the mixture has the same species as the manually build CRN object
        self.assertEqual(CRN.species, crn_from_mixture.species)
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture.reactions)
