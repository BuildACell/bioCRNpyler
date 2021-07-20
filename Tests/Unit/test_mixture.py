#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Mixture, Species, DNA, Reaction, ChemicalReactionNetwork, Component, SimpleTranscription, SimpleTranslation, GlobalMechanism


class TestMixture(TestCase):

    def test_add_species(self):
        species = Species('test_species')
        mixture = Mixture()
        mixture.add_species(species)

        # test that the new mixture only contain one species
        self.assertEqual([species], mixture.added_species)

        # adding invalid species
        with self.assertRaisesRegex(AssertionError,'only Species type is accepted!'):
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
        with self.assertRaisesRegex(ValueError,f'{component} of the same type and name already in Mixture!'):
            mixture.add_component(component)


        #use the constructor the other way
        mixture = Mixture(components = [component])
        # test that the new components was added to the mixture
        self.assertTrue(mixture.get_component(component = component) is not None)
        #test that it was added by copying
        self.assertTrue(component not in mixture._components)


        species = Species('test_species')
        # species are invalid components
        with self.assertRaisesRegex(ValueError,f'add_components expected a list of Components.'):
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
        with self.assertRaisesRegex(IndexError,""):
            mixture.get_component(index = 3)

        with self.assertRaisesRegex(ValueError,""):
            mixture.get_component(index = 3, component = c1)

        with self.assertRaisesRegex(ValueError,""):
            mixture.get_component(component = "c1")

        with self.assertRaisesRegex(ValueError,""):
            mixture.get_component(name = c1)

        with self.assertRaisesRegex(ValueError,""):
            mixture.get_component(index = c1)


    def test_add_mechanisms(self):

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

    def test_add_global_mechanism(self):
        mixture = Mixture()
        GM = GlobalMechanism(name = "gm", mechanism_type = "global")

        mixture.add_global_mechanism(GM)
        #Global Mechanisms should be copied!
        self.assertTrue(mixture.global_mechanisms["global"] != GM)
        self.assertTrue(mixture.global_mechanisms["global"].name == GM.name)

        #test setter to clear
        mixture.global_mechanisms = {}
        self.assertTrue("global" not in mixture.global_mechanisms)
        #test add via add_mechanisms with different key
        mixture.add_mechanisms({"key2":GM})
        self.assertTrue("global" not in mixture.global_mechanisms)
        self.assertTrue("key2" in mixture.global_mechanisms)

    def test_add_species_to_crn(self):
        species = Species(name='H2O')
        mixture = Mixture(species=[species])
        mixture.add_species_to_crn(species, None)

        self.assertTrue(species in mixture.crn.species)

        dna = DNA(name='test_DNA')
        mixture.add_components(dna)

        mixture.add_species_to_crn(dna.update_species(), dna)

        for s_dna in dna.update_species():
            self.assertTrue(s_dna in mixture.crn.species)

        # Currently, there is no global mechanism that creates new species

        # dilution_mechanism = Dilution()
        # global_mechanisms = {"dilution": dilution_mechanism}
        #
        # mixture = Mixture(global_mechanisms=global_mechanisms)
        # mixture.update_species()

    """
    update mechanisms has been removed - keeping here just in case it comes back
    def test_update_reactions(self):
        mixture = Mixture()
        with self.assertRaisesRegex(AttributeError, 'Mixture.crn_species not defined.'):
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
        self.assertEqual(crn_rxn, crn_rxn_mock)"""

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
        self.assertEqual(set(CRN.species), set(crn_from_mixture.species))
        #test that species are copied
        self.assertTrue(all([not s is ss for s in CRN.species for ss in crn_from_mixture.species]))
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture.reactions)
        #test that species are copied
        self.assertTrue(all([not r is rr for r in CRN.reactions for rr in crn_from_mixture.reactions]))

    def test_compile_crn_directives(self):
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

        #All the directives used below should not change the CRN. 
        #They just remove some checks and safegaurds, but compilation should work the same in this simple case.
        #directives are best used in specific cases to compile very large models where speed is essential

        crn_from_mixture1 = mixture.compile_crn(copy_objects = False)
        # test that the mixture has the same species as the manually build CRN object
        self.assertEqual(set(CRN.species), set(crn_from_mixture1.species))
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture1.reactions)


        crn_from_mixture2 = mixture.compile_crn(add_reaction_species = False)
        # test that the mixture has the same species as the manually build CRN object
        self.assertEqual(set(CRN.species), set(crn_from_mixture2.species))
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture2.reactions)

        crn_from_mixture3 = mixture.compile_crn(initial_concentrations_at_end = True)
        # test that the mixture has the same species as the manually build CRN object
        self.assertEqual(set(CRN.species), set(crn_from_mixture3.species))
        # test that the mixture has the same reactions as the manually build CRN object
        self.assertEqual(CRN.reactions, crn_from_mixture3.reactions)



    def test_compoents_in_multiple_mixtures(self):
        C = Component("comp")
        M1 = Mixture(components = [C])

        #M1 should have a copy of C, not C itself
        self.assertTrue(type(M1.get_component(component = C)) == Component)
        self.assertTrue(C not in M1._components)

        #C.mixture should not be set, only it's copy
        self.assertTrue(C.mixture is None)
        self.assertTrue(M1.get_component(component = C).mixture is M1)

        #Add C to another Mixture
        M2 = Mixture(components = [C])
        self.assertTrue(type(M2.get_component(component = C))  == Component)
        self.assertTrue(M2.get_component(component = C).mixture is M2)
        self.assertTrue(M1.get_component(component = C) != M2.get_component(component = C))
