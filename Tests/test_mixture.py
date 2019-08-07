
#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase


class TestMixture(TestCase):
    def test_add_components(self):
        from biocrnpyler import Mixture
        from biocrnpyler import Component
        # from biocrnpyler import Species

        mixture = Mixture()
        self.assertTrue(len(mixture.components) == 0)
        component = Component('test_comp')
        mixture.add_components(component)

        self.assertTrue(component in mixture.components)

        # species = Species('test_species')
        # mixture.add_components(species)
        #
        # self.assertTrue(species in mixture.components)

    def test_update_species(self):
        from biocrnpyler import Mixture
        from biocrnpyler import Species
        from biocrnpyler import DNA
        # from biocrnpyler import Dilution

        species = Species(name='H2O')
        mixture = Mixture(components=[species])

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
        from biocrnpyler import Mixture
        from biocrnpyler import Reaction
        from biocrnpyler import Component

        mixture = Mixture()
        with self.assertRaises(AttributeError):
            mixture.update_reactions()

        component = Component(name='test_component')

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
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction
        from biocrnpyler import Mixture

        a = Species(name='a')
        b = Species(name='b')

        species_list = [a, b]

        def mock_update_reactions():
            rxn = Reaction(inputs=[a], outputs=[b], k=0.1)
            return [rxn]

        rxn = Reaction(inputs=[a], outputs=[b], k=0.1)

        CRN = ChemicalReactionNetwork(species_list, [rxn])

        mixture = Mixture(components=species_list)
        mixture.update_reactions = mock_update_reactions

        crn_from_mixture = mixture.compile_crn()
        self.assertEqual(CRN.species, crn_from_mixture.species)
        self.assertEqual(CRN.reactions, crn_from_mixture.reactions)


