
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import mock_open, patch
from biocrnpyler import ChemicalReactionNetwork, Species, Reaction, ComplexSpecies, OrderedComplexSpecies
import libsbml
import warnings


class TestChemicalReactionNetwork(TestCase):

    def setUp(self) -> None:
        """this method gets executed before every test"""
        self.s1 = Species(name='test_species1')
        self.s2 = Species(name='test_species2')
        self.s3 = Species(name='test_species3')
        self.s4 = Species(name='test_species4')

        self.species_list = [self.s1, self.s2]
        # creating a valid reaction two species
        self.rx1 = Reaction.from_mass_action(inputs=[self.s1], outputs=[self.s2], k_forward=0.1)
        self.rxn_list = [self.rx1]

        self.crn = ChemicalReactionNetwork(species=self.species_list, reactions=self.rxn_list)

    def test_check_crn_validity(self):

        checked_reactions, checked_species = ChemicalReactionNetwork.check_crn_validity(reactions=self.rxn_list,
                                                                                        species=self.species_list)
        # test that the returned species list is the same as the species list supplied
        self.assertEqual(self.species_list, checked_species)
        # test that the returned reaction list is the same as the reaction list supplied
        self.assertEqual(self.rxn_list, checked_reactions)

        species_list_with_none = self.species_list.copy()
        # injecting a None to the species list
        species_list_with_none.append(None)
        # test whether a non-species object is detected and Value error has been raised
        #                                         A non-species object was used as a species: [test_species1, test_species2, None]!"'
        #                                         A non-species object was used as a species: [test_species1, test_species2, None]!"
        with self.assertRaisesRegexp(ValueError, "A non-species object was used as a species!"):
            ChemicalReactionNetwork.check_crn_validity(reactions=self.rxn_list, species=species_list_with_none)

        rxn_list_with_none = self.rxn_list.copy()
        # injecting a None to the reaction list
        rxn_list_with_none.append(None)
        # test whether a non-reaction object is detected and Value Error has been raised
        with self.assertRaisesRegexp(ValueError, 'A non-reaction object was used as a reaction!'):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list_with_none, species=self.species_list)

        rxn2 = Reaction.from_mass_action(inputs=[self.s1], outputs=[self.s3], k_forward=0.1)
        # test warning raised if a species (in the reaction outputs) is detected which is not part of the species list
        with self.assertWarnsRegex(Warning, f'are not part of any reactions in the CRN'):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn2], species=self.species_list, show_warnings=True)

        rxn3 = Reaction.from_mass_action(inputs=[self.s4], outputs=[self.s2], k_forward=0.1)
        # test warning raised if a species (in the reaction inputs) is detected which is not part of the species list
        with self.assertWarnsRegex(Warning, f'are not part of any reactions in the CRN'):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn3], species=self.species_list, show_warnings=True)

        # test warning if reaction has unlisted species
        rxn4 = Reaction.from_mass_action(inputs=[self.s4, self.s3], outputs=[self.s2], k_forward=0.1)
        with self.assertWarnsRegex(Warning, f'are not listed in the Species list, but part of the reactions'):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn4], species=[self.s4, self.s2], show_warnings=True)

        # test duplicate reactions are both added
        rxn_list = [self.rx1, self.rx1]

        CRN = ChemicalReactionNetwork(species = [self.s1, self.s2], reactions = rxn_list)
        self.assertTrue(CRN.reactions.count(self.rx1) == 2)

        with self.assertWarnsRegex(Warning, 'may be duplicated in CRN definitions'):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list, species=self.species_list, show_warnings=True)

        # test warning suppression
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list, species=self.species_list, show_warnings=False)

        assert not w

    def test_initial_condition_vector(self):

        # no initial value is supplied for s2
        init_cond = {self.s1: 5, self.s3: 10}

        x0 = self.crn.initial_condition_vector(init_cond_dict=init_cond)

        # test that value for s3 is not in the initial value list because s3 is not part of any reaction
        not_in_the_list = False
        for key, value in init_cond.items():
            if value not in x0 and key == self.s3:
                not_in_the_list = True
        self.assertTrue(not_in_the_list)

    def test_get_all_species_containing(self):

        # test that the species arument must be Species object
        with self.assertRaisesRegexp(ValueError, 'species argument must be an instance of Species!'):
            self.crn.get_all_species_containing(species=self.species_list)

        # s3 is not part of the CRN
        rtn_species_list = self.crn.get_all_species_containing(species=self.s3)
        # test that empty list is returned
        self.assertEqual(rtn_species_list, [])

        # test that s1 is returned only once
        rtn_species_list = self.crn.get_all_species_containing(species=self.s1)
        self.assertEqual(rtn_species_list, [self.s1])

        # test that s1 is returned only once as a string
        rtn_species_list = self.crn.get_all_species_containing(species=self.s1, return_as_strings=True)
        self.assertEqual(rtn_species_list, [repr(self.s1)])

    def test_replace_species(self):

        s_old = Species("s_old")
        s_new = Species("s_new")

        #Test replace species in a Species
        self.assertTrue(self.s1.replace_species(self.s2, s_new) == self.s1)
        self.assertTrue(s_old.replace_species(s_old, s_new) == s_new)

        #test replace Species in a ComplexSpecies
        c1 = ComplexSpecies([self.s1, s_old])
        c2 = ComplexSpecies([self.s1, c1])
        self.assertTrue(c1.replace_species(self.s2, s_new) == c1)
        self.assertTrue(c1.replace_species(s_old, s_new) == ComplexSpecies([self.s1, s_new]))
        self.assertTrue(c2 == ComplexSpecies([self.s1, c1]))
        self.assertTrue(c2.replace_species(s_new, s_old) == ComplexSpecies([self.s1, ComplexSpecies([self.s1, s_old])]))

        #test in OrderedComplexSpecies
        oc1 = OrderedComplexSpecies([self.s1, s_old])
        self.assertTrue(oc1.replace_species(s_old, s_new) == OrderedComplexSpecies([self.s1, s_new]))


        #Test replace species in a reaction
        c1 = ComplexSpecies([self.s1, s_old])
        c2 = ComplexSpecies([self.s1, s_new])
        r1 = Reaction([self.s1, s_old], [c1], k = 1)
        self.assertTrue(r1.replace_species(s_old, s_new) == Reaction([self.s1, s_new], [c2], k = 1))

        #test replace species with a non-massaction reaction
        r1 = Reaction([self.s1, s_old], [c1], propensity_type = "proportionalhillpositive", propensity_params = {"k":1., "d":s_old, "s1":self.s1, "n":2, "K":10})
        r1_new = Reaction([self.s1, s_new], [c1.replace_species(s_old, s_new)], propensity_type = "proportionalhillpositive", propensity_params = {"k":1., "d":s_new, "s1":self.s1, "n":2, "K":10})
        self.assertTrue(r1.replace_species(s_old, s_new) == r1_new)

        #test replace in a chemical reaction network
        c1 = ComplexSpecies([self.s1, s_old])
        c2 = ComplexSpecies([self.s1, c1])
        r1 = Reaction([self.s1, s_old], [c1], k = 1)
        species = [self.s1, s_old, c1, c2]
        r1 = Reaction([self.s1, s_old], [c1], k = 1)
        crn = ChemicalReactionNetwork(species = species, reactions = [r1])
        new_crn = crn.replace_species(s_old, s_new)

        self.assertTrue(self.s1 in new_crn.species)
        self.assertFalse(s_old in new_crn.species)

        self.assertTrue(s_new in new_crn.species)
        self.assertFalse(c1 in new_crn.species)
        self.assertFalse(c2 in new_crn.species)
        c1_new = ComplexSpecies([self.s1, s_new])
        c2_new = ComplexSpecies([self.s1, c1_new])
        self.assertTrue(c1_new in new_crn.species)
        self.assertTrue(c2_new in new_crn.species)
        r1_new = Reaction([self.s1, s_new], [c1_new], k = 1)
        self.assertFalse(r1 in new_crn.reactions)
        self.assertTrue(r1_new in new_crn.reactions)


    def test_generate_sbml_model(self):

        # generate an sbml model
        document, model = self.crn.generate_sbml_model()
        # all species from the CRN are accounted for
        self.assertEqual(len(model.getListOfSpecies()), len(self.crn.species))
        # all reactions from the CRN are accounted for
        self.assertEqual(len(model.getListOfReactions()), len(self.crn.reactions))

        # test a reversible reaction
        rx1 = Reaction.from_mass_action(inputs=[self.s1], outputs=[self.s2], k_forward=0.1, k_reverse=0.1)
        rxn_list = [rx1]
        crn = ChemicalReactionNetwork(species=self.species_list, reactions=rxn_list)

        # generate an sbml model
        document, model = crn.generate_sbml_model()
        # all species from the CRN are accounted for
        self.assertEqual(len(model.getListOfSpecies()), len(crn.species))
        # all reactions from the CRN are accounted for
        # the sbml represents a reverisble reaction with to separate reactions
        self.assertEqual(len(model.getListOfReactions()), 2*len(crn.reactions))

    def test_write_sbml_file(self):

        document, _ = self.crn.generate_sbml_model(model_id = 'test_model')
        sbml_string = libsbml.writeSBMLToString(document)

        file_name = 'test_sbml.xml'
        with patch("builtins.open", new=mock_open()) as _file:
            self.crn.write_sbml_file(file_name, model_id = 'test_model')

            _file.assert_called_once_with(file_name, 'w')
            _file().write.assert_called_once_with(sbml_string)
