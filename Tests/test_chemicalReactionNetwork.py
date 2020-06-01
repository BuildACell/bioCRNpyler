
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from unittest.mock import mock_open, patch


class TestChemicalReactionNetwork(TestCase):
    def test_check_crn_validity(self):
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction

        s1 = Species(name='test_species1')
        s2 = Species(name='test_species2')
        species_list = [s1, s2]
        # creating a valid reaction two species
        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        checked_species, checked_reactions = ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list,
                                                                                        species=species_list)
        # test that the returned species list is the same as the species list supplied
        self.assertEqual(species_list, checked_species)
        # test that the returned reaction list is the same as the reaction list supplied
        self.assertEqual(rxn_list, checked_reactions)

        species_list_with_none = species_list.copy()
        # injecting a None to the species list
        species_list_with_none.append(None)
        # test whether a non-species object is detected and Value error has been raised
        with self.assertRaisesRegexp(ValueError, 'A non-species object was used as a species!'):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list, species=species_list_with_none)

        rxn_list_with_none = rxn_list.copy()
        # injecting a None to the reaction list
        rxn_list_with_none.append(None)
        # test whether a non-reaction object is detected and Value Error has been raised
        with self.assertRaisesRegexp(ValueError, 'A non-reaction object was used as a reaction!'):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list_with_none, species=species_list)

        s3 = Species(name='test_species3')
        s4 = Species(name='test_species4')

        rxn2 = Reaction(inputs=[s1], outputs=[s3], k=0.1)
        # test warning raised if a species (in the reaction outputs) is detected which is not part of the species list
        with self.assertWarnsRegex(Warning, f'contains a species {s3.name} which is not in the CRN'):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn2], species=species_list, warnings=True)

        rxn3 = Reaction(inputs=[s4], outputs=[s2], k=0.1)
        # test warning raised if a species (in the reaction inputs) is detected which is not part of the species list
        with self.assertWarnsRegex(Warning, f'contains a species {s4.name} which is not in the CRN'):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn3], species=species_list, warnings=True)

        # test duplicate reactions
        rxn_list = [rx1,rx1]
        with self.assertWarnsRegex(Warning, 'may be duplicated in CRN definitions'):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list, species=species_list, warnings=True)

    def test_species_index(self):
        # TODO add test if we actually use this function
        pass

    def test_initial_condition_vector(self):
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction

        s1 = Species(name='test_species1')
        s2 = Species(name='test_species2')

        species_list = [s1, s2]

        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        crn = ChemicalReactionNetwork(species=species_list, reactions=rxn_list)

        s3 = Species(name='test_species3')
        # no initial value is supplied for s2
        init_cond = {s1: 5, s3: 10}

        x0 = crn.initial_condition_vector(init_cond_dict=init_cond)

        # test that value for s3 is not in the initial value list because s3 is not part of any reaction
        not_in_the_list = False
        for key, value in init_cond.items():
            if value not in x0 and key == s3:
                not_in_the_list = True
        self.assertTrue(not_in_the_list)

    def test_get_all_species_containing(self):
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction

        s1 = Species(name='test_species1')
        s2 = Species(name='test_species2')

        species_list = [s1, s2]

        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        crn = ChemicalReactionNetwork(species=species_list, reactions=rxn_list)

        # test that the species arument must be Species object
        with self.assertRaisesRegexp(ValueError, 'species argument must be an instance of Species!'):
            crn.get_all_species_containing(species=species_list)

        # s3 is not part of the CRN
        s3 = Species(name='test_species3')
        rtn_species_list = crn.get_all_species_containing(species=s3)
        # test that empty list is returned
        self.assertEqual(rtn_species_list, [])

        # test that s1 is returned only once
        rtn_species_list = crn.get_all_species_containing(species=s1)
        self.assertEqual(rtn_species_list, [s1])

        # test that s1 is returned only once as a string
        rtn_species_list = crn.get_all_species_containing(species=s1, return_as_strings=True)
        self.assertEqual(rtn_species_list, [repr(s1)])

    def test_generate_sbml_model(self):
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction

        s1 = Species(name='test_species1')
        s2 = Species(name='test_species2')

        species_list = [s1, s2]

        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        crn = ChemicalReactionNetwork(species=species_list, reactions=rxn_list)
        # generate an sbml model
        document, model = crn.generate_sbml_model()
        # all species from the CRN are accounted for
        self.assertEqual(len(model.getListOfSpecies()), len(crn.species))
        # all reactions from the CRN are accounted for
        self.assertEqual(len(model.getListOfReactions()), len(crn.reactions))

        # test a reversible reaction
        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1, k_rev=0.1)
        rxn_list = [rx1]
        crn = ChemicalReactionNetwork(species=species_list, reactions=rxn_list)

        # generate an sbml model
        document, model = crn.generate_sbml_model()
        # all species from the CRN are accounted for
        self.assertEqual(len(model.getListOfSpecies()), len(crn.species))
        # all reactions from the CRN are accounted for
        # the sbml represents a reverisble reaction with to separate reactions
        self.assertEqual(len(model.getListOfReactions()), 2*len(crn.reactions))


    def test_write_sbml_file(self):

        import libsbml
        from biocrnpyler import ChemicalReactionNetwork
        from biocrnpyler import Species
        from biocrnpyler import Reaction

        s1 = Species(name='test_species1')
        s2 = Species(name='test_species2')

        species_list = [s1, s2]

        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        crn = ChemicalReactionNetwork(species=species_list, reactions=rxn_list)
        document, _ = crn.generate_sbml_model()
        sbml_string = libsbml.writeSBMLToString(document)

        file_name = 'test_sbml.xml'
        with patch("builtins.open", new=mock_open()) as _file:
            crn.write_sbml_file(file_name)

            _file.assert_called_once_with(file_name, 'w')
            _file().write.assert_called_once_with(sbml_string)

