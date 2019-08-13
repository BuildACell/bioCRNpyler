
#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
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

        rx1 = Reaction(inputs=[s1], outputs=[s2], k=0.1)
        rxn_list = [rx1]

        checked_species, checked_reactions = ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list,
                                                                                        species=species_list)

        self.assertEqual(species_list, checked_species)

        self.assertEqual(rxn_list, checked_reactions)

        species_list_with_none = species_list.copy()
        species_list_with_none.append(None)
        with self.assertRaises(ValueError):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list, species=species_list_with_none)

        rxn_list_with_none = rxn_list.copy()
        rxn_list_with_none.append(None)
        with self.assertRaises(ValueError):
            ChemicalReactionNetwork.check_crn_validity(reactions=rxn_list_with_none, species=species_list)

        s3 = Species(name='test_species3')
        s4 = Species(name='test_species4')

        rxn2 = Reaction(inputs=[s1], outputs=[s3], k=0.1)
        with self.assertWarns(Warning):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn2], species=species_list, warnings=True)

        rxn3 = Reaction(inputs=[s4], outputs=[s2], k=0.1)
        with self.assertWarns(Warning):
            ChemicalReactionNetwork.check_crn_validity(reactions=[rxn3], species=species_list, warnings=True)

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

        init_cond = {s1: 5, s3: 10}

        x0 = crn.initial_condition_vector(init_cond_dict=init_cond)

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

        s3 = Species(name='test_species3')

        with self.assertRaises(ValueError):
            crn.get_all_species_containing(species=species_list)

        rtn_species_list = crn.get_all_species_containing(species=s3)
        self.assertEqual(rtn_species_list, [])

        rtn_species_list = crn.get_all_species_containing(species=s1)
        self.assertEqual(rtn_species_list, [s1])

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

        document, model = crn.generate_sbml_model()

        self.assertEqual(len(model.getListOfSpecies()), len(crn.species))
        self.assertEqual(len(model.getListOfReactions()), len(crn.reactions))

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

