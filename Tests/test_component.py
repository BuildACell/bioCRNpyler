
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Component, DNA, ParameterDatabase, Mixture, Mechanism, Species
from biocrnpyler import SimpleTranscription, SimpleTranslation


class TestComponent(TestCase):

    def setUp(self) -> None:
        """This method gets called before every test"""
        self.comp_name = 'test_component'
        self.default_concentration = 0
        self.component = Component(name=self.comp_name, mechanisms={}, parameters={}, parameter_file=None,
                                   mixture=None, attributes=[], initial_concentration=self.default_concentration)

    def test_initial_concentration(self):

        # test that the default initial concentration is zero
        self.assertEqual(self.component.initial_concentration, self.default_concentration)
        #test there is one entry in the parameter database (the initial concentration)
        self.assertTrue(len(self.component.parameter_database.parameters) == 1)
        param = self.component.parameter_database.find_parameter(mechanism = 'initial concentration', part_id = None, param_name = self.comp_name)
        self.assertTrue(param.value == self.default_concentration)

        new_value = 5
        self.component.initial_concentration = new_value
        # test that the initial concentration has been modified
        self.assertEqual(self.component.initial_concentration, new_value)

        #test the value in the param dictionary has changes
        param = self.component.parameter_database.find_parameter(mechanism = 'initial concentration', part_id = None, param_name = self.comp_name)
        self.assertTrue(param.value == new_value)

        not_valid_value = -1
        with self.assertRaisesRegex(ValueError, f'Initial concentration must be non-negative, this was given: {not_valid_value}'):
            self.component.initial_concentration = not_valid_value

        
        

    def test_get_species(self):

        # components has no valid get_species function, it returns None
        self.assertTrue(self.component.get_species() is None)

    def test_set_attributes(self):

        attr_list = ['attr1', 'attr2']
        
        with self.assertRaisesRegex(Warning,f'Component {self.component.name} has no internal species and therefore no attributes'):
            self.component.set_attributes(attr_list)

        # DNA is inherited from component and has valid internal species
        dna = DNA("dna")
        # it calls compoents.set_attributes 
        dna.set_attributes(attr_list)
        # test that the correct attributes are set
        self.assertTrue(len(dna.attributes) == len(attr_list))

    def test_update_parameters(self):
        kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
        parameters = {"kb": kb, "ku": ku, "ktx": ktx, "ktl": ktl, "kdeg": kdeg}

        # test that the custom parameters dictionary is empty
        component = Component(self.comp_name)
        self.assertTrue(isinstance(component.parameter_database, ParameterDatabase)
                        and len(component.parameter_database.parameters) == 0)

        

        component.update_parameters(parameters=parameters)

        # test that the component has all the parameters
        self.assertTrue(len(component.parameter_database.parameters) == len(parameters))

        # test overwriting parameters
        new_val = 111
        one_param = {"kb": new_val}
        component.update_parameters(parameters=one_param, overwrite_parameters = True)
        self.assertEqual(component.parameter_database[(None, None, "kb")].value, new_val)
        # test that the parameter dictionary is still the same length as before
        self.assertTrue(len(component.parameter_database.parameters) == len(parameters))


    def test_add_mechanism(self):
        tx = SimpleTranscription()
        tl = SimpleTranslation()

        #test adding a single mechanism instead of a list still works
        self.component.add_mechanisms(tx)
        self.assertTrue(tx.mechanism_type in self.component.mechanisms)

        #add a non-mechanism
        with self.assertRaisesRegex(TypeError, 'mechanism must be a Mechanism.'):
            self.component.add_mechanism(None)

        with self.assertRaisesRegex(ValueError, 'add_mechanisms expected a list of Mechanisms.'):
            self.component.add_mechanisms(None)

        #add same mechanism, new type
        self.component.add_mechanism(tx, mech_type = "new")
        self.assertTrue("new" in self.component.mechanisms)
        self.assertTrue(type(self.component.mechanisms["new"]) == SimpleTranscription)

        #Add invalid mech_type
        with self.assertRaisesRegex(TypeError, 'mechanism keys must be strings.'):
            self.component.add_mechanism(tx, mech_type = 1234)

        #add a mechanism already in the Component
        with self.assertRaisesRegex(ValueError, f"mech_type {tx.mechanism_type} already in component"):
            self.component.add_mechanism(tx)

        #add mechanism with optional_mechanism - does not overwrite!
        self.component.add_mechanism(tl, mech_type = tx.mechanism_type, optional_mechanism = True)
        self.assertTrue(tx.mechanism_type in self.component.mechanisms)
        self.assertTrue(type(self.component.mechanisms[tx.mechanism_type]) == SimpleTranscription)

        #add mechanism with overwrite
        self.component.add_mechanism(tl, mech_type = tx.mechanism_type, overwrite = True)
        self.assertTrue(tx.mechanism_type in self.component.mechanisms)
        self.assertTrue(type(self.component.mechanisms[tx.mechanism_type]) == SimpleTranslation)

    def test_update_mechanisms(self):

        tx = SimpleTranscription()
        tl = SimpleTranslation()

        test_mech = {tx.mechanism_type: tx, tl.mechanism_type: tl}

        # test that component has no mechanism
        self.assertTrue(isinstance(self.component.mechanisms, dict) and len(self.component.mechanisms) == 0)

        #test mechanism setter
        self.component.mechanisms = test_mech
        self.assertEqual(self.component.mechanisms.keys(), test_mech.keys())

        #test mechanisms are copied
        self.assertEqual(type(self.component.mechanisms[tx.mechanism_type]),  type(test_mech[tx.mechanism_type]))
        self.assertFalse(self.component.mechanisms[tx.mechanism_type] == test_mech[tx.mechanism_type])

        #remove all mechanisms
        self.component.mechanisms = {}
        self.assertEqual(self.component.mechanisms, {})

        #test add_mechanism
        self.component.add_mechanism(tx, tx.mechanism_type)
        self.component.add_mechanism(tl)
        self.assertEqual(self.component.mechanisms.keys(), test_mech.keys())

        #test add_mechanisms with list
        self.component.mechanisms = {}
        test_mech_list = list(test_mech.values())
        self.component.add_mechanisms(test_mech_list)
        self.assertEqual(self.component.mechanisms.keys(), test_mech.keys())

    def test_get_parameter(self):

        # testing an invalid parameter
        with self.assertRaisesRegex(ValueError, 'No parameters can be found that match the'):
            self.component.get_parameter(param_name='kb')

        # Create Param Dict
        kb, ku, ktx, ktl, kdeg, cooperativity = 100, 10, 3, 2, 1, 4
        p_id = 'p10'

        parameters = {"kb": kb, "kdeg":kdeg,
                      ("transcription", None, "ktx"): ktx,
                      ("transcription", p_id, 'ku'): ku,
                      (None, p_id, "ktl"): ktl }
        # update the component parameters
        self.component.update_parameters(parameters=parameters)

        tx = SimpleTranscription()
        # testing the different parameter definitions
        self.assertEqual(self.component.get_parameter(param_name='kb').value, kb)
        self.assertEqual(self.component.get_parameter(mechanism=tx, param_name='ktx').value, ktx)
        self.assertEqual(self.component.get_parameter(part_id=p_id, param_name='ktl').value, ktl)
        self.assertEqual(self.component.get_parameter(mechanism=tx, part_id=p_id, param_name='ku').value, ku)

        # testing parameter with return_numerical = True
        self.assertEqual(self.component.get_parameter(param_name='kb', return_numerical = True), kb)

        one_param = {"kb": kb}
        self.component.update_parameters(parameters=one_param)
        # testing that one_param was registered
        self.assertEqual(self.component.get_parameter(param_name="kb").value, one_param["kb"])

        #testing a parameter which can't be found
        with self.assertRaisesRegex(ValueError, "No parameters can be found"):
            self.component.get_parameter("not a parameter")

    def test_update_species(self):
        # warning if update_species on a component object
        with self.assertWarnsRegex(Warning, f'Unsubclassed update_species called for {self.component}'):
            self.component.update_species()

    def test_update_reactions(self):
        # warning if update_reaction on a component object
        with self.assertWarnsRegex(Warning, f'Unsubclassed update_reactions called for {self.component}'):
            self.component.update_reactions()

    def test_get_mechanism(self):
        M1_comp = Mechanism(name = "m1_comp", mechanism_type = "shared")
        M1_mix = Mechanism(name = "m1_mix", mechanism_type = "shared")
        M2_comp = Mechanism(name = "m2_comp", mechanism_type = "comp")
        M2_mix = Mechanism(name = "m2_mix", mechanism_type = "mixture")

        #Create a Mixture and Component with the above mechanisms
        C = Component(name = "comp", mechanisms = [M1_comp, M2_comp])
        M = Mixture(mechanisms = [M1_mix, M2_mix], components = [C])

        #Get the copy of C in M
        C_copy = M.get_component(component = C)

        self.assertTrue(C_copy.get_mechanism("shared").name == "m1_comp")
        self.assertTrue(M.get_mechanism("shared").name == "m1_mix")
        self.assertTrue(C_copy.get_mechanism("comp").name == "m2_comp")
        self.assertTrue(C_copy.get_mechanism("mixture").name == "m2_mix")

        #Make sure the Mixture get_mechanism works as well, just in case.
        self.assertTrue(M.get_mechanism("comp") is None)

        #test get Mechanism with no_key_error = False (Default)
        with self.assertRaisesRegex(KeyError, "Unable to find mechanism of type"):
            C_copy.get_mechanism("DNE")

        #test get_mechanism with no_key_error = True
        self.assertTrue(C_copy.get_mechanism("DNE", optional_mechanism = True) is None)

    def test_set_species(self):
        C = Component(name = "comp")

        dna_string = "dna_S"
        dna_species = Species(name = "S", material_type = "dna")
        dna_comp = DNA(name = "S")

        #add_species can take a string, Species, or Component
        s1 = C.set_species(dna_species)
        s2 = C.set_species(dna_string)
        s3 = C.set_species(dna_comp)
        self.assertEqual(str(s1), str(s2)) #A Species derived from a string will have no material_type
        self.assertFalse(s1 == s2)
        self.assertEqual(s1, s3)
        self.assertEqual(s1, dna_species)

        #add_species can also take a list
        s_list = C.set_species([dna_string, dna_species, dna_comp])
        self.assertTrue(len(s_list) == 3)
        self.assertTrue(str(s_list[0]) == str(s_list[1]) == str(s_list[2]))

        #The following cases should raise errors
        with self.assertRaisesRegex(ValueError, ""):
            C.set_species(None)

        with self.assertRaisesRegex(ValueError, ""):
            C.set_species([dna_species, None])

