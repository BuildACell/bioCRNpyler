
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Component, DNA, ParameterDatabase, Mixture, Mechanism
from biocrnpyler import SimpleTranscription, SimpleTranslation


class TestComponent(TestCase):

    def setUp(self) -> None:
        """This method gets called before every test"""
        self.comp_name = 'test_component'
        self.default_concentration = 0
        self.component = Component(name=self.comp_name, mechanisms={}, parameters={}, parameter_file=None,
                                   mixture=None, attributes=[], initial_conc=self.default_concentration)

    def test_initial_concentration(self):

        # test that the default initial concentration is zero
        self.assertEqual(self.component.initial_concentration, self.default_concentration)
        new_value = 5
        self.component.initial_concentration = new_value
        # test that the initial concentration has been modified
        self.assertEqual(self.component.initial_concentration, new_value)

        not_valid_value = -1
        with self.assertRaisesRegexp(ValueError, f'Initial concentration must be non-negative, this was given: {not_valid_value}'):
            self.component.initial_concentration = not_valid_value

    def test_get_species(self):

        # components has no valid get_species function, it returns None
        self.assertTrue(self.component.get_species() is None)

    def test_set_attributes(self):

        attr_list = ['attr1', 'attr2']
        
        with self.assertRaisesRegexp(Warning,f'Component {self.component.name} has no internal species and therefore no attributes'):
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
        self.assertTrue(isinstance(self.component.parameter_database, ParameterDatabase)
                        and len(self.component.parameter_database.parameters) == 0)

        

        self.component.update_parameters(parameters=parameters)

        # test that the component has all the parameters
        self.assertTrue(len(self.component.parameter_database.parameters) == len(parameters))

        # test overwriting parameters
        new_val = 111
        one_param = {"kb": new_val}
        self.component.update_parameters(parameters=one_param, overwrite_parameters = True)
        self.assertEqual(self.component.parameter_database[(None, None, "kb")].value, new_val)
        # test that the parameter dictionary is still the same length as before
        self.assertTrue(len(self.component.parameter_database.parameters) == len(parameters))

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
        with self.assertRaisesRegexp(ValueError, 'No parameters can be found that match the'):
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
        with self.assertRaisesRegexp(ValueError, "No parameters can be found"):
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
        with self.assertRaisesRegexp(KeyError, "Unable to find mechanism of type"):
            C_copy.get_mechanism("DNE")

        #test get_mechanism with no_key_error = True
        self.assertTrue(C_copy.get_mechanism("DNE", optional_mechanism = True) is None)


