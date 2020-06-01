
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Component, DNA
from biocrnpyler import Transcription_MM, Translation_MM, Degredation_mRNA_MM


class TestComponent(TestCase):

    def setUp(self) -> None:
        """This method gets called before every test"""
        self.comp_name = 'test_component'
        self.default_concentration = 0
        self.component = Component(name=self.comp_name, mechanisms={}, parameters={}, parameter_file=None,
                                   mixture=None, attributes=[], initial_conc=self.default_concentration, parameter_warnings=True)

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

        # components has no valid get_species function, it raise warning when it gets called
        with self.assertWarnsRegex(Warning, f'get_species is not defined for component {self.component.name}, None returned.'):
            self.component.get_species()

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
        self.assertTrue(isinstance(self.component.custom_parameters, dict)
                        and len(self.component.custom_parameters) == 0)

        one_param = {"kb": kb}
        self.component.custom_parameters = one_param

        self.component.update_parameters(mixture_parameters={}, parameters=parameters, overwrite_custom_parameters=False)
        # test that the new parameter is the only parameter custom parameter in the component
        self.assertTrue(len(self.component.custom_parameters) == len(one_param))

        self.component.update_parameters(mixture_parameters={}, parameters=parameters)
        self.assertEqual(self.component.custom_parameters, one_param)

        self.component.update_parameters(mixture_parameters=parameters, parameters={})
        # test that the parameter dictionary is still the same as before
        self.assertEqual(self.component.parameters, parameters)

    def test_update_mechanisms(self):

        tx = Transcription_MM()
        tl = Translation_MM()
        deg = Degredation_mRNA_MM()

        test_mech = {tx.mechanism_type: tx, tl.mechanism_type: tl}

        default_test_mech = {deg.mechanism_type: deg}

        # test that component has no mechanism
        self.assertTrue(isinstance(self.component.mechanisms, dict) and len(self.component.mechanisms) == 0)
        self.component.update_mechanisms(mixture_mechanisms=test_mech)
        # test that the test_mech is registered as the only mechanism
        self.assertEqual(self.component.mechanisms, test_mech)

        self.component.default_mechanisms = default_test_mech
        self.component.update_mechanisms(mixture_mechanisms=test_mech)
        test_mech.update(default_test_mech)
        self.assertEqual(self.component.mechanisms, test_mech)

        # testing that the custom mechanism gets updated
        self.assertTrue(isinstance(self.component.custom_mechanisms, dict) and len(self.component.custom_mechanisms) == 0)
        self.component.update_mechanisms(mechanisms=test_mech)
        self.assertEqual(self.component.custom_mechanisms, test_mech)

        # testing that custom mechanism is protected by the overwrite_custom_mechanisms=False flag
        self.component.update_mechanisms(mechanisms=default_test_mech, overwrite_custom_mechanisms=False)
        self.assertEqual(self.component.custom_mechanisms, test_mech)

        # multiple mechanisms can be supplied with a list
        test_mech_list = list(test_mech.values())
        self.component.update_mechanisms(mechanisms=test_mech_list)
        self.assertEqual(self.component.custom_mechanisms, test_mech)

        # testing an invalid mechanism format
        with self.assertRaisesRegexp(ValueError, 'Mechanisms must be passed as a list of instantiated objects or a '
                                                 'dictionary {type:mechanism}'):
            self.component.update_mechanisms(mechanisms=(tx, tl))

    def test_get_parameter(self):

        # testing an invalid parameter
        with self.assertRaisesRegexp(ValueError, 'No parameters can be found that match the'):
            self.component.get_parameter(param_name='kb')

        # Create Param Dict
        kb, ku, ktx, ktl, kdeg, cooperativity = 100, 10, 3, 2, 1, 1
        p_id = 'p10'

        parameters = {"kb": kb, "ku": ku, "ktx": ktx, "ktl": ktl, "kdeg": kdeg, "cooperativity": cooperativity,
                      # default params
                      ("transcription", "ktx"): ktx,
                      ("transcription", p_id, 'ku'): ku
                      }

        one_param = {"kb": kb}
        self.component.update_parameters(mixture_parameters={}, parameters=one_param)
        # testing that one_param was registered
        self.assertEqual(self.component.get_parameter(param_name="kb"), one_param["kb"])
        # update the component parameters
        self.component.update_parameters(parameters=parameters)
        # testing the different parameter definitions
        tx = Transcription_MM()
        self.assertEqual(self.component.get_parameter(mechanism=tx, param_name='ktx'), ktx)

        self.assertEqual(self.component.get_parameter(mechanism=tx, part_id=p_id, param_name='ku'), ku)

    def test_update_species(self):
        # warning if update_species on a component object
        with self.assertWarnsRegex(Warning, f'Unsubclassed update_species called for {self.component}'):
            self.component.update_species()

    def test_update_reactions(self):
        # warning if update_reaction on a component object
        with self.assertWarnsRegex(Warning, f'Unsubclassed update_reactions called for {self.component}'):
            self.component.update_reactions()
