
#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase


class TestComponent(TestCase):
    def test_initial_concentration(self):
        from biocrnpyler import Component

        comp_name = 'test_component'

        parameters = {"mock_param": 1}
        component = Component(name=comp_name, mechanisms={}, parameters=parameters, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        self.assertTrue(component.initial_concentration == 0)
        component.initial_concentration = 5
        self.assertTrue(component.initial_concentration == 5)

        with self.assertRaises(ValueError):
            component.initial_concentration = -1

    def test_get_species(self):
        from biocrnpyler import Component

        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        with self.assertWarns(Warning):
            component.get_species()

    def test_set_attributes(self):
        from biocrnpyler import Component

        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        attr_list = ['attr1','attr2']
        component.set_attributes(attr_list)
        self.assertTrue(len(component.attributes) == len(attr_list))

        with self.assertRaises(RuntimeError):
            attr_list = [5, 2]
            component.set_attributes(attr_list)

    def test_add_attribute(self):

        from biocrnpyler import Component
        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        with self.assertRaises(AssertionError):
            component.add_attribute(2)

        attr = 'atr1'
        with self.assertRaises(Warning):
            component.add_attribute(attr)

    def test_update_parameters(self):
        from biocrnpyler import Component
        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
        parameters = {"kb": kb, "ku": ku, "ktx": ktx, "ktl": ktl, "kdeg": kdeg}

        self.assertTrue(isinstance(component.custom_parameters,dict) and len(component.custom_parameters) == 0)

        one_param = {"kb": kb}
        component.custom_parameters = one_param

        component.update_parameters(mixture_parameters={}, parameters=parameters, overwrite_custom_parameters=False)
        self.assertTrue(len(component.custom_parameters) == len(one_param))

        component.update_parameters(mixture_parameters={}, parameters=parameters)
        self.assertEqual(component.custom_parameters, one_param)

        component.update_parameters(mixture_parameters=parameters, parameters={})
        self.assertEqual(component.parameters, parameters)

    def test_update_mechanisms(self):
        from biocrnpyler import Component
        from biocrnpyler import Transcription_MM
        from biocrnpyler import Translation_MM
        from biocrnpyler import Degredation_mRNA_MM

        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        tx = Transcription_MM()
        tl = Translation_MM()
        deg = Degredation_mRNA_MM()

        test_mech = {tx.mechanism_type: tx, tl.mechanism_type: tl}

        default_test_mech = {deg.mechanism_type : deg}

        self.assertTrue(isinstance(component.mechanisms, dict) and len(component.mechanisms) == 0)
        component.update_mechanisms(mixture_mechanisms=test_mech)
        self.assertEqual(component.mechanisms, test_mech)

        component.default_mechanisms = default_test_mech
        component.update_mechanisms(mixture_mechanisms=test_mech)
        test_mech.update(default_test_mech)
        self.assertEqual(component.mechanisms, test_mech)

        self.assertTrue(isinstance(component.custom_mechanisms, dict) and len(component.custom_mechanisms) == 0)
        component.update_mechanisms(mechanisms=test_mech)
        self.assertEqual(component.custom_mechanisms, test_mech)

        component.update_mechanisms(mechanisms=default_test_mech,overwrite_custom_mechanisms=False)
        self.assertEqual(component.custom_mechanisms, test_mech)

        test_mech_list = list(test_mech.values())
        component.update_mechanisms(mechanisms=test_mech_list)
        self.assertEqual(component.custom_mechanisms, test_mech)

        with self.assertRaises(ValueError):
            component.update_mechanisms(mechanisms=(tx, tl))

    def test_get_parameter(self):
        from biocrnpyler import Component
        from biocrnpyler import Transcription_MM

        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        with self.assertRaises(ValueError):
            component.get_parameter(param_name='kb')

        # Create Param Dict
        kb, ku, ktx, ktl, kdeg, cooperativity = 100, 10, 3, 2, 1, 1
        p_id = 'p10'

        parameters = {"kb": kb, "ku": ku, "ktx": ktx, "ktl": ktl, "kdeg": kdeg, "cooperativity": cooperativity,
                      # default params
                      ("transcription", "ktx"): ktx,
                      ("transcription", p_id, 'ku'): ku
                      }

        one_param = {"kb": kb}
        component.update_parameters(mixture_parameters={}, parameters=one_param)

        self.assertEqual(component.get_parameter(param_name="kb"), one_param["kb"])

        component.update_parameters(parameters=parameters)

        tx = Transcription_MM()
        self.assertEqual(component.get_parameter(mechanism=tx, param_name='ktx'), ktx)

        self.assertEqual(component.get_parameter(mechanism=tx, part_id=p_id, param_name='ku'), ku)

    def test_update_species(self):
        from biocrnpyler import Component
        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        with self.assertWarns(Warning):
            component.update_species()

    def test_update_reactions(self):
        from biocrnpyler import Component
        comp_name = 'test_component'

        component = Component(name=comp_name, mechanisms={}, parameters={}, parameter_file=None,
                              mixture=None, attributes=[], initial_conc=0, parameter_warnings=True)

        with self.assertWarns(Warning):
            component.update_reactions()
