#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import warnings
from unittest import TestCase

from biocrnpyler import (
    DNA,
    Component,
    DNAassembly,
    Mixture,
    Module,
    SimpleTranscription,
    SimpleTranslation,
    Species,
)


def txtl_mechanisms():
    """Mechanisms needed to compile a DNAassembly."""
    return {
        'transcription': SimpleTranscription(),
        'translation': SimpleTranslation(),
    }


def transcription_rate(crn, name):
    """Rate constant of the transcription reaction producing rna_<name>."""
    for reaction in crn.reactions:
        inputs = [str(w.species) for w in reaction.inputs]
        outputs = [str(w.species) for w in reaction.outputs]
        if f"rna_{name}" in outputs and f"rna_{name}" not in inputs:
            return reaction.propensity_type.k_forward
    raise AssertionError(f"no transcription reaction found for {name}")


class TestModule(TestCase):
    def test_add_components(self):
        module = Module('test_module')
        # a new module has no components
        self.assertTrue(len(module.components) == 0)

        component = Component('test_comp')
        module.add_component(component)
        self.assertTrue(len(module.components) == 1)

        # test that it was added by copying
        self.assertTrue(component not in module.components)
        self.assertTrue(module.components[0].name == component.name)

        # test that the same component cannot be added again
        with self.assertRaisesRegex(
            ValueError, 'of same type and name already in Module!'
        ):
            module.add_component(component)

        # components of the same name but different type are distinct
        module.add_component(DNA('test_comp'))
        self.assertTrue(len(module.components) == 2)

        # species are invalid components
        with self.assertRaisesRegex(
            ValueError, 'add_components expected a list of Components.'
        ):
            module.add_components(Species('test_species'))

        # use the constructor the other way
        module = Module('test_module', components=[component])
        self.assertTrue(len(module.components) == 1)
        self.assertTrue(component not in module.components)

    def test_add_species(self):
        module = Module('test_module')
        self.assertTrue(len(module.species) == 0)

        species = Species('test_species')
        module.add_species(species)
        self.assertEqual([species], module.species)

        # adding invalid species
        with self.assertRaisesRegex(
            ValueError, 'add_species expected a list of Species.'
        ):
            module.add_species(['not_a_species'])

        # use the constructor the other way
        module = Module('test_module', species=[species])
        self.assertEqual([species], module.species)

    def test_update_species_and_reactions(self):
        species = Species('test_species')
        module = Module(
            'test_module', components=[Component('c1')], species=[species]
        )

        # the module contributes its added species and no reactions
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter('always')
            self.assertEqual([species], module.update_species())
            self.assertEqual([], module.update_reactions())
        # unlike a bare Component, a Module does not warn about being
        # unsubclassed
        self.assertEqual(0, len(recorded))

    def test_apply_context_mechanisms(self):
        own_mechanism = SimpleTranscription()
        component = Component(
            'c1', mechanisms={'transcription': own_mechanism}
        )
        module = Module('test_module', mechanisms=txtl_mechanisms())

        contextualized = module.apply_context(component)

        # the module fills in a mechanism the component does not define
        self.assertTrue('translation' in contextualized.mechanisms)
        # the component's own mechanism is not overwritten
        self.assertEqual(
            own_mechanism.name,
            contextualized.mechanisms['transcription'].name,
        )
        # the original component is left alone
        self.assertTrue('translation' not in component.mechanisms)

    def test_apply_context_parameters(self):
        component = Component('c1', parameters={'ktx': 3.0})
        module = Module('test_module', parameters={'ktx': 2.0, 'ktl': 0.01})

        contextualized = module.apply_context(component)
        database = contextualized.parameter_database

        # the module fills in a parameter the component does not define
        self.assertEqual(
            0.01, database.find_parameter(None, None, 'ktl').value
        )
        # the component's own parameter is not overwritten
        self.assertEqual(
            3.0, database.find_parameter(None, None, 'ktx').value
        )
        # the original component is left alone
        self.assertTrue(
            component.parameter_database.find_parameter(None, None, 'ktl')
            is None
        )

    def test_apply_context_reaches_sub_components(self):
        # A DNAassembly builds its promoter and RBS when it is constructed,
        # so the module's context has to reach them after the fact.
        assembly = DNAassembly('A', promoter='P', rbs='B')
        module = Module(
            'test_module',
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0, 'ktl': 0.01},
        )

        contextualized = module.apply_context(assembly)

        for part in [contextualized.promoter, contextualized.rbs]:
            self.assertTrue('transcription' in part.mechanisms)
            self.assertEqual(
                2.0,
                part.parameter_database.find_parameter(
                    None, None, 'ktx'
                ).value,
            )

    def test_apply_context_shares_parameter_entries(self):
        # Entries are shared rather than rebuilt so that units and
        # parameter origins survive.
        module = Module('test_module', parameters={'ktx': 2.0})
        key = list(module.parameter_database.parameters)[0]
        entry = module.parameter_database.parameters[key]

        contextualized = module.apply_context(Component('c1'))
        self.assertIs(
            entry, contextualized.parameter_database.parameters[key]
        )

    def test_enumerate_components(self):
        component = Component('c1')
        module = Module(
            'test_module',
            components=[component, Component('c2')],
            mechanisms=txtl_mechanisms(),
        )

        enumerated = module.enumerate_components()

        self.assertEqual(2, len(enumerated))
        # enumeration returns copies carrying the module's context
        self.assertTrue(all(c not in module.components for c in enumerated))
        self.assertTrue(
            all('transcription' in c.mechanisms for c in enumerated)
        )
        # the module's stored components are unchanged
        self.assertTrue(
            all(
                'transcription' not in c.mechanisms for c in module.components
            )
        )

    def test_module_addition(self):
        m1 = Module(
            'm1',
            components=[Component('c1')],
            mechanisms={'transcription': SimpleTranscription()},
            species=[Species('s1')],
        )
        m2 = Module(
            'm2',
            components=[Component('c2')],
            mechanisms={'translation': SimpleTranslation()},
            species=[Species('s2')],
        )

        combined = m1 + m2

        self.assertTrue(isinstance(combined, Module))
        self.assertEqual('m1_m2', combined.name)
        self.assertEqual(
            ['c1', 'c2'], sorted(c.name for c in combined.components)
        )
        self.assertEqual(
            ['s1', 's2'], sorted(s.name for s in combined.species)
        )

        # each module's mechanisms are baked into its own components, so
        # the two modules keep their separate mechanisms after merging
        by_name = {c.name: c for c in combined.components}
        self.assertEqual(['transcription'], list(by_name['c1'].mechanisms))
        self.assertEqual(['translation'], list(by_name['c2'].mechanisms))

        # the originals are unchanged
        self.assertEqual(1, len(m1.components))
        self.assertEqual(1, len(m2.components))

    def test_module_addition_duplicate_components(self):
        # components shared by two modules are how modules wire together,
        # so a duplicate is a warning rather than an error
        m1 = Module('m1', components=[Component('shared'), Component('c1')])
        m2 = Module('m2', components=[Component('shared')])

        with self.assertWarns(UserWarning):
            combined = m1 + m2

        self.assertEqual(
            ['c1', 'shared'], sorted(c.name for c in combined.components)
        )

    def test_module_addition_invalid(self):
        module = Module('m1')
        with self.assertRaises(TypeError):
            module + Component('c1')

    def test_mixture_addition(self):
        mixture = Mixture('test_mixture', components=[Component('c1')])
        module = Module('test_module', components=[Component('c2')])

        system = mixture + module

        # a new Mixture is returned; the original is untouched
        self.assertTrue(isinstance(system, Mixture))
        self.assertTrue(system is not mixture)
        self.assertEqual(1, len(mixture.components))
        self.assertEqual(2, len(system.components))
        self.assertTrue(system.get_component(name='test_module') is not None)

        # chaining works, which is the main way modules are used
        system = (
            mixture + module + Module('other', components=[Component('c3')])
        )
        self.assertEqual(3, len(system.components))
        self.assertEqual(1, len(mixture.components))

    def test_mixture_addition_clears_stale_crn(self):
        mixture = Mixture('test_mixture', components=[Component('c1')])
        mixture.compile_crn()
        self.assertTrue(mixture.crn is not None)

        # the compiled CRN does not describe the new system
        system = mixture + Module('test_module')
        self.assertTrue(system.crn is None)

    def test_module_compiles_in_mixture(self):
        module = Module(
            'test_module',
            components=[DNAassembly('A', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0, 'ktl': 0.01},
            species=[Species('extra')],
        )
        crn = (Mixture('test_mixture') + module).compile_crn()

        species = [str(s) for s in crn.species]
        self.assertTrue('dna_A' in species)
        self.assertTrue('rna_A' in species)
        self.assertTrue('protein_A' in species)
        # species added directly by the module reach the CRN
        self.assertTrue('extra' in species)

    def test_parameter_precedence(self):
        # component > module > mixture
        mixture = Mixture(
            'test_mixture', parameters={'ktx': 1.0, 'ktl': 0.01}
        )

        module = Module(
            'm',
            components=[DNAassembly('A', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0},
        )
        crn = (mixture + module).compile_crn()
        self.assertEqual(2.0, transcription_rate(crn, 'A'))

        # a parameter on the component itself wins over the module's
        module = Module(
            'm',
            components=[
                DNAassembly(
                    'A', promoter='P', rbs='B', parameters={'ktx': 3.0}
                )
            ],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0},
        )
        crn = (mixture + module).compile_crn()
        self.assertEqual(3.0, transcription_rate(crn, 'A'))

        # the mixture is used for parameters the module does not define
        module = Module(
            'm',
            components=[DNAassembly('A', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
        )
        crn = (mixture + module).compile_crn()
        self.assertEqual(1.0, transcription_rate(crn, 'A'))

    def test_modules_keep_separate_parameters(self):
        # two modules in one mixture compile with their own parameters
        mixture = Mixture(
            'test_mixture', parameters={'ktx': 1.0, 'ktl': 0.01}
        )
        m1 = Module(
            'm1',
            components=[DNAassembly('A', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0},
        )
        m2 = Module(
            'm2',
            components=[DNAassembly('B', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 9.0},
        )

        crn = (mixture + m1 + m2).compile_crn()
        self.assertEqual(2.0, transcription_rate(crn, 'A'))
        self.assertEqual(9.0, transcription_rate(crn, 'B'))

    def test_module_mechanism_override(self):
        # overriding a mechanism on one module does not affect another
        m1 = Module(
            'm1',
            components=[Component('c1')],
            mechanisms={'transcription': SimpleTranscription()},
        )
        m2 = Module(
            'm2',
            components=[Component('c2')],
            mechanisms={'transcription': SimpleTranscription()},
        )

        replacement = SimpleTranslation()
        m1.add_mechanism(replacement, 'transcription', overwrite=True)

        self.assertEqual(
            replacement.name, m1.mechanisms['transcription'].name
        )
        self.assertEqual(
            SimpleTranscription().name, m2.mechanisms['transcription'].name
        )
        # the override reaches the components at enumeration time
        self.assertEqual(
            replacement.name,
            m1.enumerate_components()[0].mechanisms['transcription'].name,
        )

    def test_nested_modules(self):
        # a Module is a Component, so modules can contain modules
        inner = Module(
            'inner',
            components=[DNAassembly('A', promoter='P', rbs='B')],
            mechanisms=txtl_mechanisms(),
            parameters={'ktx': 2.0},
        )
        outer = Module(
            'outer',
            components=[inner],
            parameters={'ktl': 0.01},
        )

        crn = (Mixture('test_mixture') + outer).compile_crn()
        self.assertEqual(2.0, transcription_rate(crn, 'A'))

    def test_nested_module_parameter_precedence(self):
        # inner module > outer module > mixture
        mixture = Mixture('test_mixture', parameters={'ktx': 1.0})

        def system(inner_parameters, outer_parameters):
            inner = Module(
                'inner',
                components=[DNAassembly('A', promoter='P', rbs='B')],
                mechanisms=txtl_mechanisms(),
                parameters=dict(inner_parameters, ktl=0.01),
            )
            outer = Module(
                'outer', components=[inner], parameters=outer_parameters
            )
            return (mixture + outer).compile_crn()

        # the inner module is closest to the component, so it wins
        crn = system({'ktx': 5.0}, {'ktx': 3.0})
        self.assertEqual(5.0, transcription_rate(crn, 'A'))

        # the outer module fills in what the inner one does not define
        crn = system({}, {'ktx': 3.0})
        self.assertEqual(3.0, transcription_rate(crn, 'A'))

        # the mixture is used when neither module defines the parameter
        crn = system({}, {})
        self.assertEqual(1.0, transcription_rate(crn, 'A'))
