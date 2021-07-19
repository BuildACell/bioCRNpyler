
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Mechanism


class TestMechanism(TestCase):

    def setUp(self) -> None:
        self.mech_name = 'test_mechanism'
        self.mechanism = Mechanism(name=self.mech_name, mechanism_type='dummy')

    def test_mechanism_initialization(self):

        # test a warning for instantiation without a mechanism_type
        with self.assertWarnsRegex(Warning, f'Mechanism {self.mech_name} instantiated without a type.'):
            Mechanism(name=self.mech_name, mechanism_type=None)

        # test a warning for instantiation without a mechanism_type
        with self.assertWarnsRegex(Warning, f'Mechanism {self.mech_name} instantiated without a type.'):
            Mechanism(name=self.mech_name)

    def test_update_species(self):
        # warning if update_species on a mechanism object
        with self.assertWarnsRegex(Warning, f"Default Update Species Called for Mechanism = {self.mech_name}."):
            self.mechanism.update_species()

    def test_update_reactions(self):
        # warning if update_reaction on a mechanism object
        with self.assertWarnsRegex(Warning, f'Default Update Reactions Called for Mechanism = {self.mech_name}.'):
            self.mechanism.update_reactions()

    def test_repr(self):
        # test that the __repr__ returns the mechanism name
        self.assertEqual(self.mech_name, str(self.mechanism))
