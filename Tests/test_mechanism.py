
#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase


class TestMechanism(TestCase):
    def test_mechanism_initialization(self):
        from biocrnpyler import Mechanism

        # test a warning for instantiation without a mechanism_type
        with self.assertWarns(Warning):
            Mechanism(name='test_mech', mechanism_type=None)

        # test a warning for instantiation without a mechanism_type
        with self.assertWarns(Warning):
            Mechanism(name='test_mech')

    def test_update_species(self):
        from biocrnpyler import Mechanism

        mech = Mechanism(name='test_mechanism', mechanism_type='dummy')
        # warning if update_species on a mechanism object
        with self.assertWarns(Warning):
            mech.update_species()

    def test_update_reactions(self):
        from biocrnpyler import Mechanism

        mech = Mechanism(name='test_mechanism', mechanism_type='dummy')

        # warning if update_reaction on a mechanism object
        with self.assertWarns(Warning):
            mech.update_reactions()
