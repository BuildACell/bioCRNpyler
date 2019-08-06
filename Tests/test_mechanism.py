from unittest import TestCase


class TestMechanism(TestCase):
    def test_mechanism_initialization(self):
        from biocrnpyler import Mechanism

        with self.assertWarns(Warning):
            Mechanism(name='test_mech', mechanism_type=None)

        with self.assertWarns(Warning):
            Mechanism(name='test_mech')

    def test_update_species(self):
        from biocrnpyler import Mechanism

        mech = Mechanism(name='test_mechanism', mechanism_type='dummy')

        with self.assertWarns(Warning):
            mech.update_species()

    def test_update_reactions(self):
        from biocrnpyler import Mechanism

        mech = Mechanism(name='test_mechanism', mechanism_type='dummy')

        with self.assertWarns(Warning):
            mech.update_reactions()
