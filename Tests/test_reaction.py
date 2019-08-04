from unittest import TestCase


class TestReaction(TestCase):

    def test_reaction_initialization(self):
        from biocrnpyler import Reaction

        with self.assertWarns(Warning):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="massaction", propensity_params=0)

        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="not_massaction", k_rev=1)

        with self.assertRaises(ValueError):
            prop_types = ["hillpositive", "hillnegative", "proportionalhillpositive", "proportionalhillnegative"]
            for prop_type in prop_types:
                Reaction(inputs=[], outputs=[], k=0.1, propensity_type=prop_type, propensity_params={'s1': 0})

        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="dummy")

        with self.assertRaises(ValueError):
            Reaction(inputs=['a'], outputs=[], k=0.1)
            Reaction(inputs=[], outputs=['b'], k=0.1)

        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0)
            Reaction(inputs=[], outputs=[], k=-1)

        rxn = Reaction(inputs=[], outputs=[], k=0.1, k_rev=1)
        self.assertTrue(rxn.reversible)

        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, input_coefs=[1,2])

        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, output_coefs=[1, 2])

        # TODO add test for the equality operator

    def test_complex_set_equality(self):
        from biocrnpyler import Reaction
        from biocrnpyler import Species

        rxn1 = Reaction(inputs=[], outputs=[], k=0.1)
        rxn2 = Reaction(inputs=[], outputs=[], k=0.1)

        rtn = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs= rxn1.input_coefs,
                                            c2=rxn2.inputs, c2_coefs=rxn2.output_coefs)
        self.assertTrue(rtn)

        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')

        rxn1 = Reaction(inputs=[sp1], outputs=[], k=0.1)
        rxn2 = Reaction(inputs=[sp2], outputs=[], k=0.1)

        rtn2 = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs=rxn1.input_coefs,
                                             c2=rxn2.inputs, c2_coefs=rxn2.output_coefs)
        self.assertFalse(rtn2)
