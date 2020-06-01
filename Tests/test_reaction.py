
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Reaction, Species


class TestReaction(TestCase):

    def test_reaction_initialization(self):
        # warns if both input and output species are empty
        with self.assertWarns(Warning):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="massaction", propensity_params=None)

        # non-massaction propensities require propensity_params dict
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="not_massaction")

        # non-massaction propensities cannot be reversible
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="not_massaction", propensity_params={}, k_rev=1)

        # test under specified propensity parameters
        propensity_types = ["hillpositive", "hillnegative", "proportionalhillpositive", "proportionalhillnegative"]
        for propensity_type in propensity_types:
            with self.assertRaises(ValueError):
                Reaction(inputs=[], outputs=[], k=0.1, propensity_type=propensity_type, propensity_params={'s1': 0})

        # test when rate is missing from the propensity parameter dictionary for a general propensity
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type='general', propensity_params={})

        # test unknown propensity type
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, propensity_type="dummy", propensity_params={})

        # input must be a valid species object
        with self.assertRaises(ValueError):
            Reaction(inputs=['a'], outputs=[], k=0.1)
        # output must be a valid species object
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=['b'], k=0.1)

        # reaction rate coefficient must be larger than zero
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0)
        # reaction rate coefficient must be larger than zero
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=-1)

        rxn = Reaction(inputs=[], outputs=[], k=0.1, k_rev=1)
        # test whether the reaction is registered as reversible
        self.assertTrue(rxn.reversible)
        # test whether the reaction is registered as massaction
        self.assertTrue(rxn.propensity_type == 'massaction')

        # test overspecified mass action
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')
        with self.assertWarns(Warning):
            Reaction(inputs=[sp1], outputs=[sp2], propensity_type="massaction", propensity_params={}, k=0.1)

        # test whether the number of input and output coefficients match
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, input_coefs=[1,2])

        # test whether the number of input and output coefficients match
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], k=0.1, output_coefs=[1, 2])

        # TODO add test for the equality operator

    def test_complex_set_equality(self):
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')

        # test whether two reactions with the same species are equal
        rxn1 = Reaction(inputs=[sp1], outputs=[], k=0.1, input_coefs=[1])

        rtn = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs=rxn1.input_coefs,
                                            c2=rxn1.inputs, c2_coefs=rxn1.input_coefs)
        self.assertTrue(rtn)

        # test that two reactions have the two species with equal coefficients are equal
        rxn2 = Reaction(inputs=[sp1], outputs=[sp2], k=0.1, input_coefs=[1], output_coefs=[1])
        rtn = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs=rxn1.input_coefs,
                                            c2=rxn2.inputs, c2_coefs=rxn2.input_coefs)
        self.assertTrue(rtn)

        # test that two reactions have the two species with different coefficients are not equal
        rxn1 = Reaction(inputs=[sp1], outputs=[], k=0.1, input_coefs=[1])
        rxn2 = Reaction(inputs=[sp2], outputs=[], k=0.1, input_coefs=[2])

        rtn2 = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs=rxn1.input_coefs,
                                             c2=rxn2.inputs, c2_coefs=rxn2.input_coefs)
        self.assertFalse(rtn2)

        # test that two reactions with different species are not equal
        rxn1 = Reaction(inputs=[sp1,sp2], outputs=[], k=0.1, input_coefs=[1,2])
        rxn2 = Reaction(inputs=[sp2], outputs=[], k=0.1, input_coefs=[2])

        rtn3 = Reaction.complex_set_equality(c1=rxn1.inputs, c1_coefs=rxn1.input_coefs,
                                             c2=rxn2.inputs, c2_coefs=rxn2.input_coefs)
        self.assertFalse(rtn3)
