
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import Reaction, Species, MassAction, ChemicalComplex


class TestReaction(TestCase):

    def test_reaction_initialization(self):
        # warns if both input and output species are empty
        mak = MassAction(k_forward=0.1)
        with self.assertWarns(Warning):
            Reaction(inputs=[], outputs=[], propensity_type=mak)

        # test for invalid propensity type
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=[], propensity_type=Species)

        # input must be a valid species object
        with self.assertRaises(ValueError):
            Reaction(inputs=['a'], outputs=[], propensity_type=mak)
        # output must be a valid species object
        with self.assertRaises(ValueError):
            Reaction(inputs=[], outputs=['b'], propensity_type=mak)

        rxn = Reaction.from_mass_action(inputs=[], outputs=[], k_forward=0.1, k_reverse=1)
        # test whether the reaction is registered as reversible
        self.assertTrue(rxn.is_reversible)
        # test whether the reaction is registered as massaction
        self.assertTrue(isinstance(MassAction, rxn.propensity_type))

        # test ChemicalComplex inputs
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')
        chem_com_sp1 = ChemicalComplex(species=sp1, stoichiometry=2)
        chem_com_sp2 = ChemicalComplex(species=sp2, stoichiometry=1)
        Reaction(inputs=[chem_com_sp1], outputs=[chem_com_sp2], propensity_type=MassAction(k_forward=1))

        # test different input and output lists
        Reaction(inputs=[chem_com_sp1], outputs=[sp2], propensity_type=MassAction(k_forward=1))

        # mixing ChemicalComplex and Species is not allowed
        with self.assertRaises(ValueError):
            Reaction(inputs=[chem_com_sp1, sp2], outputs=[sp1], propensity_type=MassAction(k_forward=1))

    def test_reaction_equality(self):
        """test for the_equality operator"""
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')
        rxn1 = Reaction(inputs=[sp1, sp2], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=1))
        rxn2 = Reaction(inputs=[sp2, sp1], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=1))
        rxn3 = Reaction(inputs=[sp2, sp1], outputs=[sp2, sp2], propensity_type=MassAction(k_forward=10))
        rxn4 = Reaction(inputs=[sp2, sp1], outputs=[sp2], propensity_type=MassAction(k_forward=1))
        self.assertTrue(rxn1 == rxn2)
        self.assertFalse(rxn1 == rxn3)
        self.assertFalse(rxn1 == rxn4)

    def test_reaction_list_flattening(self):
        sp1 = Species(name='test_species_a')
        sp2 = Species(name='test_species_b')
        k_f = 1
        mak = MassAction(k_forward=k_f)
        rxn1 = Reaction.from_nested_list(inputs=[sp1, [sp1, sp2]], outputs=[[sp2, sp2], sp1], propensity_type=mak)
        rxn2 = Reaction(inputs=[sp1, sp1, sp2], outputs=[sp1, sp2, sp2], propensity_type=mak)
        self.assertTrue(rxn1 == rxn2)

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
