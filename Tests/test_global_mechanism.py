
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import GlobalMechanism,Species, Complex

class TestGlobalMechanism(TestCase):

    def setUp(self) -> None:
        self.mech_name = 'test_mechanism'
        self.mechanism = GlobalMechanism(name=self.mech_name, mechanism_type='dummy')

    def test_update_species(self):
        s = Species("s")
        params = {}
        # GlobalMechanism returns empty list of species by default
        self.assertTrue(len(self.mechanism.update_species(s, params)) == 0)

    def test_update_reactions(self):
        s = Species("s")
        params = {}
        # GlobalMechanism returns empty list of reactions by default
        self.assertTrue(len(self.mechanism.update_reactions(s, params)) == 0)

    def test_default_filtering(self):
        #No filter dictionary used in these tests
        s1 = Species("s1", material_type = "m1", attributes = ["a1"])
        s2 = Species("s2", material_type = "m2", attributes = ["a2"])
        c1 = Complex([s1, s2], name = "c1")
        c2 = Complex([c1, s2], name = "c2")

        #Always ON for all species
        mech_default_on = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", default_on = True)
        self.assertTrue(False not in [mech_default_on.apply_filter(s) for s in [s1, s2, c1, c2]])

        #Always OFF for all species
        mech_default_off = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", default_on = False)
        self.assertTrue(True not in [mech_default_off.apply_filter(s) for s in [s1, s2, c1, c2]])


    def test_filter_dictionary(self):
        #Test default filters via name, attributes, and materials with recursive_species_filtering = False

        s1 = Species("s1", material_type = "m1", attributes = ["a1"])
        s2 = Species("s2", material_type = "m2", attributes = ["a2"])
        c1 = Complex([s1, s2], name = "c1")
        c2 = Complex([c1, s2], name = "c2")


        fd = {"s1":False} #Filter based on name
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = True, filter_dict = fd, recursive_species_filtering = False)
        self.assertFalse(mech_default_on_fds1.apply_filter(s1))
        self.assertTrue(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1))
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))

        fd = {"m1":False} #Filter based on material
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = True, filter_dict = fd, recursive_species_filtering = False)
        self.assertFalse(mech_default_on_fds1.apply_filter(s1))
        self.assertTrue(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1))
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))

        fd = {"a1":False} #Filter based on attribute
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = True, filter_dict = fd, recursive_species_filtering = False)
        self.assertFalse(mech_default_on_fds1.apply_filter(s1))
        self.assertTrue(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1)) #attributes are NOT inherited through ComplexSpecies
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))


    def test_recursive_filtering(self):
        #Test default off via name, attributes, and materials with recursive_species_filtering = True

        s1 = Species("s1", material_type = "m1", attributes = ["a1"])
        s2 = Species("s2", material_type = "m2", attributes = ["a2"])
        c1 = Complex([s1, s2], name = "c1")
        c2 = Complex([c1, s2], name = "c2")

        fd = {"s1":True} #Filter based on name
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = False, filter_dict = fd, recursive_species_filtering = True)
        self.assertTrue(mech_default_on_fds1.apply_filter(s1))
        self.assertFalse(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1))
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))

        fd = {"m1":True} #Filter based on material
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = False, filter_dict = fd, recursive_species_filtering = True)
        self.assertTrue(mech_default_on_fds1.apply_filter(s1))
        self.assertFalse(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1))
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))

        fd = {"a1":True} #Filter based on attribute
        mech_default_on_fds1 = GlobalMechanism(name = self.mech_name, mechanism_type = "dummy", 
            default_on = False, filter_dict = fd, recursive_species_filtering = True)
        self.assertTrue(mech_default_on_fds1.apply_filter(s1))
        self.assertFalse(mech_default_on_fds1.apply_filter(s2))
        self.assertTrue(mech_default_on_fds1.apply_filter(c1)) #attributes are not inherited through ComplexSpecies, but contained inside
        self.assertTrue(mech_default_on_fds1.apply_filter(c2))


