
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import GlobalMechanism,Species, Complex, Mixture, MichaelisMenten, Deg_Tagged_Degredation, Degredation_mRNA_MM

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



    def test_rna_degredation_mm(self):

        M = Mixture(parameters = {"kdeg":1, "kb":1, "ku":1})

        rnaase = Species("P")
        rna = Species("X", material_type = "rna")

        #compare deg_Tagged_Degredation to MichaelisMenten
        MM = MichaelisMenten(name = "name", mechanism_type = "type")
        
        rdmm = Degredation_mRNA_MM(rnaase)
        #test default global mechanism parameters
        self.assertTrue(rdmm.recursive_species_filtering is True)
        self.assertTrue(rdmm.default_on is False)
        self.assertTrue("rna" in rdmm.filter_dict)
        self.assertTrue(rdmm.filter_dict["rna"] is True)
        self.assertTrue(len(rdmm.filter_dict)==2)

        #test default functionality on a tagged species
        species = rdmm.update_species(rna, M)
        rxns = rdmm.update_reactions(rna, M)
        species_mm = MM.update_species(Enzyme = rnaase, Sub = rna, Prod = None)
        rxns_mm = MM.update_reactions(Enzyme = rnaase, Sub = rna, Prod = None, kb=1, ku=1, kcat=1)

        #species
        self.assertTrue(species_mm == species)
        #reactions
        self.assertTrue(all([str(rxns[i]) == str(rxns_mm[i]) for i in range(len(rxns))]))

        #Overwrite default_on, recursive_species_filtering, and filter_dict
        rdmm2 = Degredation_mRNA_MM(rnaase, default_on = True, recursive_species_filtering = False, filter_dict = {"test":True})
        self.assertTrue(rdmm2.recursive_species_filtering is False)
        self.assertTrue(rdmm2.default_on is True)
        self.assertTrue("rna" not in rdmm2.filter_dict)
        self.assertTrue("test" in rdmm2.filter_dict and rdmm2.filter_dict["test"] is True)
        self.assertTrue(len(rdmm2.filter_dict)==1)

    def test_deg_tagged_degredation(self):

        M = Mixture(parameters = {"kdeg":1, "kb":1, "ku":1})

        protease = Species("P")
        tagged_protein = Species("X", attributes = ["degtagged"])
        untagged_protein = Species("Y")

        #compare deg_Tagged_Degredation to MichaelisMenten
        MM = MichaelisMenten(name = "name", mechanism_type = "type")
        
        dtd = Deg_Tagged_Degredation(protease)
        #test default global mechanism parameters
        self.assertTrue(dtd.recursive_species_filtering is False)
        self.assertTrue(dtd.default_on is False)
        self.assertTrue("degtagged" in dtd.filter_dict)
        self.assertTrue(dtd.filter_dict["degtagged"] is True)
        self.assertTrue(len(dtd.filter_dict)==1)

        #test default functionality on a tagged species
        species = dtd.update_species(tagged_protein, M)
        rxns = dtd.update_reactions(tagged_protein, M)
        species_mm = MM.update_species(Enzyme = protease, Sub = tagged_protein, Prod = None)
        rxns_mm = MM.update_reactions(Enzyme = protease, Sub = tagged_protein, Prod = None, kb=1, ku=1, kcat=1)

        #species
        self.assertTrue(species_mm == species)
        #reactions
        self.assertTrue(all([str(rxns[i]) == str(rxns_mm[i]) for i in range(len(rxns))]))

        #Overwrite default_on, recursive_species_filtering, and filter_dict
        dtd2 = Deg_Tagged_Degredation(protease, default_on = True, recursive_species_filtering = True, filter_dict = {"test":True})
        self.assertTrue(dtd2.recursive_species_filtering is True)
        self.assertTrue(dtd2.default_on is True)
        self.assertTrue("degtagged" not in dtd2.filter_dict)
        self.assertTrue("test" in dtd2.filter_dict and dtd2.filter_dict["test"] is True)
        self.assertTrue(len(dtd2.filter_dict)==1)