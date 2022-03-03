
from unittest import TestCase
from biocrnpyler import  Mixture, Species, DNA, Reaction, ChemicalReactionNetwork, Component, SimpleTranscription, SimpleTranslation, GlobalMechanism, Compartment, MultiMixtureGraph

class TestMultiMixtureGraph(TestCase):

    def test_add_mixture(self):
        mmg = MultiMixtureGraph("name") 
        mixture1 = Mixture('test_mixture1')
        mixture1_cpy, compartment1_name, compartment1 = mmg.add_mixture(mixture1)

        # Checking that the internal list of mixtures is only the added mixture 
        self.assertEqual([mixture1_cpy], mmg.mixtures)
        self.assertEqual({compartment1_name: compartment1}, mmg.compartment_name_map)
        self.assertEqual({compartment1_name: mixture1_cpy}, mmg.compartment_mixture_map)
        self.assertEqual(mixture1_cpy.compartment, compartment1)
        # Checking that the mixture has been added to the dictionary correctly 
        self.assertEqual(False, mixture1 in mmg.mixture_graph)
        self.assertEqual(True, mixture1_cpy in mmg.mixture_graph)
        # Checking that the string representation of the mixture and its copy are the same 
        self.assertTrue(mixture1_cpy.name == mixture1.name)
        self.assertEqual(mmg.mixture_graph[mixture1_cpy], [])
        mmg2 = MultiMixtureGraph("name") 
        # Testing behavior with adding compartments 
        mixture2 = Mixture('test_mixture2')
        test_compartment = Compartment("c1")
        mixture2_cpy, compartment2_name, compartment2 = mmg2.add_mixture(mixture2, test_compartment)
        self.assertEqual(mixture2_cpy.compartment, test_compartment)
        mixture3 = Mixture('test_mixture3')
        test_compartment_string = "c2"
        mixture3_cpy, compartment3_name, compartment3 = mmg2.add_mixture(mixture3, compartment = test_compartment_string)
        self.assertEqual(mixture3_cpy.compartment.name[0:2], test_compartment_string)
        
        mixture_list = [mixture2, mixture3]
        mixture_copy_list, comp_name_list, comp_list = mmg.add_mixture(mixture_list)
        self.assertTrue(mixture_copy_list[0] in mmg.mixtures)
        self.assertTrue( mixture_copy_list[1] in mmg.mixtures)
        self.assertEqual(False, mixture2 in mmg.mixture_graph and mixture2 in mmg.mixture_graph)
        self.assertEqual(True, mixture_copy_list[0] in mmg.mixture_graph and mixture_copy_list[1] in mmg.mixture_graph)
        
   
                
        bad_compartment = mixture2
        with self.assertRaisesRegex(ValueError,"You did not input a valid compartment. You need to input a Compartment object or string name, or nothing so MultiMixtureGraph can self-generate"):
             mmg.add_mixture(mixture3, compartment = bad_compartment)
        bad_compartment2  = [Compartment("bc1"),Compartment("bc2")]
        with self.assertRaisesRegex(ValueError,"You provided a list for compartment when mixture is one item"):
             mmg.add_mixture(mixture3, compartment = bad_compartment2)
        
        bad_mixture = "I'm not a mixture!"
        with self.assertRaisesRegex(ValueError,"You did not input a Mixture or list of Mixtures"):
            mmg.add_mixture(bad_mixture)
            
    def test_add_mixtures(self):
        mmg = MultiMixtureGraph("name") 
        mixture1 = Mixture('test_mixture1')
        mixture1_cpy, compartment1_name, compartment1 = mmg.add_mixtures(mixture1)
        
        # Testing behavior with a single mixture 
        self.assertEqual([mixture1_cpy], mmg.mixtures)
        self.assertEqual({compartment1_name: compartment1}, mmg.compartment_name_map)
        self.assertEqual({compartment1_name: mixture1_cpy}, mmg.compartment_mixture_map)
        self.assertEqual(mixture1_cpy.compartment, compartment1)
        # Checking that the mixture has been added to the dictionary correctly 
        self.assertEqual(False, mixture1 in mmg.mixture_graph)
        self.assertEqual(True, mixture1_cpy in mmg.mixture_graph)
        # Checking that the string representation of the mixture and its copy are the same 
        self.assertTrue(mixture1_cpy.name == mixture1.name)
        self.assertEqual(mmg.mixture_graph[mixture1_cpy], [])
        
        # Testing behavior with multiple mixtures
        mmg = MultiMixtureGraph("name") 
        mixture2 = Mixture('test_mixture2')
        mixture3 = Mixture('test_mixture3')
        mixture_list = [mixture2, mixture3]
        
        mixture_list, compartment_name_list, compartment_list = mmg.add_mixture(mixture_list)
        for item in mixture_list:
            self.assertTrue(item in mmg.mixtures)
        self.assertTrue(compartment_name_list[0] in mmg.compartment_name_map)
        self.assertTrue(compartment_name_list[1] in mmg.compartment_name_map)
        self.assertTrue(compartment_name_list[0] in mmg.compartment_mixture_map)
        self.assertTrue(compartment_name_list[1] in mmg.compartment_mixture_map)
        self.assertTrue(mixture_list[0].compartment is compartment_list[0])
        self.assertTrue(mixture_list[1].compartment is compartment_list[1])
        self.assertFalse(mixture2 in mmg.mixture_graph and mixture2 in mmg.mixture_graph)
        self.assertTrue(mixture_list[0] in mmg.mixture_graph and mixture_list[1] in mmg.mixture_graph)
        
        # Testing error messages
        
        bad_compartment_list = [Compartment("c1"),Compartment("c2") ]
        bad_mixture_list = [Mixture("m1"), Mixture("m2"),Mixture("m3")]
        bad_mixture_list2 = Mixture("m4")
        
        with self.assertRaisesRegex(ValueError,"You did not input a list of compartments for each item in your list of mixtures"):
            mmg.add_mixtures(bad_mixture_list,bad_compartment_list)
            
        with self.assertRaisesRegex(ValueError,"You did not input a Mixture or list of Mixtures, or for 1 mixture, you had more than 1 compartment."):
            mmg.add_mixtures(bad_mixture_list2,bad_compartment_list) 
            
    def test_connect(self):

        mmg = MultiMixtureGraph("name") 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        mixture1_cpy, compartment1_name, compartment1 = mmg.add_mixture(mixture1)
        mixture2_cpy, compartment2_name, compartment2 = mmg.add_mixture(mixture2)
        mmg.connect(compartment1_name, compartment2_name, "internal", "external")
    
        # test adding with no mixtures in there
        self.assertEqual(True, mixture1_cpy in mmg.mixture_graph[mixture2_cpy])
        self.assertEqual(True, mixture2_cpy in mmg.mixture_graph[mixture1_cpy])
        self.assertTrue(compartment2 in  mixture1_cpy.compartment.get_compartment("internal"))
        self.assertTrue(compartment1 in mixture2_cpy.compartment.get_compartment("external"))
        
        # testing error messages
        
        with self.assertRaisesRegex(ValueError,"The first compartment you inputted was not added to the graph!"):
            mmg.connect(mixture1_cpy, mixture2_cpy, "internal", "external")
        with self.assertRaisesRegex(ValueError,"The second compartment you inputted was not added to the graph!"):
            mmg.connect(compartment1_name, mixture2_cpy, "internal", "external")
            
    def test_duplicate_structure(self):
        mmg = MultiMixtureGraph("name") 
        c1 = Compartment("test1")
        m1 = Mixture("mixture1")
        c2 = Compartment("test2")
        m2 = Mixture("mixture2")
        with self.assertRaisesRegex(ValueError,"The compartment you are trying to duplicate is not in the graph yet!"):
            mmg.duplicate_structure(c1.name, 1)
        
        added_mixt1, c1_name, c1_returned = mmg.add_mixture(m1, c1)
        added_mixt2, c2_name, c2_returned = mmg.add_mixture(m2, c2)
        
        mmg.connect(c1_name, c2_name, "internal","external") 
        self.assertTrue(len(mmg.mixtures) is 2)
        
        mix_list, comp_name_list, comp_list= mmg.duplicate_structure(c1_name, 2)
        
        self.assertTrue(len(mmg.mixtures) is 6)
        
        self.assertTrue(len(mmg.compartment_name_map) is 6)
       
        
        self.assertTrue(len(mix_list[0].compartment.compartment_dict.keys())is 1)
        for item in mix_list:
            self.assertTrue(len(item.compartment.compartment_dict.keys()) is 1)
            self.assertTrue("internal" in item.compartment.compartment_dict or "external" in item.compartment.compartment_dict)
        
        mmg2 = MultiMixtureGraph("name2") 
        c3 = Compartment("test3")
        m3 = Mixture("mixture3")
        c4 = Compartment("test4")
        m4 = Mixture("mixture4")
        
        added_mixt3, c3_name, c3_returned = mmg2.add_mixture(m3, c3)
        added_mixt4, c4_name, c4_returned = mmg2.add_mixture(m4, c4)
        mmg2.connect(c3_name, c4_name, "internal","external") 
        
        # do something where you have shared 
        mix_list2, comp_name_list2, comp_list2= mmg2.duplicate_structure(c3_name, 1, shared_compartments =["internal"])
       
        self.assertTrue(len(mmg2.mixtures) is 3)
        self.assertTrue(len(mmg2.compartment_name_map) is 3)
        for item in mmg2.mixtures:
            if "internal" in item.compartment.get_compartment_dict().keys():
                self.assertTrue(not item.compartment.compartment_dict["internal"] is c4 or not item.compartment.compartment_dict["internal"] is c3 )
           
      
        
    def test_remove_mixture(self): 
        pass
    
    def test_get_mixtures(self):
        mmg = MultiMixtureGraph("name") 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        
        mixtures_list, c_name_list, c_list = mmg.add_mixtures([mixture1, mixture2])
        
        self.assertTrue(mixtures_list[0] in mmg.get_mixtures())
        self.assertTrue(mixtures_list[1] in mmg.get_mixtures())
    def test_get_graph(self):
        mmg = MultiMixtureGraph("name") 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        g = {}
        
        
        mixtures_list, c_name_list, c_list= mmg.add_mixtures([mixture1, mixture2])
        g[mixtures_list[0]] = []
        g[mixtures_list[1]] = []
        self.assertEqual(mmg.get_graph(), g)
    
    def test_print_graph(self):
        pass 
    
    def test_compile_crn(self):
        mmg = MultiMixtureGraph("test")
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        
        # This is just to appease the test while we are figuring out the best way to check the species 
        s = Species("s")
        t = Species("t")
        r = Species("r")

        mixture1.add_species([s, t])
        mixture2.add_species([s, r])
        
        mx_cpy, comp_name, comp = mmg.add_mixture(mixture1)
        mx_cpy2, comp_name2, comp2 = mmg.add_mixture(mixture2)
        
        mmg.connect(comp_name, comp_name2, "diffusion1", "diffusion2")
        
        crn = mmg.compile_crn()

        #True answer
        s_m1 = Species("s", compartment = comp_name)
        s_m2 = Species('s', compartment = comp_name2)
        t_m1 = Species("t", compartment = comp_name)
        r_m2 = Species("r", compartment = comp_name2)

        assert len(crn.species) == 4
        assert s_m1 in crn.species
        assert s_m2 in crn.species
        assert t_m1 in crn.species
        assert r_m2 in crn.species

        #Assert proper renaming
        assert s not in crn.species
        assert t not in crn.species
        assert r not in crn.species
            
        
        