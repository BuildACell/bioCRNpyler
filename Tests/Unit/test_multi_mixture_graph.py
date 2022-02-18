from unittest import TestCase
from biocrnpyler import Mixture, Species, Reaction, ChemicalReactionNetwork, Compartment

class TestMultiMixtureGraph(TestCase):
    
    def test_name(self):
        # This need not be tested well; tested in other classes. 
        # Here, we just test the name is being added correctly
        mmg = MultiMixtureGraph() 
        name = "name"
        with self.assertRaisesRegex(ValueError,"MultiMixtureGraph name must be a string."):
            mmg.name(None)
        mmg.name(name)
        self.assertEqual(name, mmg.name)

    def test_add_mixture(self):
        mmg = MultiMixtureGraph() 
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
        self.assertTrue(mixture1_cpy.name == mixture.name)
        self.assertEqual(mmg.mixture_graph[mixture_1_cpy], [])
        mmg2 = MultiMixtureGraph() 
        # Testing behavior with adding compartments 
        mixture2 = Mixture('test_mixture2')
        test_compartment = Compartment("c1")
        mixture2_cpy, compartment2_name, compartment2 = mmg.add_mixture(mixture2, test_compartment)
        self.assertEqual(mixture2_cpy.compartment, test_compartment)
        mixture3 = Mixture('test_mixture3')
        test_compartment_string = "c2"
        mixture3_cpy, compartment3_name, compartment3 = mmg.add_mixture(mixture3, compartment = test_compartment_string)
        self.assertEqual(mixture3_cpy.compartment.name[0:2], test_compartment_string)
        
        mixture_copy_list = mmg.add_mixture(mixture_list)
        self.assertEqual([mixture1_cpy, mixture_copy_list[0], mixture_copy_list[1]], mmg.mixtures)
        self.assertEqual(False, mixture2 in mmg.mixture_graph and mixture2 in mmg.mixture_graph)
        self.assertEqual(True, mixture_copy_list[0] in mmg.mixture_graph and mixture_copy_list[1] in mmg.mixture_graph)
        
        with self.assertRaisesRegex(ValueError,"This compartment has a name that is already associated with a mixture. Rename the compartment"):
             mmg.add_mixture(mixture3, compartment = test_compartment)
                
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
        mmg = MultiMixtureGraph() 
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
        self.assertTrue(mixture1_cpy.name == mixture.name)
        self.assertEqual(mmg.mixture_graph[mixture_1_cpy], [])
        
        # Testing behavior with multiple mixtures
        mmg = MultiMixtureGraph() 
        mixture2 = Mixture('test_mixture2')
        mixture3 = Mixture('test_mixture3')
        mixture_list = [mixture2, mixture3]
        
        mixture_list, compartment_name_list, compartment_list = mmg.add_mixture(mixture_list)
        self.assertEqual(mixture_list, mmg.mixtures)
        self.assertTrue(compartment_name_list[0] in compartment_name_map)
        self.assertTrue(compartment_name_list[1] in compartment_name_map)
        self.assertTrue(compartment_name_list[0] in compartment_mixture_map)
        self.assertTrue(compartment_name_list[1] in compartment_mixture_map)
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

        mmg = MultiMixtureGraph() 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        mixture1_cpy, compartment1_name, compartment1 = mmg.add_mixture(mixture_1)
        mixture2_cpy, compartment2_name, compartment2 = mmg.add_mixture(mixture_2)
        mmg.connect(compartment1_name, compartment2_name, "internal", "external")
    
        # test adding with no mixtures in there
        [m1, m2, transp_m] = mmg.add_transport(mixture1, transport_mixture, mixture_2)
        self.assertEqual(True, mixture1_cpy in mmg.mixture_graph[mixture2_cpy])
        self.assertEqual(True, mixture2_cpy in mmg.mixture_graph[mixture1_cpy])
        self.assertEqual(True, mixture1_cpy.compartment.get_compartment("internal") is compartment2)
        self.assertEqual(True, mixture2_cpy.compartment.get_compartment("external") is compartment1)
        
        # testing error messages
        
        with self.assertRaisesRegex(ValueError,"The first compartment you inputted was not added to the graph!"):
            mmg.connect(mixture1_cpy, mixture2_cpy, "internal", "external")
        with self.assertRaisesRegex(ValueError,"The second compartment you inputted was not added to the graph!"):
            mmg.connect(compartment1_name, mixture2_cpy, "internal", "external")
            
    def test_remove_mixture(self): 
        pass
    
    def test_get_mixtures(self):
        mmg = MultiMixtureGraph() 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        
        mixtures_list, c_name_list, c_list = mmg.add_mixtures([mixture1, mixture2])
        
        assertEqual(mmg.get_mixtures,[mixtures_list[0], mixtures_list[1]])
    
    def test_get_graph(self):
        mmg = MultiMixtureGraph() 
        mixture1 = Mixture('test_mixture1')
        mixture2 = Mixture('test_mixture2')
        g = {}
        g[mixture1] = []
        g[mixture2] = []
        
        mixtures_list, c_name_list, c_list= mmg.add_mixtures[mixture1, mixture2]
        assertEqual(mmg.get_graph, g)
    
    def test_print_graph(self):
        pass 
    
    def test_compile_crn(self):
        pass
    # TODO -- will write once my compile_crn is finished