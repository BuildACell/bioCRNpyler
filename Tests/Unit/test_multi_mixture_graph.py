
from unittest import TestCase
from biocrnpyler import  Mixture, Species, DNA, Reaction, ChemicalReactionNetwork, Component, SimpleTranscription, SimpleTranslation, GlobalMechanism, Compartment, MultiMixtureGraph

class TestMultiMixtureGraph(TestCase):
    def test_name(self):
        mmg = MultiMixtureGraph("test_my_name")
        self.assertEqual(mmg.name, "test_my_name")

    def test_check_consistency(self):
        """
        The point of this function is checking if any mixture or compartment
        is "unregistered". This is expecially important when merging two 
        multimixture graphs. 
        """
        mixture1 = Mixture("m1")
        mixture2 = Mixture("m2")
        mmg = MultiMixtureGraph("name")
        m1, cn1, c1 = mmg.add_mixture(mixture1)

        try:
            mmg.check_consistency()
        except ValueError as exc:
            assert False,  f"{exc}"
        
        mmg.mixtures.append(mixture2)
        try:
            mmg.check_consistency()
        except ValueError as exec:
            assert True, f"{exc}"
        mmg.mixtures.remove(mixture2)
        c = Compartment("dummy_compartment")
        mmg.compartment_mixture_map[c.name] = mixture2
        try:
            mmg.check_consistency()
        except ValueError as exec:
            assert True, f"{exc}"

        mmg2 = MultiMixtureGraph("combine_mmg1")
        m1, cn1, c1 = mmg2.add_mixture(mixture1)
        mmg3 = MultiMixtureGraph("combine_mmg2")
        m2, cn2, c2 = mmg3.add_mixture(mixture2)

        combined = MultiMixtureGraph.combine_multi_mixture_graph([mmg2, mmg3], shared = {}, new_name = "combined")
        try:
            combined.check_consistency()
        except ValueError as exc:
            assert False,  f"{exc}"

    def test_add_mixture(self):
        mmg = MultiMixtureGraph("test_add_mixture") 
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

        mmg2 = MultiMixtureGraph("test_add_mixture_2") 
        
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
        
        #Test naming convention
        m1 = Mixture("A")
        m2 = Mixture("B")
        m1_cpy1, name1, c1 = mmg.add_mixture(m1)
        m1_cpy2, name2, c2 = mmg.add_mixture(m1)
        m2_cpy1, name3, c3 = mmg.add_mixture(m2)
        m2_cpy2, name4, c4 = mmg.add_mixture(m2)
        
        assert name1 == m1.name + "_1"
        assert name2 == m1.name + "_2"
        assert name3 == m2.name + "_1"
        assert name4 == m2.name + "_2"
                
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
        mixture = Mixture('M')

        mixture_cpy1, compartment1_name, compartment1 = mmg.add_mixture(mixture)
        mixture_cpy2, compartment2_name, compartment2 = mmg.add_mixture(mixture)
        mixture_cpy3, compartment3_name, compartment3 = mmg.add_mixture(mixture)

        #Bidirectional connection
        mmg.connect(compartment1_name, compartment2_name, "internal", "external")
        
        # test adding with no mixtures in there
        self.assertTrue(len(mmg.mixture_graph[mixture_cpy1]) == len(mmg.mixture_graph[mixture_cpy2]) == 1)
        self.assertEqual(True, mixture_cpy1 in mmg.mixture_graph[mixture_cpy2])
        self.assertEqual(True, mixture_cpy2 in mmg.mixture_graph[mixture_cpy1])
        self.assertTrue(compartment2 ==  mixture_cpy1.compartment.get_compartment("internal"))
        self.assertTrue(compartment1 == mixture_cpy2.compartment.get_compartment("external"))

        #1 directional connection
        mmg.connect(compartment3_name, compartment1_name, "oneway")
        self.assertTrue(len(mmg.mixture_graph[mixture_cpy1]) == len(mmg.mixture_graph[mixture_cpy2]) == len(mmg.mixture_graph[mixture_cpy3]) == 1)
        self.assertTrue(mixture_cpy1 in mmg.mixture_graph[mixture_cpy3])
        self.assertFalse(mixture_cpy3 in mmg.mixture_graph[mixture_cpy1])

        
        # testing error messages
        with self.assertRaisesRegex(ValueError,"The first compartment you inputted was not added to the graph!"):
            mmg.connect(mixture_cpy1, mixture_cpy2, "internal", "external")
        with self.assertRaisesRegex(ValueError,"The second compartment you inputted was not added to the graph!"):
            mmg.connect(compartment1_name, mixture_cpy2, "internal", "external")

    def test_combine_multi_mixture_graph(self):
        mixture1 = Mixture("m1")
        mixture2 = Mixture("m2")

        mmg1 = MultiMixtureGraph("combine_mmg1")
        m1, cn1, c1 = mmg1.add_mixture(mixture1)
        mmg2 = MultiMixtureGraph("combine_mmg2")
        m2, cn2, c2 = mmg2.add_mixture(mixture2)
        m3, cn3, c3  = mmg2.add_mixture(mixture1)
        mmg2.connect(cn2, cn3, "internal", "external")
        shared = {"shared_mixture": [cn1, cn3]}
        combined = MultiMixtureGraph.combine_multi_mixture_graph([mmg1, mmg2], shared, new_name = "combined")

        self.assertTrue(combined.name == "combined")
        self.assertTrue(len(combined.mixtures) == 2)
        self.assertTrue(m1 in combined.mixtures )
        self.assertTrue(m2 in combined.mixtures )
        self.assertTrue(len(m1.compartment.get_compartment_dict().keys()) == 1)
        self.assertTrue("external" in m1.compartment.get_compartment_dict().keys())


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
        
        crn = mmg.compile_crn(add_passive_diffusion = False)

        #True answer
        s_m1 = Species("s", compartment = comp_name)
        s_m2 = Species('s', compartment = comp_name2)
        t_m1 = Species("t", compartment = comp_name)
        r_m2 = Species("r", compartment = comp_name2)

        assert len(crn.species) == 4
        assert len(crn.reactions) == 0
        assert s_m1 in crn.species
        assert s_m2 in crn.species
        assert t_m1 in crn.species
        assert r_m2 in crn.species

        #Assert proper renaming
        assert s not in crn.species
        assert t not in crn.species
        assert r not in crn.species
        


        # Testing adding passive diffusion 
        mmg2 = MultiMixtureGraph("test2")
        m1 = Mixture('m1')
        m2= Mixture('m2')

        m1.add_species([s, t])
        m2.add_species([s, r])

        # default is that passive diffusion is true
        crn2 = mmg.compile_crn(passive_diffusion_dict = {s: 0.2})
 
        assert len(crn2.reactions) > 0
        assert crn2.reactions[0].propensity_type.propensity_dict['parameters']['k_forward'] == 0.2
    
def test_create_diffusion_lattice():
    lattice = MultiMixtureGraph.create_diffusion_lattice(n = 3, 
    diffusion_species = [Species("spec")], mmg_name = "lattice")

    assert len(lattice.mixtures) is 9

    connections  = {}
    for mixture in lattice.mixtures:
        n = len(lattice.mixture_graph[mixture])
        if n in connections.keys():
            connections[n] +=1 
        else: 
            connections[n] = 1
    assert connections[2] is 4
    assert connections[3] is 4
    assert connections[4] is 1 








        




        
       



            
        
        