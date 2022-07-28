import copy
import numpy as np
from typing import List, Union
from warnings import resetwarnings, warn

from .chemical_reaction_network import ChemicalReactionNetwork
from .component import Component
from .global_mechanism import GlobalMechanism
from .mechanism import Mechanism
from .parameter import ParameterDatabase
from .reaction import Reaction
from .species import Species
from .utils import remove_bindloc
from .compartments import Compartment
from .mixture import Mixture


class MultiMixtureGraph(object):
    def __init__(self, name="",
                  mixtures=None, 
                  parameters=None, 
                  parameter_file=None,
                  compartment_mixture_map=None,
                  compartment_name_map=None,
                  mixture_graph_start=None,
                  mixtures_no_copy=None,
                  **kwargs):
        self.name = name
        self.idCounter = {} 
        self.mixtures = []
        # if mixtures is None:
        #     self.mixtures = []
        if compartment_mixture_map is None:
            self.compartment_mixture_map = {}
        else: 
            self.compartment_mixture_map = compartment_mixture_map
        if compartment_name_map is None:
            self.compartment_name_map = {}
        else:
            self.compartment_name_map = compartment_name_map
        if mixture_graph_start is None:
            self.mixture_graph = {}
        else:
            self.mixture_graph = mixture_graph_start
        if mixtures_no_copy is None:
            mixtures_no_copy = []

        # mixtures to add without creating a copy 
        if not mixtures_no_copy is None :
            if not isinstance(mixtures_no_copy, List):
                raise TypeError("mixtures_no_copy must be a list!")
            for item in mixtures_no_copy:
                self.mixtures.append(item)
        if mixtures:
            self.add_mixture(mixtures)
        # else:
        #     self.mixtures = []
        for item in self.mixtures:
            if item not in self.mixture_graph.keys():
                self.mixture_graph[item] = []
        self.check_consistency()
    
    
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name: str):
        if name is None:
            raise TypeError("MultiMixtureGraph name must be a string.")
        elif type(name) is str:
            no_underscore_string = name.replace("_", "")
            if no_underscore_string.isalnum() and "__" not in name and name[len(name)-1] != "_" and name[0].isalpha():
                self._name = name
            else:
                raise ValueError(f"name attribute {name} must consist of letters, numbers, or underscores and cannot contained double underscores or begin/end with a special character.")
        else:
            raise ValueError('MultiMixtureGraph name must be a string.')
            
    


    def get_mixture_id_counter(self, mixture_name):
        if mixture_name in self.idCounter:
            self.idCounter[mixture_name] += 1
        else:
            self.idCounter[mixture_name] = 1

        return self.idCounter[mixture_name]
    
    def check_consistency(self):
        """
        Checks that across the different substructures in the multimixture graph,
        all things that should be consistent are consistent. 
        """

        for mixture in self.mixture_graph.keys():
            if mixture not in self.mixtures:
                # print(mixture)
                raise ValueError("mixture in mixture graph that isn't in mixtures")
        for mixture in self.mixtures:
            if mixture not in self.mixture_graph.keys():
                raise ValueError("mixture in mixtures and not in mixture graph")
    
        for compartment in self.compartment_mixture_map.keys():
            if self.compartment_mixture_map[compartment] not in self.mixtures:
                raise ValueError("compartment not in compartment_mixture_map")
        for compartment_name in self.compartment_mixture_map.keys():
            if compartment_name not in self.compartment_name_map.keys():
                raise ValueError("compartment not in compartment_name_map")
    

    def add_mixture(self, mixture, compartment = None):
        """
        Makes the MultiMixtureGraph aware of a mixture. It copies the inputted
        mixture so that it may be reused by the user. An associated compartment 
        is added unless passed in. It is added to relevant data structures in 
        the MultiMixtureGraph. 
        """
        if isinstance(mixture, Mixture): 
            mixture_copy = copy.deepcopy(mixture)
            if not isinstance(compartment, Compartment):
                if compartment == None:
                    compartment = Compartment(name = mixture.name + "_" + str(self.get_mixture_id_counter(mixture.name)))
                elif isinstance(compartment, str):
                    compartment = Compartment(name = compartment + "_" + str(self.get_mixture_id_counter(compartment)))
                elif isinstance(compartment, List):
                    raise ValueError("You provided a list for compartment when mixture is one item")
                else:
                    raise ValueError("You did not input a valid compartment. You need to input a Compartment object or string name, or nothing so MultiMixtureGraph can self-generate")
            elif compartment.name in self.compartment_mixture_map:
                raise ValueError(f"A compartment called {compartment.name} is already part of the MultiMixtureGraph.")
            
            self.compartment_mixture_map[compartment.name] = mixture_copy
            self.compartment_name_map[compartment.name] = compartment
            mixture_copy.compartment = compartment 
            self.mixture_graph[mixture_copy] = [] # mixture will not have connections already. 
            self.mixtures.append(mixture_copy)

            return mixture_copy, mixture_copy.compartment.name, mixture_copy.compartment
        elif isinstance(mixture, List):
            return self.add_mixtures(mixture)
        else:
            raise ValueError("You did not input a Mixture or list of Mixtures")
            
    def add_mixtures(self, mixtures, compartments= None):
        """
        Calls add_mixture. 
        """
        if isinstance(mixtures, Mixture) and (compartments is None or isinstance(compartments, Compartment)):
            return self.add_mixture(mixtures, compartments)
        elif isinstance(mixtures, List):
            if compartments is None or (isinstance(compartments, List) and len(compartments) is len(mixtures)):
                if compartments is None:
                    compartments = [None] * len(mixtures)
                mixtures_list = []
                compartments_list = []
                compartments_name_list = []
                for i in range(len(mixtures)):
                    mix, comp_name, comp = self.add_mixture(mixtures[i], compartments[i])
                    mixtures_list.append(mix)
                    compartments_list.append(comp)
                    compartments_name_list.append(comp_name)
                return mixtures_list, compartments_name_list, compartments_list 
            else: 
                raise ValueError("You did not input a list of compartments for each item in your list of mixtures") 
        else:
            raise ValueError("You did not input a Mixture or list of Mixtures, or for 1 mixture, you had more than 1 compartment.")
                
    
    def connect(self, compartment_name_1, compartment_name_2, relationship_1, relationship_2 = None):
        """
        Connects two things in the internal representation of the graph. 
        relationship_2 being None represents a one-directional edge. 
        """
        if not compartment_name_1 in self.compartment_name_map:
            raise ValueError("The first compartment you inputted was not added to the graph!") 
        if not compartment_name_2 in self.compartment_name_map:
            raise ValueError("The second compartment you inputted was not added to the graph!") 
        self.mixture_graph[self.compartment_mixture_map[compartment_name_1]].append(self.compartment_mixture_map[compartment_name_2])
        self.compartment_name_map[compartment_name_1].add_relationship(relationship_1, self.compartment_name_map[compartment_name_2])
        if relationship_2 is not None:
            self.compartment_name_map[compartment_name_2].add_relationship(relationship_2, self.compartment_name_map[compartment_name_1])
            self.mixture_graph[self.compartment_mixture_map[compartment_name_2]].append(self.compartment_mixture_map[compartment_name_1])

    
    @classmethod
    def combine_multi_mixture_graph(cls, mmgs = None, shared =None, new_name = "defaultName"):
        '''
        Combines two MultiMixtureGraphs into one and returns the new version. 
        Inputs:
            mmgs: list of MultiMixtureGraphs you want to combine
            shared: dictionary with the form {new_compartment_name_of_combined: 
            [list of compartment names, matching indicies of mmgs. Lists may 
            contain None to match indices ]}
        '''
        if mmgs is None:
            mmgs = []
        if shared is None:
            shared = {}

        new_mixtures =[]
        new_compartment_mixture_map = {}
        new_compartment_name_map = {} 
        new_mixture_graph = {}
        compartment_collision_counter = {}
        compartments_removed = {} 
        
        for new_compartment_name in shared.keys():
            keys_counter = {} 
            compartment_collision_counter[new_compartment_name] = 1
            comp_list = shared[new_compartment_name]
       
            base_compartment = mmgs[0].get_compartment(comp_list[0])
            base_dict = base_compartment.get_compartment_dict().copy()
            base_keys = base_dict.keys()
            base_mixture = (mmgs[0].get_compartment_mixture_map())[comp_list[0]]
            for base_key in base_keys:
                keys_counter[base_key] = 1
            to_remove = []
            for i, val in enumerate(comp_list[1::]):
                
                compartment = mmgs[i+1].get_compartment(val)
                to_remove.append(compartment.name)
                comp_dict = compartment.get_compartment_dict()
                for key in comp_dict.keys():
                  
                    if key in keys_counter.keys(): 
                        keys_counter[key] += 1
                        num = keys_counter[key]
                        new_key = key + "_" + str(num)
                        base_dict[new_key] = comp_dict[key]
                        keys_counter[new_key] = 1
                    else: 
                        base_dict[key] = comp_dict[key]
                        keys_counter[key] = 1
            c = Compartment(name = new_compartment_name, compartment_dict = base_dict) 
            new_compartment_name_map[new_compartment_name] = c 
            new_compartment_mixture_map[c.name] = base_mixture
            new_mixtures.append(base_mixture)
            base_mixture.compartment = c
            for thing in to_remove:
                compartments_removed[thing] = c
                
        # now we want to add the non-shared things. must check if they are in shared, but otherwise, add
        mmg_shared_comps = {}
    
        for k, mmg in enumerate(mmgs):
            mmg_shared_comps[k] = []
            for shared_key in shared.keys():
                mmg_shared_comps[k].append(shared[shared_key][k]) 
            temp_cmm = mmg.get_compartment_mixture_map()
            temp_cnm = mmg.get_compartment_name_map()
            for comp_name in temp_cnm.keys():
                if comp_name not in mmg_shared_comps[k]: 
                    if comp_name in compartment_collision_counter.keys():
                        num = compartment_collision_counter[comp_name]
                        compartment_collision_counter[comp_name] += 1
                        new_comp_name = comp_name + "_" + str(num)
               
                    else:
                        new_comp_name = comp_name
                    compartment_collision_counter[new_comp_name] = 1
                    new_compartment_name_map[new_comp_name] = temp_cnm[comp_name]
                    new_compartment_mixture_map[new_comp_name] = temp_cmm[comp_name]
                    new_mixtures.append(temp_cmm[comp_name]) 
        #  redirect 

  
        for name in new_compartment_name_map.keys():
            cmpt = new_compartment_name_map[name]
            d = cmpt.get_compartment_dict()
            for link in d.keys(): # this could be problematic because two could have the same name 
                if (d[link].name in compartments_removed.keys()):
                    d[link] = compartments_removed[d[link].name]
        for name in new_compartment_name_map:
            compartment = new_compartment_name_map[name]
            d = compartment.get_compartment_dict()
            mixture = new_compartment_mixture_map[name]
            new_mixture_graph[mixture] =[]
            for link in d.keys():
                link_mixture = new_compartment_mixture_map[(d[link].name)]
                new_mixture_graph[mixture].append(link_mixture)
                if(link_mixture in new_mixture_graph.keys()):
                    new_mixture_graph[link_mixture].append(mixture)
                else:
                    new_mixture_graph[link_mixture] = [mixture]
        
        return cls(name = new_name, mixtures_no_copy = new_mixtures,
         compartment_mixture_map = new_compartment_mixture_map, 
         compartment_name_map = new_compartment_name_map, 
         mixture_graph= new_mixture_graph)
    
    def remove_mixture(self,mixture):
        '''
        This should be a "backdoor" funciton that is not to be used regularly. 
        '''
        # TODO 
        pass
    def get_mixtures(self):
        return self.mixture_graph.keys()
    def get_graph(self):
        return self.mixture_graph
    def get_compartment_mixture_map(self):
        return self.compartment_mixture_map
    def get_compartment_name_map(self):
        return self.compartment_name_map
    def get_compartment(self, name):
        return self.compartment_name_map[name]
    def print_graph(self):
        '''
        This should give a way to visualize the graph. 
        '''
        # TODO 
        pass 


    def compile_crn(self,recursion_depth = None, initial_concentration_dict = None, 
    return_enumerated_components = False, initial_concentrations_at_end = False, 
    copy_objects = True, add_reaction_species = True, add_passive_diffusion = True, 
    df = 0.1, passive_diffusion_dict = {}) -> ChemicalReactionNetwork:
        crn_species = []
        crn_reactions = []
        mixture_crns = {}
        for mixture in self.mixtures:
            
            #TODO: Move this to a debug UTIL function?
            #Seems like below warnings might be too much as models are being built iteratively...
        
            # Checking for shared species
            #for connection in self.mixture_graph[mixture]:
                # added_species may not be the best method for this!!! REVISIT  
                #mixture_species = mixture.added_species 
                #connection_species = connection.added_species
                #intersection = [ms.comp_ind_eq(cs) for ms in mixture_species for cs in connection_species]
                #if len(intersection)<1:
                #    raise warn(f"There are no shared species between the connected Mixtures: {mixture.name} and {connection.name}.")
                
            temp = mixture.compile_crn()
            crn_species += temp.species
            crn_reactions+= temp.reactions
            mixture_crns[mixture.name] = temp
                                 
        if add_passive_diffusion:
            for key in self.mixture_graph.keys():
                values = self.mixture_graph[key]
                for value in values:
                    mixture_1_species = mixture_crns[key.name].species
                    mixture_2_species = mixture_crns[value.name].species
                    for species1 in mixture_1_species:
                        for species2 in mixture_2_species:
                            if species1.comp_ind_eq(species2):
                                diff_rate = df
                                for item in passive_diffusion_dict:
                                    if item.comp_ind_eq(species1):
                                        diff_rate = passive_diffusion_dict[item]
                                rxn = Reaction.from_massaction(inputs = [species1], 
                                outputs= [species2], k_forward = diff_rate, k_reverse = diff_rate)
                                crn_reactions.append(rxn)


        self.crn = ChemicalReactionNetwork(crn_species, crn_reactions)
        return self.crn
    
    @classmethod
    def create_diffusion_lattice(cls, n: int, diffusion_species: List[Species],
     mmg_name: str = "mmg"): 
        '''
        Generates a lattice of basic mixtures with diffusive species. Diffusion
        rates can be specified when the CRN is compiled. 
        '''
        new_mixtures = []
        compartment_names= []
        for i in range(n * n):
            new_mixtures.append(Mixture("m" + str(i+1)))
        mmg = cls(mmg_name)
        mmg.mixture_graph = {}
        mmg.compartment_mixture_map = {}
        mmg.mixtures = []
        mmg.compartment_name_map = {}
        for mixture in new_mixtures:
            for specs in diffusion_species:
                mixture.add_species(specs) # this should be added in zoila's component when integrated
            m, cn, c = mmg.add_mixture(mixture)
            compartment_names.append(cn)

        name_counter = np.zeros((n, n)) # counter for names of diffusion compartments. 
        # connecting horizontally 
        for i in range(n):
            for j in range(n-1):
                mmg.connect(compartment_names[n * i + j], 
                compartment_names[n * i + j + 1], "diffusion" + 
                str(name_counter[i, j]), "diffusion" + str(name_counter[i, j+1]))
                name_counter[i, j] += 1
                name_counter[i, j + 1] += 1
        # connecting vertically 
        for i in range(n):
            for j in range(n-1):
                mmg.connect(compartment_names[n * j + i], 
                compartment_names[n * (j+1) + i ], "diffusion" + 
                str(name_counter[j, i]), "diffusion" + str(name_counter[j +1, i]))
                name_counter[j, i] += 1
                name_counter[j + 1, i] += 1
        return mmg

   


    

        
        


        