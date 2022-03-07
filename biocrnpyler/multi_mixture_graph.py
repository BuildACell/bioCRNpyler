import copy
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
                  **kwargs):
        self.name = name
        self.idCounter = {}
        
        # adjacency dictionary for graph 
        self.mixture_graph = {}
        
        # maps mixtures to compartments 
        self.compartment_mixture_map = {} 
        self.compartment_name_map = {} 
        
        # list of mixtures in the graph 
        if mixtures is None:
            self.mixtures = []
        else:
            self.add_mixture(mixtures)
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
    

    def add_mixture(self, mixture, compartment = None):
        if isinstance(mixture, Mixture): 
            # Copying the input mixture 
            mixture_copy = copy.deepcopy(mixture)
        
            # Using given compartment, or creating a new one for the mixture.  
            if not isinstance(compartment, Compartment):
                if compartment == None:
                    compartment = Compartment(name = mixture.name + str(self.get_mixture_id_counter(mixture.name)))
                elif isinstance(compartment, str):
                    compartment = Compartment(name = compartment + str(self.get_mixture_id_counter(compartment)))
                elif isinstance(compartment, List):
                    raise ValueError("You provided a list for compartment when mixture is one item")
                else:
                    raise ValueError("You did not input a valid compartment. You need to input a Compartment object or string name, or nothing so MultiMixtureGraph can self-generate")
            elif compartment.name in self.compartment_mixture_map:
                raise ValueError(f"A compartment called {compartment.name} is already part of the MultiMixtureGraph.")
            
            # Add compartment to compartment map 
            self.compartment_mixture_map[compartment.name] = mixture_copy
            self.compartment_name_map[compartment.name] = compartment
            mixture_copy.compartment = compartment 
            
            # Updating mixture_graph and list of mixtures
            self.mixture_graph[mixture_copy] = [] # mixture will not have connections already. 
            self.mixtures.append(mixture_copy)
            
            # Returns new mixture and it's compartment name
            return mixture_copy, mixture_copy.compartment.name, mixture_copy.compartment
        
        elif isinstance(mixture, List):
            return self.add_mixtures(mixture)
        
        else:
            raise ValueError("You did not input a Mixture or list of Mixtures")
            
    def add_mixtures(self, mixtures, compartments= None):
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
    
        if not compartment_name_1 in self.compartment_name_map:
            raise ValueError("The first compartment you inputted was not added to the graph!") 
        if not compartment_name_2 in self.compartment_name_map:
            raise ValueError("The second compartment you inputted was not added to the graph!") 

        self.mixture_graph[self.compartment_mixture_map[compartment_name_1]].append(self.compartment_mixture_map[compartment_name_2])
        
        # adding relationships in compartments
        self.compartment_name_map[compartment_name_1].add_relationship(relationship_1, self.compartment_name_map[compartment_name_2])

        if relationship_2 is not None:
            self.compartment_name_map[compartment_name_2].add_relationship(relationship_2, self.compartment_name_map[compartment_name_1])
            self.mixture_graph[self.compartment_mixture_map[compartment_name_2]].append(self.compartment_mixture_map[compartment_name_1])


    def duplicate_structure(self, compartment_name, n, shared_compartments= [] ):
        
        if not compartment_name in self.compartment_name_map:
            raise ValueError("The compartment you are trying to duplicate is not in the graph yet!")
        added_mixtures = []
        added_compartment_names = []
        added_compartments = []
        for i in range(n):
            added_mixture, added_compartment_name, added_compartment = self.add_mixture(self.compartment_mixture_map[compartment_name])
            added_mixtures.append(added_mixture)
            added_compartment_names.append(added_compartment_name)
            added_compartments.append(added_compartment)
            
            # When we add_mixture above, the corresponding compartments and their mixtures are not created, so we create new ones for each here. The issue with this is that if these newly created mixtures also have other mixtures that need to be copied, they won't be. Should we implement this recusively so all layers get taken care of?
            new_compartment_dict = self.compartment_name_map[added_compartment_name].get_compartment_dict()
            old_compartment_dict = self.compartment_name_map[compartment_name].get_compartment_dict()
            
            for item in old_compartment_dict.keys() :
                compartment = old_compartment_dict[item];
                if not item in shared_compartments:
                    mixture_to_copy = self.compartment_mixture_map[compartment.name]
                    comp_to_cpy = copy.deepcopy(compartment)
                    comp_to_cpy.name = "temp"
                    mx, cmp_name, cmp = self.add_mixture(mixture_to_copy, comp_to_cpy)
                    cmp.name = cmp.name + str(self.get_mixture_id_counter(compartment_name))
                    self.compartment_name_map[cmp.name] = cmp
                    self.compartment_mixture_map[cmp.name] = mx
                    del self.compartment_mixture_map["temp"]
                    del self.compartment_name_map["temp"]
                else: 
                    cmp= compartment
                # cmp is the compartment that is related that is to be added 
            
                # Finds they relationships between each pair and adds connections 
                related_compartment = compartment
                other_comp_dict = related_compartment.get_compartment_dict()
                key_of_interest = ""
                for key in other_comp_dict.keys():
                    if self.compartment_name_map[compartment_name] is other_comp_dict[key]:
                        key_of_interest = key 
                self.connect(added_compartment.name, cmp.name, item + str(self.get_mixture_id_counter(str(added_compartment.name) + str(item))), str(self.get_mixture_id_counter(str(cmp.name) + str(key_of_interest))))

        return added_mixtures, added_compartment_names, added_compartments
    
    def remove_mixture(self,mixture):
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
    def print_graph(self):
        pass 
    def compile_crn(self,recursion_depth = None, initial_concentration_dict = None, return_enumerated_components = False,
        initial_concentrations_at_end = False, copy_objects = True, add_reaction_species = True) -> ChemicalReactionNetwork:
        crn_species = []
        crn_reactions = []
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
         
        self.crn = ChemicalReactionNetwork(crn_species, crn_reactions)
        return self.crn
        