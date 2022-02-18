import copy
from typing import List, Union
from warnings import resetwarnings, warn
from .chemical_reaction_network import ChemicalReactionNetwork
from .reaction import Reaction
from .species import Species
from .compartments import Compartment

class MultiMixtureGraph(object):
    def __init__(self, name="",
                  mixtures=None, 
                  parameters=None, 
                  parameter_file=None,
                  **kwargs):
        self.name = name
        self.idCounter = 0
        
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
            
            
    def add_mixture(self, mixture, compartment = None):
        if isinstance(mixture, Mixture): 
            # Copying the input mixture 
            mixture_copy = copy.deepcopy(mixture)
        
            # Using given compartment, or creating a new one for the mixture.  
            if not isinstance(compartment, Compartment):
                if compartment == None:
                    compartment = Compartment(name = mixture.name + str(self.idCounter))
                    self.idCounter += 1            
                elif isinstance(compartment, str):
                    compartment = Compartment(name = compartment + str(self.idCounter))
                    self.idCounter += 1
                elif isinstance(compartment, List):
                    raise ValueError("You provided a list for compartment when mixture is one item")
                else:
                    raise ValueError("You did not input a valid compartment. You need to input a Compartment object or string name, or nothing so MultiMixtureGraph can self-generate")
            else:
                if compartment.name in compartment_name_map:
                    raise ValueError("This compartment has a name that is already associated with a mixture. Rename the compartment")
            
            # set global compartments 
            
            
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
                mixtures_list = []
                compartments_list = []
                compartments_name_list = []
                for i in len(mixtures):
                    mixtures_list.append(self.add_mixture(mixtures[i], compartments[i])[0])
                    compartments_list.append(self.add_mixture(mixtures[i], compartments[i])[1])
                    compartments_name_list.append(self.add_mixture(mixtures[i], compartments[i])[2])
                return mixtures_list, compartments_name_list, compartments_list 
            else: 
                raise ValueError("You did not input a list of compartments for each item in your list of mixtures") 
        else:
            raise ValueError("You did not input a Mixture or list of Mixtures, or for 1 mixture, you had more than 1 compartment.")
                
    

    # use keys to define relationship between mix1 and mix2 
    # instead of mixtures, use strings as keys which are stored in a dictionary 
    def connect(self, compartment_name_1, compartment_name_2, key_1, key_2):
    
        if not compartment_name_1 in self.compartment_name_map:
            raise ValueError("The first compartment you inputted was not added to the graph!") 
        if not compart2ent_name_1 in self.compartment_name_map:
            raise ValueError("The second compartment you inputted was not added to the graph!") 

        self.mixture_graph[compartment_mixture_graph[compartment_name_1]].append(compartment_mixture_graph[compartment_name_2])
        self.mixture_graph[compartment_mixture_graph[compartment_name_2]].append(compartment_mixture_graph[compartment_name_1])
        
        # adding relationships in compartments
        self.compartment_name_map[compartment_name_1].add_compartment(key_1, self.compartment_name_map[compartment_name_2])
        self.compartment_name_map[compartment_name_2].add_compartment(key_2, self.compartment_name_map[compartment_name_1])

    def duplicate_structure(self, compartment_name, shared_compartments, n ):
        
        if not compartment_name in compartment_name_map:
            raise ValueError("The compartment you are trying to duplicate is not in the graph yet!")
            
        for i in range(n):
            added_mixture, added_compartment_name, added_compartment = add_mixture(self.compartment_mixture_map[compartment_name])
            
            # When we add_mixture above, the compartment connections are created, but the corresponding 
            # mixtures are not, so we create new ones for each here
            # --- 
            # The issue with this is that if these newly created mixtures also have other mixtures that need to be copied, they won't be
            # This originates from the fact that compartments don't know the mixtures they relate to. 
            #^ This would be an easy change, but I'm not sure if it defeats the purpose of having a compartment class. Why don't we just store this 
            # local map in the mmg? 
            new_compartment_dict = self.compartment_name_map[added_compartment_name].get_compartment_dict()
            old_compartment_dict = self.compartment_name_map[compartment_name].get_compartment_dict()
            for item in new_compartment_dict :
                # making a new mixture by copying 
                if not item in shared_compartments:
                    # there could be many compartments under the same label 
                    for compartment in new_compartment_dict[item]:
                        mixture_to_copy = self.compartment_mixture_map[compartment]
                        mx, cmp_name, cmp = self.add_mixture(mixture_to_copy, compartment)
                else: 
                    added_compartment.set_compartment(item, old_compartment_dict[item])
               
                
                # getting keys for connection 
                # gets the compartment of the original, and gets what it connected to with this keyword
                # then, searches for what keyword had the original compartment in that corresponding compartment
            
                for related_compartment in new_compartment_dict[item]:
                    other_comp_dict = self.compartment_name_map[related_compartment].get_compartment_dict()
                    key_of_interest = ""
                    for key in other_comp_dict:
                        if compartment_name in other_comp_dict[key]:
                            key_of_interest = key 
                            break; 
                    self.connect(added_compartment, related_compartment, item, key_of_interest)
            return added_mixture, added_compartment_name, added_compartment
    
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
        for mixture in self.mixture_graph.keys():
            
            # Checking for shared species
            for connection in self.mixture_graph[mixture]:
                # added_species may not be the best method for this 
                mixture_species = mixture.added_species 
                connection_species = connection.added_species
                if not bool((set(mixture_species)).intersection(set(connection_species))):
                    raise ValueError("There are no shared species between at least one of the mixtures you have connected in the graph!")

            
            temp = mixture.compile_crn()
            specs = temp.species
            rxns = temp.reactions
            
            for item in rxns:
                if isinstance(item, Reaction):
                    crn_reactions.append(item)
            # adding the compartment to each species
            for species in specs:
                species.compartment = mixture.compartment
            crn_species.append(specs) 
         
        self.crn = ChemicalReactionNetwork(crn_species, crn_reactions)
        return self.crn
        