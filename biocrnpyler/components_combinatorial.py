#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from .species import Complex, ComplexSpecies, OrderedPolymerSpecies, PolymerConformation
from .component import Component
from itertools import permutations
import copy

class CombinatorialComplex(Component):
    """
    A class to represent a Complex of many Species which can bind together in many different ways.
    """
    
    def __init__(self, final_states, initial_states = None, intermediate_states = None, name = None, **keywords):
        """
        Binding reactions will be generated to form all the ComplexSpecies in final_states from all the species in initial_states 
        (or, if initial_states is None, from all the individual species inside each ComplexSpecies). Intermediate states restricts
        the binding reactions to only form species in this list. At a high level this generates the following reactions:
        
        If just final_states are given:
            final_states_internal_species <-[Combinatorial Binding]-> final_states
        
        if initial_states are given:
            intial_states <-[Combinatorial Binding]-> final_states

        if intermediate_states are given:
            final_states_internal_species <-[Combinatorial Binding]-> intermediate_states <-[Combinatorial Binding]->final_states

        if initial_states and intermediate_states are given:
            intial_states <-[Combinatorial Binding]-> intermediate_states <-[Combinatorial Binding]->final_states
        

        :param final_states: a single ComplexSpecies or a list of ComplexSpecies. 
        :param initial_states: a list of initial Species which are bound together to form the ComplexSpecies in final_states. 
                                If None defaults to the members of the ComplexSpecies in final_states.
        :param intermediate_states: a list of intermediate species formed when converting initial_states to final_states. 
                                    If None, all possible intermediate ComplexSpecies are enumerated.


        Example 1: final_states = ComplexSpecies([A, B, C]). initial_states = None, intermediate_states = None.
                   initial_states will default to A, B, C. All intermediate states [A, B], [A, C], [B, C] will be enumerated. 
                   This results in the 6 reversible reactions:
                   1. A + B <--> Complex([A, B])
                   2. A + C <--> Complex([A, C])
                   3. B + C <--> Complex([B, C])
                   4. Complex([A, B]) + C <--> Complex([A, B, C])
                   5. Complex([A, C]) + B <--> Complex([A, B, C])
                   6. Complex([B, C]) + A <--> Complex([A, B, C])

        Example 2: final_states = ComplexSpecies([A, B, C]), intial_states = [Complex([A, B]), Complex([A, C])], intermediate_states = None
                   This results in the reactions:
                   1. Complex([A, B]) + C <--> Complex([A, B, C])
                   2. Complex([A, C]) + B <--> Complex([A, B, C])
        
        Example 3: final_states = ComplexSpecies([A, B, C]). initial_states = None, intermaied_states = [Complex([A, B]), Complex([A, C])].
                   This results in reactions:
                   1. A + B <--> Complex([A, B])
                   2. A + C <--> Complex([A, C])
                   3. Complex([A, B]) + C <--> Complex([A, B, C])
                   4. Complex([A, C]) + B <--> Complex([A, B, C])

        Example 4: final_states = [Complex([A, A, B], Complex([A, B, B]))], initial_states = None, intermediate_states = None
                   This results in the reactions:
                   1. A + A <--> Complex([A, A])
                   2. Complex([A, A]) + B <--> Complex ([A, A, B])
                   3. B + B <--> Complex([B, B])
                   4. Complex([B, B]) + A <--> Complex ([A, B, B])
                   5. A + B <--> Complex([A, B])
                   6. Complex([A, B]) + A <--> Complex ([A, A, B])
                   7. Complex([A, B]) + B <--> Complex ([A, B, B])
        """

        #The order these run in is important!

        #1. set final_states
        self.final_states = final_states
        #2. set intermidiate_states
        self.intermediate_states = intermediate_states
        #3. set initial_states
        self.initial_states = initial_states

        #used to store combinations of species during update
        self.combination_dict = {}

        #Call super
        if name is None:
            name = ""
            for s in self.final_states:
                name += s.name+"_"
            name = name[:-1]
        super().__init__(name, **keywords)

    #Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        if isinstance(final_states, list):
            self._final_states = [self.set_species(s) for s in final_states]
        else:
            self._final_states = [self.set_species(final_states)]

        #all final_states must be ComplexSpecies
        if not all([isinstance(s, ComplexSpecies) for s in self.final_states]):
            raise ValueError(f"final_states must be a list of {ComplexSpecies} (or subclasses thereof). Recieved: {final_states}.")

        #Then create a list of all sub-species included in final_states Complexes
        self.sub_species = []
        for s in self.final_states:
            self.sub_species += s.species_set

    #Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        #set initial states
        if initial_states is None:
            self._initial_states = self.sub_species
        else:
            self._initial_states = [self.set_species(s) for s in initial_states]

            for s in self._initial_states:
                if not (s in self.sub_species or (isinstance(s, ComplexSpecies) and all([ss in self.sub_species for ss in s.species_set]))):
                    raise ValueError(f"Invalid initial species {s}; initial_states must either be contained in the final_states or a {ComplexSpecies} made of Species in the final_states.")

    #Intermediate states allows the user to restrict the complexes formed between the intial state and final state
    @property
    def intermediate_states(self):
        return self._intermediate_states
    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            self._intermediate_states = intermediate_states

            #All intermediate_states must be ComplexSpecies or OrderdedComplexSpecies
            if not all([isinstance(s, ComplexSpecies) for s in self._intermediate_states]):
                raise ValueError(f"intermediate must be a list of {ComplexSpecies} (or subclasses thereof). Recieved: {intermediate_states}.")
            #All intermediate_states must be made of sub_species
            for s in self._intermediate_states:
                intermediate_sub_species = s.species_set
                if not all([ss in self.sub_species for ss in intermediate_sub_species]):
                    raise ValueError(f"intermediate species {s} contains subspecies not in the final_states.")

    
    def compute_species_to_add(self, s0, sf):
        #Compute Species that need to be added to s0 to get the Complex sf

        if not isinstance(sf, ComplexSpecies):
            raise ValueError(f"sf must be a ComplexSpecies. Recieved {sf}")

        species_to_add = []
        for s in sf.species_set:
            if s == s0:
                s0_count = 1
            elif s in s0:
                s0_count = s0.species.count(s)
            else:
                s0_count = 0

            #Add the correct number of s to the list
            sf_count = sf.species.count(s)
            if s0_count < sf_count:
                species_to_add += (sf_count-s0_count)*[s]
            elif s0_count > sf_count:
                species_to_add = [] #In this case, s0 contains more stuff than sf, so nothing should be returned
                break
            else:
                pass #if they have the same number, do not add it

        #s0 contains more or different species than sf, return None
        if (not isinstance(s0, ComplexSpecies)) and s0 not in sf.species:
            species_to_add = None
        elif (isinstance(s0, ComplexSpecies) and not all([sf.species.count(s) >= s0.species.count(s) for s in s0.species_set])) and s0 not in sf.species:
            species_to_add = None
        elif len(species_to_add) == 0:
            species_to_add = None

        return species_to_add

    def get_combinations_between(self, s0, sf):
        """
        Returns all combinations of Species to create the Complex sf from s0.
        """

        species_to_add = self.compute_species_to_add(s0, sf)
        
        if species_to_add is None or len(species_to_add) == 0:
            return []
        else:
            #combinations (binder, bindee, complex_species) to be returned
            combinations = []
            #get all permutations of the species for different binding orders 
            perms = permutations(species_to_add, len(species_to_add))

            #Iterate through permutations
            for perm in perms:
                bindee = s0
                for s in perm:
                    binder = s
                    if isinstance(bindee, ComplexSpecies):
                        s_list = list(bindee.species)
                    else:
                        s_list = [bindee]

                    s_list.append(s)
                    cs = Complex(s_list)

                    #append the new combination
                    combinations.append((binder, bindee, cs))
                    #update bindee
                    bindee = cs

            return combinations
            

    def update_species(self):
        mech_b = self.get_mechanism("binding")
        species = []
        species_added_dict = {} #save which combinations have already been added
        self.combination_dict = {} #this should be recomputed every updated species

        #If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    #Get combinatorial species between s0 and si
                    if (s0, si) not in self.combination_dict:
                        self.combination_dict[s0, si] = self.get_combinations_between(s0, si)

                    #iterate through combinations of species between s0 and si
                    for binder, bindee, cs in self.combination_dict[s0, si]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

            for si in self.intermediate_states:
                for sf in self.final_states:
                    #Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = self.get_combinations_between(si, sf)

                    #iterate through combinations of species between si and sf
                    for binder, bindee, cs in self.combination_dict[si, sf]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

        #If there are no intermediate restrictions, compute combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    #Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = self.get_combinations_between(s0, sf)
                        
                    #iterate through combinations of species between s0 and sf
                    for binder, bindee, cs in self.combination_dict[s0, sf]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

        return list(set(species))

    def update_reactions(self):
        mech_b = self.get_mechanism("binding")
        reactions = []
        reactions_added_dict = {} #save which combinations have already been added in order to not add duplicates

        #If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    #Get combinatorial species between s0 and si
                    if (s0, si) not in self.combination_dict:
                        self.combination_dict[s0, si] = self.get_combinations_between(s0, si)

                    #iterate through combinations of species between s0 and si
                    for binder, bindee, cs in self.combination_dict[s0, si]:
                        if (binder, bindee, cs) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

            for si in self.intermediate_states:
                for sf in self.final_states:
                    #Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = self.get_combinations_between(si, sf)

                    #iterate through combinations of species between si and sf
                    for binder, bindee, cs in self.combination_dict[si, sf]:
                        if (binder, bindee, cs) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

        #If there are no intermediate restrictions, compute combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    #Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = self.get_combinations_between(s0, sf)
                        
                    #iterate through combinations of species between s0 and sf
                    for binder, bindee, cs in self.combination_dict[s0, sf]:
                        if (binder, bindee, cs) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(binder = binder, bindee = bindee, complex_species = cs, 
                            component = self, part_id = self.name+"_"+str(binder)+"_"+str(bindee))

        return reactions



class CombinatorialConformation(CombinatorialComplex):
    """
    A class to represent a PolymerConformation (or OrderedPolymerSpecies) with many internal Complexes which can bind and unbind in many different ways.
    """
    
    def __init__(self, final_states, initial_states, intermediate_states = None, name = None, **keywords):
        """
        Binding reactions will be generated to form all PolymerConformations in final_states 
        from all the OrderedPolymerSpecies or PolymerConformations in initial_states.
        (or, if initial_states is None, from all the individual Polymers inside each PolymerConformation). 
        Intermediate states restricts the binding reactions to first form OrderdPolymerSpecies or PolymerConformations in this list. 
        At a high level this generates the following reactions:

        intial_states <-[Combinatorial Binding]-> final_states

        if  intermediate_states are given:
            intial_states <-[Combinatorial Binding]-> intermediate_states <-[Combinatorial Binding]->final_states
        

        :param final_states: a single PolymerConformation/OrderedPolymerSpecies or a list of said classes.
        :param initial_states: a list of initial OrderedPolymerSpecies which are bound together to form the final_states.
        :param intermediate_states: a list of intermediate species formed when converting initial_states to final_states. 
                                    If None, all possible intermediate OrderedPolymerSpecies and PolymerConformations are enumerated.
        """
        super().__init__(final_states, initial_states = initial_states, intermediate_states = intermaied_states, name = name, **keywords)

    #Getters and setters
    #Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        if isinstance(final_states, list):
            self._final_states = [self.set_species(s) for s in final_states]
        else:
            self._final_states = [self.set_species(final_states)]

        #all final_states must be ComplexSpecies
        if not all([isinstance(s, PolymerConformation) or isinstance(s, OrderedPolymerSpecies) for s in self.final_states]):
            raise ValueError(f"final_states must be a list of PolymerConformation or OrderedPolymerSpecies. Recieved: {final_states}.")

        #Then create a list of all sub-polymers included in final_states Conformations
        self.sub_polymers = []
        for s in self.final_states:
            if isinstance(s, PolymerConformation):
                for p in s.polymers:
                    self.sub_polymers.append(p)
            elif isinstance(s, OrderedPolymerSpecies):
                self.sub_polymers.append(s)

        self.sub_polymers = list(set(sub_polymers))

    #Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        #set initial states
        if initial_states is None:
            self._initial_states = self.sub_species
        else:
            self._initial_states = [self.set_species(s) for s in initial_states]

            for s in self._initial_states:
                if not (s in self.sub_species or (isinstance(s, ComplexSpecies) and all([ss in self.sub_species for ss in s.species_set]))):
                    raise ValueError(f"Invalid initial species {s}; initial_states must either be contained in the final_states or a {ComplexSpecies} made of Species in the final_states.")

    #Intermediate states allows the user to restrict the complexes formed between the intial state and final state
    @property
    def intermediate_states(self):
        return self._intermediate_states
    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            self._intermediate_states = intermediate_states

            #All intermediate_states must be ComplexSpecies or OrderdedComplexSpecies
            if not all([isinstance(s, ComplexSpecies) for s in self._intermediate_states]):
                raise ValueError(f"intermediate must be a list of {ComplexSpecies} (or subclasses thereof). Recieved: {intermediate_states}.")
            #All intermediate_states must be made of sub_species
            for s in self._intermediate_states:
                intermediate_sub_species = s.species_set
                if not all([ss in self.sub_species for ss in intermediate_sub_species]):
                    raise ValueError(f"intermediate species {s} contains subspecies not in the final_states.")