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
    
    def __init__(self, final_states, initial_states = None, intermediate_states = None, excluded_states = None, name = None, **keywords):
        """
        Binding reactions will be generated to form all the ComplexSpecies in final_states from all the species in initial_states 
        (or, if initial_states is None, from all the individual species inside each ComplexSpecies). Intermediate states restricts
        the binding reactions to only form species in this list. Excluded states are not allowed to be reactants or products. 
        At a high level this generates the following reactions:
        
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
        :param intermediate_states: a list of intermediate ComplexSpecies formed when converting initial_states to final_states. 
                                    If None: all possible intermediate ComplexSpecies are enumerated.
        :param excluded_states: a list of ComplexSpecies which are NOT allowed to form when converting initial states to final states.
                                If None: no ComplexSpecies are excluded.


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
        #2. set initial_states
        self.initial_states = initial_states
        #3. set intermidiate_states
        self.intermediate_states = intermediate_states
        #4. set excluded_states
        self.excluded_states = excluded_states

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
        final_states = list(self.set_species(final_states))

        #all final_states must be ComplexSpecies
        if not all([isinstance(s, ComplexSpecies) for s in final_states]):
            raise ValueError(f"final_states must be a list of {ComplexSpecies} (or subclasses thereof). Recieved: {final_states}.")

        self._final_states = final_states

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
            initial_states = list(self.set_species(initial_states))
            for s in initial_states:
                if not (s in self.sub_species or (isinstance(s, ComplexSpecies) and all([ss in self.sub_species for ss in s.species_set]))):
                    raise ValueError(f"Invalid initial species {s}; initial_states must either be contained in the final_states or a {ComplexSpecies} made of Species in the final_states.")
            self._initial_states = initial_states

    #Intermediate states allows the user to restrict the complexes formed between the intial state and final state
    @property
    def intermediate_states(self):
        return self._intermediate_states
    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            intermediate_states = list(self.set_species(intermediate_states))

            #All intermediate_states must be ComplexSpecies or OrderdedComplexSpecies
            if not all([isinstance(s, ComplexSpecies) for s in intermediate_states]):
                raise ValueError(f"intermediate must be a list of {ComplexSpecies} (or subclasses thereof). Recieved: {intermediate_states}.")
            #All intermediate_states must be made of sub_species
            for s in intermediate_states:
                intermediate_sub_species = s.species_set
                if not all([ss in self.sub_species for ss in intermediate_sub_species]):
                    raise ValueError(f"intermediate species {s} contains subspecies not in the final_states.")

            self._intermediate_states = intermediate_states

    #Excluded states allows the user to exclude specific Species from being enumerated
    @property
    def excluded_states(self):
        return self._excluded_states
    @excluded_states.setter
    def excluded_states(self, excluded_states):
        if excluded_states is None:
            self._excluded_states = []
        else:
            self._excluded_states = list(self.set_species(excluded_states))

    def compute_species_to_add(self, s0, sf):
        #Compute Species that need to be added to s0 to get the Complex sf

        if not isinstance(sf, ComplexSpecies):
            raise ValueError(f"sf must be a ComplexSpecies. Recieved {sf}")

        species_to_add = []
        for s in sf.species_set:
            if s.monomer_eq(s0): #this is used instead of == to deal with the potential for different parents
                s0_count = 1
            elif any([ss.monomer_eq(s) for ss in s0.species]): #this is used instead of in to deal with the potential for different parents
                s0_count = s0.monomer_count(s)#This is used instead of .count to deal with species with different parents
            else:
                s0_count = 0

            #Add the correct number of s to the list
            sf_count = sf.monomer_count(s)
            if s0_count < sf_count:
                species_to_add += (sf_count-s0_count)*[s]
            elif s0_count > sf_count:
                species_to_add = None #In this case, s0 contains more stuff than sf, so nothing should be returned
                break
            else:
                pass #if they have the same number, do not add it

        #s0 contains more or different species than sf, return None
        if (not isinstance(s0, ComplexSpecies)) and (sf.monomer_count(s0) == 0):
            species_to_add = None
        elif isinstance(s0, ComplexSpecies) and any([s0.monomer_count(s) > sf.monomer_count(s) for s in s0.species_set]) and sf.monomer_count(s0) == 0:
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

                    #append the new combination if it isn't excluded
                    if binder not in self.excluded_states and bindee not in self.excluded_states and cs not in self.excluded_states:
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
    
    def __init__(self, final_states, initial_states, intermediate_states = None, excluded_states = None, excluded_complexes = None, name = None, **keywords):
        """
        Binding reactions will be generated to form all PolymerConformations in final_states 
        from all the OrderedPolymerSpecies or PolymerConformations in initial_states.
        (or, if initial_states is None, from all the individual Polymers inside each PolymerConformation). 
        Intermediate states restricts the binding reactions to first form OrderdPolymerSpecies or PolymerConformations in this list. 
        At a high level this generates the following reactions:

        intial_states <-[Combinatorial Binding]-> final_states

        if  intermediate_states are given:
            intial_states <-[Combinatorial Binding]-> intermediate_states <-[Combinatorial Binding]->final_states
        
        
        Unlike CombinatorialComplex where Species are added individual, in CombinatorialConformation, groups of Species are added in single steps to produce the appropriate Complexes.

        :param final_states: a single PolymerConformation/OrderedPolymerSpecies or a list of said classes.
        :param initial_states: a list of initial OrderedPolymerSpecies which are bound together to form the final_states.
        :param intermediate_states: a list of intermediate PolymerConformation/OrderedPolymerSpecies formed when converting initial_states to final_states. 
                                    If None, all possible intermediate OrderedPolymerSpecies and PolymerConformations are enumerated.
        :param excluded_states: a list of intermediate PolymerConformation/OrderedPolymerSpecies which will not be formed during enumeration.
                                if None: no intermediates will be excluded.
        :param excluded_complexes: a list of ComplexSpecies which will not be formed inside or between any OrderedPolymerSpecies.
                                   if None: no intermediates will be excluded.
        """
        self.excluded_complexes = excluded_complexes
        super().__init__(final_states, initial_states = initial_states, intermediate_states = intermediate_states, excluded_states = excluded_states, name = name, **keywords)


    #Helper function to assert the correct class type
    def _assert_polymer_species_or_conformation(self, states, input_name = "states"):
        if not all([isinstance(s, PolymerConformation) or isinstance(s, OrderedPolymerSpecies) for s in states]):
            raise ValueError(f"{input_name} must be a list of PolymerConformation or OrderedPolymerSpecies. Recieved: {states}.")


    #Getters and setters
    #Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        final_states = list(self.set_species(final_states))

        #all final_states must be PolymerConformation or OrderedPolymerSpecies
        self._assert_polymer_species_or_conformation(final_states, "final_states")
        self._final_states = final_states

        #Then create a list of all sub-polymers included in final_states Conformations
        self.sub_polymers = []
        for s in self.final_states:
            if isinstance(s, PolymerConformation):
                for p in s.polymers:
                    self.sub_polymers.append(p)
            elif isinstance(s, OrderedPolymerSpecies):
                self.sub_polymers.append(s)

        self.sub_polymers = list(set(self.sub_polymers))

    #Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        #set initial states
        if initial_states is None:
            self._initial_states = self.sub_polymers
        else:
            initial_states = list(self.set_species(initial_states))
            #all initial_states must be PolymerConformation or OrderedPolymerSpecies
            self._assert_polymer_species_or_conformation(initial_states, "initial_states")

            #all initial states must have a valid route to at least one final state
            #TODO

            self._initial_states = initial_states

    #Intermediate states allows the user to restrict the complexes formed between the intial state and final state
    @property
    def intermediate_states(self):
        return self._intermediate_states
    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            intermediate_states = list(self.set_species(intermediate_states))

            #All intermediate_states must be OrderedPolymerSpecies or PolymerConformations
            self._assert_polymer_species_or_conformation(intermediate_states, "intermediate_states")

            #All intermediate_states must be constructable from at least one initial states
            #TODO


            self._intermediate_states = intermediate_states

    #excluded_states are OrderedPolymerSpecies or PolymerConformations which are not allowed to form
    @property
    def excluded_states(self):
        return self._excluded_states
    @excluded_states.setter
    def excluded_states(self, excluded_states):
        if excluded_states is None:
            self._excluded_states = []
        else:
            #All excluded states must be OrderedPolymerSpecies or PolymerConformations
            excluded_states = list(self.set_species(excluded_states))
            self._assert_polymer_species_or_conformation(excluded_states, "excluded_states")
            self._excluded_states = excluded_states

    #excluded_complexes are ComplexSpecies not allowed to form inside the OrderedPolymerSpecies or PolymerConformation
    @property
    def excluded_complexes(self):
        return self._excluded_complexes
    @excluded_complexes.setter
    def excluded_complexes(self, excluded_complexes):
        if excluded_complexes is None:
            self._excluded_complexes = []
        else:
            excluded_complexes = list(self.set_species(excluded_complexes))
            if not all([isinstance(s, ComplexSpecies) for s in excluded_complexes]):
                raise ValueError(f"excluded_complexes may only contain ComplexSpecies. Recievied: {excluded_complexes}.")
            self._excluded_complexes = excluded_complexes


    def compute_polymer_mapping(self, s0, sf):
        """
        Returns a mapping of the polymers in s0 to those in sf. 
        Also figures out which OrderedPolymerSpecies need to s0 be added to get a final PolymerConformation sf.

        WARNING: currently it is required that there be a unique mapping between polymers in s0 and polymers in sf. 
        If multiple polymers in s0 could become the same polymer in sf or one polymer in s0 can become multiple polymers in sf
        an error is thrown. A future version of the software will cover this more complex kind of enumeration.

        returns a tuple (polymer_mapping, polymers_to_add):
            polymer_mapping: dictionary {p0 : [list of new complexes], p0 : final polymers pf, pf : initial polymers p0} 
                             of initial polymers, the final polymer they become, and a list of complexes needed to create them
            polymers_to_add: a list of polymers in sf not in s0.
        """
        
        polymer_mapping = {} #
        polymers_to_add = [] #stores polymers not in s0 that are in sf

        #Case 1: s0 and sf are OrderedPolymerSpecies
        if isinstance(s0, OrderedPolymerSpecies) and isinstance(sf, OrderedPolymerSpecies):
            polymer_complexes_to_add = self.compute_complexes_to_add_to_polymer(s0, sf)
            if polymer_complexes_to_add is not None:
                polymer_mapping[s0] = sf
                polymer_mapping[sf] = s0
                polymer_mapping[s0, sf] = polymer_complexes_to_add

        #Case 2: s0 is an OrderedPolymerSpecies and sf is a PolymerConformation
        elif isinstance(s0, OrderedPolymerSpecies) and isinstance(sf, PolymerConformation):

            #Compare s0 to all other polymers in p
            for pf in sf.polymers:
                #p is different from s0 and cannot be created from s0
                polymer_complexes_to_add = self.compute_complexes_to_add_to_polymer(s0, pf)

                #If there is no mapping
                if polymer_complexes_to_add is None:
                    polymers_to_add.append(copy.copy(pf).remove()) #remove p from any PolymerConformation it is part of before returning.

                #p is the same as s0 or can be created from s0
                elif s0 not in polymer_mapping and pf not in polymer_mapping:
                    polymer_mapping[pf] = s0
                    polymer_mapping[s0] = pf #keep track of the polymers s0 can become
                    polymer_mapping[s0, pf] = polymer_complexes_to_add #keep track of a "recipe" to go from s0 to pf
                else:
                    #This means there are multiple copies of s0 or a polymer constructable from s0 in the conformation
                    #This will require an additional combinatorial search.
                    #TODO
                    raise NotImplementedError("Mapping between polymers is not one-to-one.")

        #Case 3: s0 and sf are PolymerConformations
        elif isinstance(s0, PolymerConformation) and isinstance(sf, PolymerConformation):
            #Compare each polymer in s0 to each polymer in sf
            for pf in sf.polymers:
                #pf is not equal to or constructable from any p0 in s0.polymers
                mapping_found = False
                for p0 in s0.polymers:
                    polymer_complexes_to_add = self.compute_complexes_to_add_to_polymer(p0, pf)

                    #Toggle that a mapping to pf has been found
                    if polymer_complexes_to_add is not None:
                        mapping_found = True

                        if p0 not in polymer_mapping and pf not in polymer_mapping:
                            polymer_mapping[pf] = p0
                            polymer_mapping[p0] = pf #keep track of the polymers s0 can become
                            polymer_mapping[p0, pf] = polymer_complexes_to_add #keep track of a "recipe" to go from s0 to pf
                        else:
                            #This means there are multiple pf in sf.polymers which can be constructed from the same p0 in s0.polymers
                            #This will require additional combinatorial search
                            #TODO
                            raise NotImplementedError("Mapping between polymers is not one-to-one.")

                #If there is no mapping, pf will need to be added as a seperate polymer
                if not mapping_found:
                    polymers_to_add.append(pf)
        else:
            raise ValueError(f"s0 and sf must be OrderedPolymerSpecies or PolymerConformations. Received {s0} and {sf}.")

        if len(polymer_mapping) == 0:
            return None #This indicates there is no mapping between s0 and sf
        else:
            return polymer_mapping, polymers_to_add


    def compute_complexes_to_add_to_polymer(self, p0, pf):
        """
        similar to CombinatorialComplex.compute_species_to_add 
        figures out which ComplexSpecies need to be added to a convert OrderedPolymerSpecies s0 to sf
        """
        if not isinstance(p0, OrderedPolymerSpecies) or not isinstance(pf, OrderedPolymerSpecies):
            raise ValueError(f"p0 and pf must be OrderedPolymerSpecies; recieved p0 = {p0}, pf = {pf}")

        #Polymers must be the same length
        if len(p0) != len(pf):
            return None
        #Nothing to add because they are teh same
        elif p0.monomer_eq(pf):
            return []
        #Compare monomer by monomer
        else:
            complexes_to_add = [] #will be a list of tuples (index, [Species to add])
            for i in range(len(p0)):
                mi = p0[i]
                mf = pf[i]

                
                if isinstance(mf, ComplexSpecies):
                    species_to_add = self.compute_species_to_add(mi, mf)

                    if species_to_add is not None and len(species_to_add)>0:
                        complexes_to_add.append((i, species_to_add))
                    elif species_to_add is None: #mi cannot become mf, no mapping exists
                        complexes_to_add = None 
                        break
                elif not mi.monomer_eq(mf): #mi cannot become mf, no mapping exists
                    complexes_to_add = None 
                    break

                print(mi, mf, complexes_to_add)

            return complexes_to_add

    def compute_complexes_to_add_to_conformation(self, s0, sf, polymer_mapping = None):
        """
        similar to CombinatorialComplex.compute_species_to_add 
        figures out which ComplexSpecies need to be added to a convert PolymerConformation/OrderedPolymerSpecies s0 to PolymerConformation sf
        polymer_mapping is a mapping between polymers in sf and polymers in s0. Computed via compute_polymer_mapping.

        returns a list of lists of tuples.
        [   
            lists of complexes in the conformation in the form of:
            [list of elements in a complex: (p0, species index in p0)]
        ]
        """

        if polymer_mapping is None:
            mapping = self.compute_polymer_mapping(s0, sf)
            if mapping is not None:
                polymer_mapping, polymers_to_add = mapping
            else: #If there is no mapping, return None
                return None

        #A list of the polymers in s0
        if isinstance(s0, OrderedPolymerSpecies):
            current_polymers = [s0]
            current_complexes = []
        elif isinstance(s0, PolymerConformation):
            current_polymers = s0.polymers
            current_complexes = s0.complexes
        else:
            raise ValueError(f"s0 must be an OrderedPolymerSpecies or a PolymerConformation; recievied {s0}.")

        #a list of polymers in sf
        if isinstance(sf, OrderedPolymerSpecies):
            return None #If sf is an OrderedPolymerSpecies, there are no Complexes to add
        elif isinstance(sf, PolymerConformation):
            final_polymers = sf.polymers
            final_complexes = sf.complexes
        else:
            raise ValueError(f"sf must be an OrderedPolymerSpecies or a PolymerConformation; recievied {sf}.")


        complexes_to_add = []
        complex_mapping = {}

        #cycle through final complexes
        for cf in final_complexes:
            #Compute complex position
            cf_position = [] #[(polymer, index)....]
            cf_additional_species = [] #[Species]
            pf_inds = sf.get_polymer_indices(cf) #computes the indices polymer index (or None) of each species in the complex
            pf_positions = sf.get_polymer_positions(cf)

            for i in range(len(cf.species)):
                if pf_inds[i] is not None:
                    if sf.polymers[pf_inds[i]]  in polymer_mapping:
                        p = polymer_mapping[sf.polymers[pf_inds[i]]] #the initial polymer mapped to the final polymer if possible
                    else:
                        p = sf.polymers[pf_inds[i]] #otherwise use the final polymer (these means that polymer is in the polymers_to_add list)

                    cf_position.append((p, pf_positions[i]))

                else: #this means the species isn't part of a Polymer
                    cf_additional_species.append(cf.species[i])

            #cycle through current_complexes - if a current complex can become a final complex
            #the species need to be added to the current_complex
            for cc in current_complexes:
                #Compute complex position
                cc_position = [] #[(polymer, index)....]
                cc_additional_species = [] #[Species]
                pc_inds = s0.get_polymer_indices(cc) #computes the indices polymer index (or None) of each species in the complex
                cc_positions = s0.get_polymer_positions(cc)

                for i in range(len(cc.species)):
                    if pc_inds[i] is not None:
                        cc_position.append((s0.polymers[pc_inds[i]], cc_positions[i]))

                    else: #this means the species isn't part of a Polymer
                        cc_additional_species.append(cc.species[i])

                if all([p in cf_position for p in cc_position]) and all([cf.monomer_count(s) >= cc.monomer_count(s) for s in cc_additional_species]):
                    #This means cc can become cf

                    if cc in complex_mapping:
                        #This means that there are multiple Complexes inside the conformation which can mutually convert
                        #This more complex enumeration is not being implemented at this time
                        raise NotImplementedError("Mapping between Conformations not one-to-one.")
                    else:
                        complex_mapping[cc] = cf
                        
                        for p in cc_position:
                            cf_position.remove(p)
                        for s in cc_additional_species:
                            cf_additional_species.remove(s)

                        #add the current complex to the list if there are things to add to it
                        if len(cf_position) > 0 or len(cf_additional_species) > 0:
                            cf_position = [(s0, cc)] + cf_position

            #Only add the Complex if it contains something new
            if len(cf_position) > 0:
                complexes_to_add.append((cf_position, cf_additional_species))

        #Ensure every complex in s0 maps to a complex in sf
        if len(current_complexes) > 0 and any([cc not in complex_mapping for cc in current_complexes]):
            return None
        else:
            return complex_mapping, complexes_to_add


    def get_combinations_between(self, s0, sf):
        """
        Returns a list of lists. Each sublist contains tuples of the form:

            ([(p0, index)...], [additional species]): to represent the formation of a complex between the monomers in the polymers p0 at the respective indices which also includes the additional species
            ([(s0, complex)...], [additional species]): to represent the formation of a complex between the complex in the confomrations s0 which also includes the additional species

        Note: if the length of the first list is 1, it means that the Complex is formed inside an individual Polymer, not between Polymers.
        """

        #A list of the polymers in s0
        if isinstance(s0, OrderedPolymerSpecies):
            current_polymers = [s0]
            current_complexes = []
        elif isinstance(s0, PolymerConformation):
            current_polymers = s0.polymers
            current_complexes = s0.complexes
        else:
            raise ValueError(f"s0 must be an OrderedPolymerSpecies or a PolymerConformation; recievied {s0}.")

        #return nothing if they are the same
        if s0 == sf:
            return None

        #Get the mapping between polymers in s0 and sf
        polymer_mapping = self.compute_polymer_mapping(s0, sf)
        if polymer_mapping is None:
            print("polymer mapping is None!")
            return None
        else:
            polymer_mapping, polymers_to_add = polymer_mapping

        #convert the mapping to a list
        complexes_to_add_to_polymers = []
        for p in current_polymers:
            pf = polymer_mapping[p]
            if (p, pf) in polymer_mapping and len(polymer_mapping[(p, pf)])>0:
                complexes_to_add_to_polymers += [([(p, index)], species_list) for (index, species_list) in polymer_mapping[(p, pf)]]

        print("polymers_to_add", polymers_to_add)
        print("complexes_to_add_to_polymers", complexes_to_add_to_polymers)
        #get the complexes to add to the conformation
        complex_mapping = self.compute_complexes_to_add_to_conformation(s0, sf, polymer_mapping)
        

        #Create a list of all complexes to add
        #the elements of this list will imply if they are complexes between polymers or complexes in polymers based upon their form
        if complex_mapping is None:
            print("complexes_to_add_to_conformation is None")
            all_complexes_to_add = complexes_to_add_to_polymers
        else:
            complex_mapping, complexes_to_add_to_conformation = complex_mapping
            all_complexes_to_add = complexes_to_add_to_conformation+complexes_to_add_to_polymers
            print("complexes_to_add_to_conformation", complexes_to_add_to_conformation)

        
        print("all_complexes_to_add", all_complexes_to_add)
        combos = [list(i) for i in permutations(all_complexes_to_add, len(all_complexes_to_add))]
        print("combos", combos)

        return combos
        




                





    
    


