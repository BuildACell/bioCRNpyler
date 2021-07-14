#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from .species import Complex, ComplexSpecies, OrderedPolymerSpecies, PolymerConformation
from .component import Component
from .components_combinatorial_complex import CombinatorialComplex
from . mechanisms_conformation import One_Step_Reversible_Conformation_Change
from .dna_part_promoter import Promoter
from itertools import permutations
import copy
import warnings

class CombinatorialConformation(Component):
    """
    A class to represent a PolymerConformation (made of one unique OrderedPolymerSpecies) with many internal Complexes which can bind and unbind in many different ways.
    """
    
    def __init__(self, final_states, initial_states = None, intermediate_states = None, excluded_states = None, state_part_ids = None, name = None, **keywords):
        """
        Binding reactions will be generated to form all PolymerConformations in final_states 
        from all the PolymerConformations in initial_states. There must be a single, unique, OrderedPolymerSpecies in all the conformations.
        Intermediate states restricts the binding reactions to first form PolymerConformations in this list. 
        At a high level this generates the following reactions:

        intial_states <-[Combinatorial Binding]-> final_states

        if  intermediate_states are given:
            intial_states <-[Combinatorial Binding]-> intermediate_states <-[Combinatorial Binding]->final_states
        
        
        Unlike CombinatorialComplex where Species are added individual, in CombinatorialConformation, groups of Species are added in single steps to produce the appropriate Complexes.

        :param final_states: one or more PolymerConformations.
        :param initial_states: a list of initial PolymerConformations which can bind/unbind to become the final_state
        :param intermediate_states: a list of intermediate PolymerConformations formed when converting initial_states to final_states. 
                                    If None, all possible intermediate PolymerConformations are enumerated.
        :param excluded_states: a list of intermediate PolymerConformations which will not be formed during enumeration.
                                if None: no intermediates will be excluded.
        :param state_part_ids: a dictionary {PolymerConformation : str} used to generate shorter part-ids for this conformation
        """

        if state_part_ids is None:
            self.state_part_ids = {}
        else:
            self.state_part_ids = state_part_ids
        self.internal_polymer = None #this will be set inside final_states setter
        self.final_states = final_states
        self.initial_states = initial_states
        self.intermediate_states = intermediate_states
        self.excluded_states = excluded_states

        if name is None:
            name = str(self.internal_polymer.name)
        Component.__init__(self, name = name, **keywords)
        self.add_mechanism(One_Step_Reversible_Conformation_Change())


    #Helper function to assert the correct class type
    def _assert_conformation(self, states, input_name = "states"):
        if not all([isinstance(s, PolymerConformation) for s in states]):
            raise ValueError(f"{input_name} must be a list of PolymerConformations. Recieved: {states}.")
        if not all([len(s.polymers) == 1 for s in states]):
            raise ValueError(f"All PolymerConformations in {input_name} must contain a single unique internal OrderedPolymerSpecies. Recieved: {states}.")

        if self.internal_polymer is None:
            self.internal_polymer = copy.deepcopy(states[0].polymers[0])
            self.internal_polymer.parent = None

        if not all([str(p) == str(self.internal_polymer) for s in states for p in s.polymers]):
            raise ValueError(f"All PolymerConformations in {input_name} must contain a single unique internal OrderedPolymerSpecies. Recieved: {states}.")

    #Getters and setters
    #Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        final_states = list(self.set_species(final_states))

        self._assert_conformation(final_states, "final_states")
        self._final_states = final_states

    #Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        #set initial states
        if initial_states is None or len(initial_states) == 0:
            self._initial_states = [PolymerConformation(polymer = self.internal_polymer)]
        else:
            initial_states = list(self.set_species(initial_states))
            #all initial_states must be PolymerConformation
            self._assert_conformation(initial_states, "initial_states")

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

            #All intermediate_states must be PolymerConformations
            self._assert_conformation(intermediate_states, "intermediate_states")

            self._intermediate_states = intermediate_states

    #excluded_states are PolymerConformations which are not allowed to form
    @property
    def excluded_states(self):
        return self._excluded_states
    @excluded_states.setter
    def excluded_states(self, excluded_states):
        if excluded_states is None:
            self._excluded_states = []
        else:
            #All excluded states must be PolymerConformations
            excluded_states = list(self.set_species(excluded_states))
            self._assert_conformation(excluded_states, "excluded_states")
            self._excluded_states = excluded_states

    def compute_species_changes(self, s0, sf):
        #print("computing species changes between", s0, "and", sf)
        #Compute Species that need to be added to s0 to get the PolymerConformation sf
        #Assumes the underlying internal polymer is the same
        #Computes a list of moves to go from s0 --> sf. Each move produces one of the Complexes in SF


        #1. if c0 contains bound locations not in cf, c0 cannot be transformed additively into cf
        if any([len(s0.get_complexes_at(0, i)) > len(sf.get_complexes_at(0, i)) for i in range(len(self.internal_polymer))]):
            species_changes = False
        else:
            species_changes = {}
            merged_complexes = {}

            #Figure out what happens to create each complex in sf
            for cf in sf.complexes:
                pf_inds = sf.get_polymer_positions(cf, 0)

                #initially, assume all external Species should be added
                #species are stored as (s, index) in order to differentiate species in different parts of the polymer
                cf_species_and_inds = [(cf.species[i], pf_inds[i])  for i in range(len(cf.species))]
                cf_species_and_inds_str = [str(i) for i in cf_species_and_inds] #string version used to check for equality
                sub_complex_species_and_inds = []
                sub_complex_species_and_inds_str = []
                merged_complexes[cf, pf_inds] = []


                #Iterate through speces in the polymer sf
                for ind in pf_inds:

                    #get complex of the polymer at this location
                    found_complexes = s0.get_complexes_at(0, ind)
                    for c0 in found_complexes:
                        p0_inds = s0.get_polymer_positions(c0, 0)
                        #Ensure each subcomplex c0 is only evaluated once
                        #then add all the Species in the complex to the sub_complex_species_and_inds list
                        if (cf, pf_inds) not in merged_complexes or c0 not in merged_complexes[cf, pf_inds]:
                            merged_complexes[cf, pf_inds].append(c0)

                            sub_complex_species_and_inds += [(c0.species[i], p0_inds[i])  for i in range(len(c0.species))]
                            sub_complex_species_and_inds_str += [str(i)  for i in sub_complex_species_and_inds]


                #cf cannot be created additively from the complexes in s0
                sub_complex_exclusive_species = [i for i in sub_complex_species_and_inds_str if i not in cf_species_and_inds_str]
                cf_complex_exclusive_species = [i for i in cf_species_and_inds_str if i not in sub_complex_species_and_inds_str]
                if len(sub_complex_exclusive_species) > 0:
                    species_changes[cf, pf_inds] = False
                #keep track of all the external species added to produce the new complex
                elif len(cf_complex_exclusive_species) > 0:
                    species_changes[cf, pf_inds] = [cf_species_and_inds[i][0] for i in range(len(cf_species_and_inds_str)) if (cf_species_and_inds_str[i] not in sub_complex_species_and_inds_str and cf_species_and_inds[i][1] is None)]
        #If any complex cannot be created, return False
        if species_changes == False or any([species_changes[k] is False for k in species_changes]) or (len(species_changes) == 0 and len(merged_complexes) == 0):
            #print("cannot convert.")
            return False
        else:
            #print("species_changes =", species_changes)
            #print("merged_complexes =", merged_complexes)
            return species_changes, merged_complexes



    def get_combinations_between(self, s0, sf):
        """
        Returns a list of 
        """
        #print("geting combinations between", s0, "and", sf)
        X = self.compute_species_changes(s0, sf)


        if X is False:
            return []
        else:
            species_changes, merged_complexes = X
            changes_list = []
            for cf in sf.complexes:
                inds = sf.get_polymer_positions(cf, 0)
                #print("checking ", cf, "inds =", inds, "in species_changes", (cf, inds) in species_changes, "in merged_complexes", (cf, inds) in merged_complexes)

                if ((cf, inds) in merged_complexes and len(merged_complexes[(cf, inds)]) == 1 and merged_complexes[(cf, inds)][0].monomer_eq(cf)) or (cf, inds) not in merged_complexes:
                    complexes_to_merge = None
                else:
                    complexes_to_merge = merged_complexes.get((cf, inds))

                if (cf, inds) not in species_changes or len(species_changes[cf, inds]) == 0:
                    species_to_add = None
                else:
                    species_to_add = species_changes.get((cf, inds))

                #print("(cf, inds, species_to_add, complexes_to_merge)=", (cf, inds, species_to_add, complexes_to_merge))
                if species_to_add is not None or complexes_to_merge is not None:
                    changes_list.append((cf, inds, species_to_add, complexes_to_merge))


            #get all permutations of the species for different binding orders 
            perms = permutations(changes_list, len(changes_list))

            #combinations (binder, bindee, complex_species) to be returned
            combinations = []

            #Iterate through permutations
            for perm in perms:
                old_state = s0
                for cf, inds, species_to_add, complexes_to_merge in perm:

                    #Determine all the species that are bound together
                    species_list = []
                    for i, p_loc in enumerate(inds):
                        if p_loc is None:
                            species_list.append(cf.species[i])
                        else:
                            species_list.append(old_state.polymers[0][p_loc])

                    new_complex = Complex(species_list)

                    #Determine which complexes to remove from the conformation
                    if complexes_to_merge is None:
                        merged_complexes = None
                    else:
                        merged_complexes = [new_complex.parent.get_complex(c) for c in complexes_to_merge]

                    #create a neew polymer conformation by removing the merged complexes
                    if merged_complexes not in [[], None]:
                        new_state = new_complex.parent.copy_remove_complexes(merged_complexes)
                    else:
                        new_state = new_complex.parent

                    #append the new combination if it isn't excluded
                    if old_state not in self.excluded_states and new_state not in self.excluded_states:
                        combinations.append((old_state, species_to_add, new_state))

                    #update bindee
                    old_state = new_state
            #print("combinations = ", combinations)
            return combinations

    def _get_part_id(self, state):
        if state in self.state_part_ids:
            return self.state_part_ids[state]
        else:
            return str(state)

    def update_species(self):
        mech_c = self.get_mechanism("conformation_change")
        species = []
        species_added_dict = {} #save which combinations have already been added
        self.combination_dict = {} #this should be recomputed every updated species

        #If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    if si != s0:
                        #Get combinatorial species between s0 and si
                        if (s0, si) not in self.combination_dict:
                            self.combination_dict[s0, si] = self.get_combinations_between(s0, si)

                        #iterate through combinations of species between s0 and si
                        for (s00, additional_species, sff) in self.combination_dict[s0, si]:
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            species += mech_c.update_species(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)

            for si in self.intermediate_states:
                for sf in self.final_states:
                    if si != sf:
                        #Get combinatorial species between si and sf
                        if (si, sf) not in self.combination_dict:
                            self.combination_dict[si, sf] = self.get_combinations_between(si, sf)

                        #iterate through combinations of species between si and sf
                        for (s00, additional_species, sff) in self.combination_dict[si, sf]:
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            species += mech_c.update_species(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)
            
        #If there are no intermediate restrictions, compute combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    if s0 != sf:
                        #Get combinatorial species between s0 and sf
                        if (s0, sf) not in self.combination_dict:
                            self.combination_dict[s0, sf] = self.get_combinations_between(s0, sf)
                            
                        #iterate through combinations of species between s0 and sf
                        for (s00, additional_species, sff) in self.combination_dict[s0, sf]:
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            species += mech_c.update_species(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)

        return list(set(species))

    def update_reactions(self):
        mech_c = self.get_mechanism("conformation_change")
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
                    for (s00, additional_species, sff) in self.combination_dict[s0, si]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            reactions += mech_c.update_reactions(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)

            for si in self.intermediate_states:
                for sf in self.final_states:
                    #Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = self.get_combinations_between(si, sf)

                    #iterate through combinations of species between si and sf
                    for (s00, additional_species, sff) in self.combination_dict[si, sf]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            reactions += mech_c.update_reactions(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)

        #If there are no intermediate restrictions, compute combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    #Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = self.get_combinations_between(s0, sf)
                        
                    #iterate through combinations of species between s0 and sf
                    for (s00, additional_species, sff) in self.combination_dict[s0, sf]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = f"{self._get_part_id(s00)}-{self._get_part_id(sff)}"
                            reactions += mech_c.update_reactions(s0 = s00, sf = sff, additional_species = additional_species, 
                            component = self, part_id = part_id)

        return reactions



class CombinatorialConformationPromoter(CombinatorialConformation, Promoter):
    """
    A combinatorial conformaiton with an additional set of states "expressing_states" which can transcribe/express rna/protein products.
    This class merges Promoter and CombinatorialConformation.

        :param promoter_states: one or more PolymerConformations which are used by the promoter class
        :param promoter_states_on: True/False if True all promoter_states are transcribable. If False all states except promoter_states are transcribable
        :param promoter_location: the index of the monomer in the PolymerConformation which represents the promoter
        :param final_states: one or more PolymerConformations.
        :param initial_states: a list of initial PolymerConformations which can bind/unbind to become the final_state
        :param intermediate_states: a list of intermediate PolymerConformations formed when converting initial_states to final_states. 
                                    If None, all possible intermediate PolymerConformations are enumerated.
        :param excluded_states: a list of intermediate PolymerConformations which will not be formed during enumeration.
                                if None: no intermediates will be excluded.
        :param state_part_ids: a dictionary {PolymerConformation : str} used to generate shorter part-ids for this conformation
        :param activating_complexes: a list of ComplexSpecies which activate PolymerConformations allowing them to be transcribed.
        :param inactivating_complexes: a list of ComplexSpecies which innactive the PolymerConformation, preventing transcription.

    """
    def __init__(self, promoter_states, promoter_location, promoter_states_on = True, activating_complexes = None, inactivating_complexes = None,
        name = "CombinatorialConformationPromoter", **keywords):
        
        Promoter.__init__(self, name, **keywords)
        CombinatorialConformation.__init__(self, name = name, **keywords)

        self.promoter_states = promoter_states
        self.promoter_states_on = promoter_states_on

        if activating_complexes is None:
            self.activating_complexes = []
        else:
            self.activating_complexes = activating_complexes

        if inactivating_complexes is None:
            self.inactivating_complexes = []
        else:
            self.inactivating_complexes = inactivating_complexes

        assert all([isinstance(c, ComplexSpecies) for c in self.activating_complexes])
        assert all([isinstance(c, ComplexSpecies) for c in self.inactivating_complexes])

        if promoter_location not in range(len(self.internal_polymer)):
            raise ValueError(f"promoter_location must be an index of the polymer {self.internal_polymer}. Recieved {promoter_location}.")
        else:
            self.promoter_location = promoter_location


    #promoter_states are PolymerConformations which are used by the Promoter class
    #These can be ON or OFF
    @property
    def promoter_states(self):
        return self._promoter_states

    @promoter_states.setter
    def promoter_states(self, promoter_states):
        if promoter_states is None:
            self._promoter_states = []
        else:
            #All excluded states must be PolymerConformations
            promoter_states = list(self.set_species(promoter_states))
            self._assert_conformation(promoter_states, "promoter_states")
            self._promoter_states = promoter_states


    def update_species(self):
        self.conformation_species = CombinatorialConformation.update_species(self)
        promoter_species = []
        old_name = self.name
        for s in self.conformation_species:
            if isinstance(s, PolymerConformation):
                active_state = False
                if  (s in self.promoter_states and self.promoter_states_on) or (s not in self.promoter_states and not self.promoter_states_on):
                    active_state = True

                active_complex = False
                if any([s.get_complex(c) is not None for c in self.activating_complexes]):
                    active_complex = True

                innactive_complex = False
                if any([s.get_complex(c) is not None for c in self.inactivating_complexes]):
                    innactive_complex = True
                    if (active_state and len(self.promoter_states)>0) or active_complex:
                        warnings.warn("innactive_complex conflicts with active_complex or active_state in CombinatorialConformationPromoter. Defaulting to innactive.")
                        
                if (active_state or active_complex) and not innactive_complex:
                    self.dna_to_bind = s.polymers[0][self.promoter_location]
                    self.name = self._get_part_id(s) #Reset name for unique addressable promoter states
                    promoter_species += Promoter.update_species(self)
        
        return list(set(promoter_species+self.conformation_species))

    def update_reactions(self):
        if not hasattr(self, "conformation_species"):
            self.update_species

        conformation_reactions = CombinatorialConformation.update_reactions(self)
        promoter_reactions = []
        old_name = self.name
        for s in self.conformation_species:
            if isinstance(s, PolymerConformation):
                active_state = False
                if  (s in self.promoter_states and self.promoter_states_on) or (s not in self.promoter_states and not self.promoter_states_on):
                    active_state = True

                active_complex = False
                if any([s.get_complex(c) is not None for c in self.activating_complexes]):
                    active_complex = True

                innactive_complex = False
                if any([s.get_complex(c) is not None for c in self.inactivating_complexes]):
                    innactive_complex = True
                    if (active_state and len(self.promoter_states)>0) or active_complex:
                        warnings.warn("innactive_complex conflicts with active_complex or active_state in CombinatorialConformationPromoter. Defaulting to innactive.")
                        
                if (active_state or active_complex) and not innactive_complex:
                    self.dna_to_bind = s.polymers[0][self.promoter_location]
                    self.name = self._get_part_id(s) #Reset name for unique addressable promoter states
                    promoter_reactions += Promoter.update_reactions(self)
        self.name = old_name
        return promoter_reactions + conformation_reactions