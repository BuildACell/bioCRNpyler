################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 12/21/2020
#
################################################################


# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from .component import Component
from .components_basic import DNA, RNA
from .dna_part import DNA_part
from .dna_part_misc import IntegraseSite
from .species import ComplexSpecies, OrderedPolymer, OrderedPolymerSpecies
from .utils import all_comb, remove_bindloc
import logging

class Construct(Component,OrderedPolymer):
    def __init__(self,
                parts_list,
                name=None,
                circular=False,
                mechanisms=None,  # custom mechanisms
                parameters=None,  # customized parameters
                attributes=None,
                initial_concentration=None, 
                component_enumerators = None,
                **keywords):
        """this represents a bunch of parts in a row.
        A parts list has [[part,direction],[part,direction],...]
        Each part must be an OrderedMonomer"""

        if(component_enumerators is None):
            component_enumerators = []
        self.component_enumerators = component_enumerators
        OrderedPolymer.__init__(self,parts_list, default_direction="forward")
        self.circular=circular
        if(name is None):
            name = self.make_name() #automatic naming
        self.name = name
        Component.__init__(self=self,name=name,length = len(parts_list),
                    mechanisms=mechanisms,parameters=parameters,
                    attributes=attributes,initial_concentration = initial_concentration,
                     **keywords)
        self.update_parameters()
        self.transcripts = []
        if(not hasattr(self,"material_type")):
            self.material_type=None #set this when you inherit this class
        self.out_components = None
        self.predicted_rnas = None
        self.predicted_proteins = None

    @property
    def parts_list(self):
        return self.polymer
    def make_name(self):
        output = ""
        outlst = []
        for part in self.parts_list:
            pname = part.name
            if(part.direction=="reverse"):
                pname+="_r"
            outlst += [pname]
        output = '_'.join(outlst)
        if(self.circular):
            output+="_o"
        return output

    def get_part(self,part = None, part_type=None, name = None, index = None):
        """
        Function to get parts from Construct.parts_list.

        One of the 3 keywords must not be None.

        part: an instance of a DNA_part. Searches Construct.parts_list for a DNA_part with the same type and name.
        part_type: a class of DNA_part. For example, Promoter. Searches Construct.parts_list for a DNA_part with the same type.
        name: str. Searches Construct.parts_list for a DNA_part with the same name
        index: int. returns Construct.parts_list[index]

        if nothing is found, returns None.
        """

        if [part, name, index,part_type].count(None) != 3:
            raise ValueError(f"get_component requires a single keyword. Recieved component={part}, name={name}, index={index}.")
        if not (isinstance(part, DNA_part) or part is None):
            raise ValueError(f"component must be of type DNA_part. Recieved {part}.")
        if not (type(part_type) == type or part_type is None):
            raise ValueError(f"part_type must be a type. Recieved {part_type}.")
        if not (isinstance(name, str) or name is None):
            raise ValueError(f"name must be of type str. Recieved {name}.")
        if not (isinstance(index, int) or index is None):
            raise ValueError(f"index must be of type int. Recieved {index}.")

        matches = []
        if index is not None:
            matches.append(self.parts_list[index])
        else:
            for comp in self.parts_list:
                if part is not None:
                    if type(part) == type(comp) and comp.name == part.name:
                        matches.append(comp)
                elif name is not None:
                    if comp.name == name:
                        matches.append(comp)
                elif part_type is not None:
                    if(isinstance(comp,part_type)):
                        matches.append(comp)
        if len(matches) == 0:
            return None
        elif len(matches) == 1:
            return matches[0]
        else:
            warn("get_part found multiple matching components. A list has been returned.")
            return matches 


    def reverse(self):
        """reverses everything, without actually changing the DNA.
        also updates the name and stuff, since this is now a different Construct"""
        OrderedPolymer.reverse(self)
        self.reset_stored_data()
        self.name = self.make_name()
        #self.update_base_species()
        return self
    def get_reversed(self):
        """returns a reversed version of this construct without changing this construct"""
        outcon = copy.deepcopy(self)
        outcon.reverse()
        return outcon
    def get_circularly_permuted(self,new_first_position):
        """returns a new construct which has the first position changed to new_first_position"""
        if(not self.circular):
            return ValueError("cannot circularly permute linear construct")
        else:
            return DNA_construct(self.parts_list[new_first_position:]+self.parts_list[:new_first_position], circular=True, 
                            component_enumerators = self.component_enumerators, attributes=self.attributes,mechanisms=self.mechanisms, mixture=self.mixture)
    def set_mixture(self, mixture):
        self.mixture = mixture
        for part in self.parts_list:
            part.set_mixture(mixture)

    def update_base_species(self, base_name=None, attributes = None):
        if base_name is None:
            self.base_species = self.set_species(self.name, material_type = self.material_type, attributes = attributes)
        else:
            self.base_species = self.set_species(base_name, material_type = self.material_type, attributes = attributes)

    def update_parameters(self, overwrite_parameters = True):
        """update parameters of all parts in the construct"""
        Component.update_parameters(self = self,parameter_database=self.parameter_database)
        for part in self.parts_list:
            part.update_parameters(parameter_database = self.parameter_database,
                                    overwrite_parameters = overwrite_parameters)

    def add_mechanism(self, mechanism, mech_type = None, overwrite = False, optional_mechanism = False):
        Component.add_mechanism(self, mechanism, mech_type = mech_type, \
                                overwrite = overwrite, optional_mechanism = optional_mechanism)
        for part in self.parts_list:
            part.add_mechanism( mechanism, mech_type = mech_type, \
                                overwrite = overwrite, optional_mechanism = optional_mechanism)
    
    def __repr__(self):
        """this is just for display purposes"""
        return "Construct = "+ self.make_name()

    def __contains__(self,obj2):
        """checks if this construct contains a certain part, or a copy of a certain part"""
        if(isinstance(obj2,DNA_part)):
            #if we got a DNA part it could mean one of two things:
            #1 we want to know if a dna part is anywhere
            #2 we want to know if a specific DNA part is in here
            #this is complicated by the fact that we want to have the same DNA part be reusable
            #in many locations
            if(obj2.parent==self):
                #the object should already know if it's a part of me
                return True
            elif(obj2.parent==None):
                #this object has been orphaned. 
                #that means we are looking for matching objects in any position
                new_obj2 = copy.copy(obj2).unclone()
                uncloned_list = [copy.copy(a).unclone() for a in self.parts_list]
                return new_obj2 in uncloned_list
            else:
                return False
        elif(isinstance(obj2,str)):
            #if we get a string, that means we want to know if the name exists anywhere
            return obj2 in str(self)

    def get_species(self):
        ocomplx = []
        for part in self.parts_list:
            partspec = copy.copy(part.dna_species)
            #partspec.material_type = self.material_type+"_part"
            ocomplx += [partspec.set_dir(part.direction)]
        out_species = OrderedPolymerSpecies(ocomplx,circular = self.circular,material_type=self.material_type)
        
        return out_species
    
    def located_allcomb(self,spec_list):
        """recursively trace all paths through a list
        [[[part1,1],[part2,5]],[[part3,1]],[[part4,5],[part5,12]]]
        ====================>
        compacted_indexes = [1,5,12]
        prototype_list = [[part1,part3],[part2,part4],[part5]]
        comb_list = [[1],[5],[12],[1,5],[1,12],[5,12],[1,5,12]]
        ===========================
        then, take the lists from comb_list and create all possible lists
        out of prototype_list that includes those elements"""
        #first we have to construct the list we are tracing paths through
        
        spec_indexes = [a.position for a in spec_list] #extract all indexes
        #the following takes apart the lists because i don't yet know how to deal
        #with multiple binders at the same time
        compacted_indexes = sorted(list(set(spec_indexes)))

        prototype_list = [None]*len(compacted_indexes)
        for spec in spec_list:
            #now, spec is a list which contains all the elements which are defined for each variant.
            #
            #go through every element and put it in the right place
            proto_ind = compacted_indexes.index(spec.position) #where to put it?
            if(prototype_list[proto_ind] is None):
                #if nothing's been placed here, then create a list
                prototype_list[proto_ind] = [spec]
            else:
                #if something is already here, then add to the list
                prototype_list[proto_ind]+= [spec]
        # at this point we have a list that looks like this:
        # [[[part1,0],[part2,0]],[[part2,3]],[[part3,12],[part5,12]]
        # next step is to pick one of the first list (either [part1,0] or [part2,0])
        # one of the second list (only [part2,3] is our option), one of the third list, etc
        # for all possible choices made this way
        comb_list = all_comb(compacted_indexes)
        def recursive_path(in_list):
            """this function takes every possible "path" through a list of lists.
            For example:
            input:  [[A,B],[C],[D,E]]
            we have two options for the first position, one option for the second, and two options for the third.
            so, output:[[A,C,D],[B,C,D],[A,C,E],[B,C,E]]
            """
            if(len(in_list)==1):
                out_list = []
                for a in in_list[0]:
                    out_list+= [[a]]
                return out_list
            elif(len(in_list)==0):
                return []
            else:
                out_list = []
                for a in in_list[0]:
                    out_list += [[a] + z for z in recursive_path(in_list[1:])]
                return out_list
        outlist = []
        for combo in comb_list:
            combo_sublists = []
            for combo_index in combo:
                combo_sublists += [prototype_list[compacted_indexes.index(combo_index)]]
            outlist+= recursive_path(combo_sublists)
        outdict_list = []
        for combo in outlist:
            replacedict = {}
            for spec in combo:
                replacedict.update({spec.position:[spec,spec.direction]})
            outdict_list += [replacedict]
        return outdict_list

    def make_polymers(self,species_lists,backbone):
        """makes polymers from lists of species
        inputs:
        species_lists: list of species which are to be assembled into a polymer
        backbone: the base_species which all these polymers should have"""
        polymers = []
        for combo in species_lists:
            #to make combinatorially bound versions of a species, we take the unbound species (backbone)
            #and replace monomers of it with bound versions (combo)
            #combo is a dictionary of the form: {number:OrderedMonomer,...}
            #where the number is which element should be replaced and the OrderedMonomer is what will be replacing that element
            new_backbone = OrderedPolymerSpecies.from_polymer_species(backbone,combo)
            polymers += [new_backbone] #we make a new OrderedComplexSpecies
        return polymers

    def update_combinatorial_complexes(self,active_components):
        """given an input list of components, we produce all complexes
        yielded by those components, mixed and matched to make all possible combinatorial
        complexes, where each component is assumed to only care about binding to one spot.
        First, the components are asked what species they make, then these species are sifted to
        reveal only the ones which are versions of the same polymer, just with different locations bound.
        Then, combinatorial combinations are made.
        for example:
            construct: <A,B,C>
        two new species are possible: <[A:RNAP],B,C>; <A,[B:RNAP],C>
        combinatorial species is also possible (since A and B are assumed to act independantly)
                                        <[A:RNAP],[B:RNAP],C>
        complexes, where each component is assumed to only care about binding to one spot"""
        species = []
        for part in active_components:
            #first we make binary complexes
            sp_list =  part.update_species()
            species+=sp_list
        unique_complexes = []
        #species need to be uniqueified
        unique_species = list(set(species)) 
        for specie in unique_species:
            #in this list we extract all the variants of the complexes from possible_backbones
            #that exist in our species list.
            if(isinstance(specie,ComplexSpecies) and specie.position is not None):
                unique_complexes+=[specie]
        #unique_complexes now has a list of all the non-combinatorial complexes we can make
        combinatorial_complexes = unique_complexes
        allcomb = self.located_allcomb(combinatorial_complexes) #all possible combinations of binders are made here
        allcomb += [{}] #unbound dna should also be used
        #now, all possibilities have been enumerated.
        #we construct the OrderedPolymerSpecies
        out_polymers = self.make_polymers(allcomb,self.get_species())
        return out_polymers

    #Overwrite Component.enumerate_components 
    def enumerate_constructs(self,previously_enumerated=None):
        """Runs all our component enumerators to generate new constructs"""
        new_constructs = []
        for enumerator in self.component_enumerators:
            new_comp = enumerator.enumerate_components(component=self,previously_enumerated=previously_enumerated)
            new_constructs += new_comp
        return new_constructs

    def combinatorial_enumeration(self):
        """returns a list of new components that are copies of
        existing components, but with a different species placed inside. This different
        species represents different combinatorial states of the polymer.
        for example:
            construct: <A,B,C>
        two new species are possible: <[A:RNAP],B,C>; <A,[B:RNAP],C>
        combinatorial species is also possible (since A and B are assumed to act independantly)
                                        <[A:RNAP],[B:RNAP],C>
        Thus, this function returns A which binds to <A,B,C> (creating <[A:RNAP],B,C>) AND also
        A which binds to <A,[B:RNAP],C> (creating <[A:RNAP],[B:RNAP],C>). Likewise for B
        In total two A components are returned, and two B components are returned.
        """
        #Looks at combinatorial states of constructs to generate DNA_parts
        #my_polymer = self.get_species()
        self.update_parameters()

        #Go through parts
        active_components = []
        for part in self.parts_list:
            if(hasattr(part,"update_component")):
                dummy_species = self.get_species()[part.position]
                updated_components = part.update_component(dummy_species,practice_run=True)
                if(updated_components is not None):
                    active_components += [updated_components]
        #this next part creates "combinatorial" bound complexes, given singly bound complexes generated above.
        #for example, let's say we got <[A:RNAP],B,C> and <A,[B:RNAP],C> from A and B individually, then a combinatorial
        #construct would be <[A:RNAP],[B:RNAP],C>
        combinatorial_complexes = self.update_combinatorial_complexes(active_components)

        combinatorial_components = []
        for comb_specie in combinatorial_complexes:
            #after making all combinatorially bound parts, we must seed components with the right species
            #so that the right reactions are made. For example, A cannot react with <[A:RNAP],[B:RNAP],C> because
            #that's the species you get after A has already reacted. So here we are finding species with the
            #proper positions unbound, so they can be fed to the proper components
            if(isinstance(comb_specie,OrderedPolymerSpecies)):
                for part in active_components:
                    part_pos = part.position
                    if(not isinstance(comb_specie[part_pos],ComplexSpecies)):
                        #in this case the position of interest is not complexed. Thus, we need to update!
                        combinatorial_components += [part.update_component(comb_specie[part_pos])]
        return combinatorial_components

    def enumerate_components(self,previously_enumerated=None):
        """returns a list of new components and constructs.
        New components are generated if:
        - a component creates a species which results in binding to part of the construct
            Example: <A,B,C> -> <[A:RNAP],B,C>
            Then, A would be returned since a new species is created
        - more than one such component exist in the same construct, for example:
            construct: <A,B,C>
            two new species are possible: <[A:RNAP],B,C>; <A,[B:RNAP],C>
            combinatorial species is also possible (since A and B are assumed to act independantly)
                                         <[A:RNAP],[B:RNAP],C>
            Thus, this function returns A which binds to <A,B,C> (creating <[A:RNAP],B,C>) AND also
            A which binds to <A,[B:RNAP],C> (creating <[A:RNAP],[B:RNAP],C>). Likewise for B
            In total two A components are returned, and two B components are returned.
        New constructs are generated if:
            self.enumerate_construcs() says so.
            For example, in <A,B,C>, A is a promoter and makes an RNA_construct containing <B,C>
        """
        #Runs component enumerator to generate new constructs
        new_constructs = self.enumerate_constructs(previously_enumerated=previously_enumerated)

        #Looks at combinatorial states of constructs to generate DNA_parts
        combinatorial_components = self.combinatorial_enumeration()

        return combinatorial_components+new_constructs


    def __hash__(self):
        return OrderedPolymer.__hash__(self)
    def __eq__(self,construct2):
        """equality means comparing the parts list in a way that is not too deep"""
        if(self.__repr__()==construct2.__repr__()):
            return True
        else:
            return False
    
    def update_species(self):
        species = [self.get_species()]
        return species
    def reset_stored_data(self):
        self.out_components = None
        self.predicted_rnas = None
        self.predicted_proteins = None
    def changed(self):
        self.reset_stored_data()
        self.name = self.make_name()
    def update_reactions(self,norna=False):
        return []
        


class DNA_construct(Construct,DNA):
    def __init__(self,
                parts_list,
                name=None,
                circular=False,
                mechanisms=None,  # custom mechanisms
                parameters=None,  # customized parameters
                attributes=None,
                initial_concentration=None,
                copy_parts=True,
                component_enumerators = None,
                **keywords):

        self.material_type = "dna"
        if component_enumerators is None:
            from .construct_explorer import TxExplorer
            component_enumerators = [TxExplorer()]
        
        Construct.__init__(self=self, parts_list =parts_list, name = name, \
                            circular=circular, mechanisms=mechanisms, \
                            parameters=parameters, attributes=attributes, \
                            initial_concentration=initial_concentration, \
                            component_enumerators = component_enumerators, **keywords)
        DNA.__init__(self=self,name=self.name)
    def __repr__(self):
        return "DNA_construct = "+ self.make_name()
    
        
class RNA_construct(Construct,RNA):
    def __init__(self,parts_list,name=None,promoter=None,\
                component_enumerators = None,\
                **keywords):
        """an RNA_construct is a lot like a DNA_construct except it can only translate, and
        can only be linear"""
        
        self.material_type = "rna"
        self.promoter = promoter
        if component_enumerators is None:
            from .construct_explorer import TlExplorer
            component_enumerators = [TlExplorer()]

        Construct.__init__(self=self,parts_list=parts_list,circular=False,name=name,\
                                component_enumerators = component_enumerators,**keywords)
        RNA.__init__(self=self,name=self.name)
    def __repr__(self):
        """the name of an RNA should be different from DNA, right?"""
        return "RNA_construct = "+self.name
        