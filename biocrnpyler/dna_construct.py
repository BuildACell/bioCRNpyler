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
from .dna_part_misc import AttachmentSite
import biocrnpyler.component_enumerator as ce

from .component_enumerator import ComponentEnumerator,LocalComponentEnumerator, GlobalComponentEnumerator
from .species import (ComplexSpecies, OrderedMonomer, OrderedPolymer,
                      OrderedPolymerSpecies)
from .utils import all_comb, remove_bindloc, rev_dir
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
                copy_parts=True,
                component_enumerators = None,
                **keywords):
        """this represents a bunch of parts in a row.
        A parts list has [[part,direction],[part,direction],...]
        Each part must be an OrderedMonomer"""
        if(component_enumerators is None):
            component_enumerators = []
        self.component_enumerators = component_enumerators
        '''
        myparts = []
        
        if(copy_parts):
            parts_list  = copy.deepcopy(parts_list)
            
        for part in parts_list:
            newpart = []
            if(type(part)==list or type(part)==tuple):
                newpart = [part[0].remove(),part[1]]
                #if(copy_parts):
                #    newpart = [part[0].remove(),part[1]]
                #else:
                #    newpart = [part[0].remove(),part[1]]
            elif(isinstance(part,OrderedMonomer)):
                npartd = None
                npart = None
                if(part.direction is not None):
                    #remember the direction that is in the part
                    npartd = copy.deepcopy(part.direction)
                else:
                    #if no direction is specified, forward is default
                    npartd = "forward"
                if(copy_parts):
                    #we have to copy the part
                    npart = copy.deepcopy(part)
                else:
                    npart = part
                #forget everything about being part of a polymer
                npart.remove()
                #for the purposes of OrderedPolymer, you still need [part,direction]
                newpart = [npart,npartd]
            myparts += [newpart]
        
        OrderedPolymer.__init__(self,myparts)
        #'''
        OrderedPolymer.__init__(self,parts_list,default_direction="forward")
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
        self.update_base_species(self.name)
        self.out_components = None
        self.predicted_rnas = None
        self.predicted_proteins = None

    @property
    def parts_list(self):
        return self._polymer

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
        self.update_base_species()
        return self

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
            partspec.material_type = self.material_type
            ocomplx += [partspec.set_dir(part.direction)]
        out_species = OrderedPolymerSpecies(ocomplx,base_species = self.base_species,circular = self.circular,\
                                                    name = self.name,material_type=self.material_type)
        
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
        spec_list = [a[0] for a in spec_list]
        spec_indexes = [a.position for a in spec_list] #extract all indexes
        #print(spec_indexes)
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
        return outlist

    def make_polymers(self,species_lists,backbone):
        """makes polymers from lists of species
        inputs:
        species_lists: list of species which are to be assembled into a polymer
        backbone: the base_species which all these polymers should have"""
        polymers = []
        for combo in species_lists:
            #members of allcomb are now OrderedMonomers, which contain direction and position
            #there could be multiple OrderedPolymerSpecies we are making combinatorial.
            #for example, RNAs
            new_backbone = copy.deepcopy(backbone)
            new_material = backbone.material_type
            for spec in combo:
                new_material = OrderedPolymerSpecies.default_material
                new_backbone.replace(spec.position,spec)
            newname = None
            self_species = self.get_species()
            if(new_backbone==self_species):
                #if we have just re-created ourselves, then make sure to call it that
                newname = self_species.name
                polymers += [copy.deepcopy(self_species)]
            else:
                new_backbone.material_type = new_material
                polymers += [new_backbone] #we make a new OrderedComplexSpecies
        return polymers

    def update_combinatorial_complexes(self,active_components):
        """given an input list of components, we produce all complexes
        yielded by those components, mixed and matched to make all possible combinatorial
        complexes, where each component is assumed to only care about binding to one spot"""
        species = [self.get_species()]
        for part in active_components:
            #first we make binary complexes
            sp_list =  part.update_species()
            species+=remove_bindloc(sp_list)
        #print("initial species")
        #print(species)
        unique_complexes = {}
        possible_backbones = {self.base_species:self.get_species()}
        #possible_backbones = {a.name:a.get_species() for a in [self]+[a for a in proteins]}
        #species need to be uniqueified
        #print(species)
        unique_species = list(set(species)) 
        for specie in unique_species:
            #in this list we extract all the variants of the complexes from possible_backbones
            #that exist in our species list.
            if(isinstance(specie,OrderedPolymerSpecies) and specie.base_species in possible_backbones):
                #we only care about OrderedPolymerSpecies made from this construct
                if(specie.base_species in unique_complexes):
                    unique_complexes[specie.base_species] += [specie]
                else:
                    unique_complexes[specie.base_species] = [specie]
        #unique_complexes now has a list of all the non-combinatorial complexes we can make
        combinatorial_complexes = []
        for bb_name in unique_complexes:
            #for each backbone, make all combinatorial combinations.
            comp_binders = []
            for unit in unique_complexes[bb_name]:
                #for each complex for this backbone, find out what is bound at which location
                comp_bound = []
                pos_i = 0
                for pos in unit:
                    #find the position that has a ComplexSpecies in it
                    if(isinstance(pos,ComplexSpecies)):
                        comp_bound += [copy.deepcopy(pos)] 
                    pos_i+=1
                if(len(comp_bound) >0):
                    comp_binders += [comp_bound] #record what is bound, and at what position
            #comp_binders is a list of lists because multiple things can be bound at the same time
            allcomb = self.located_allcomb(comp_binders) #all possible combinations of binders are made here
            allcomb += [[]] #unbound dna should also be used
            #now, all possibilities have been enumerated.
            #we construct the OrderedPolymerSpecies
            combinatorial_complexes += self.make_polymers(allcomb,possible_backbones[bb_name])
            
        return combinatorial_complexes

    #Overwrite Component.enumerate_components 
    def enumerate_constructs(self):
        #Runs component enumerator to generate new constructs
        new_constructs = []
        for enumerator in self.component_enumerators:
            new_comp = enumerator.enumerate_components(component=self)
            new_constructs += new_comp
        return new_constructs

    def combinatorial_enumeration(self):
        #Looks at combinatorial states of constructs to generate DNA_parts
        multivalent_self = self.get_species()
        self.update_parameters()


        #Go through parts
        active_components = []
        for part in self.parts_list:
            if(hasattr(part,"update_component")):
                updated_components = part.update_component(multivalent_self[part.position])
                if(updated_components is not None):
                    active_components += [updated_components]
        combinatorial_complexes = self.update_combinatorial_complexes(active_components)
        combinatorial_components = []
        for comb_specie in combinatorial_complexes:
            if(isinstance(comb_specie,OrderedPolymerSpecies) and comb_specie.base_species == self.base_species):
                for part in active_components:
                    part_pos = part.position
                    if(isinstance(comb_specie[part_pos],ComplexSpecies)):
                        #in this case the position of interest is already complexed. Skip!
                        pass
                    else:
                        combinatorial_components += [part.update_component(comb_specie[part_pos])]
        return combinatorial_components

    def enumerate_components(self):
        #Runs component enumerator to generate new constructs
        new_constructs = self.enumerate_constructs()

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
                            copy_parts=copy_parts, \
                            component_enumerators = component_enumerators, **keywords)

        pind = 0
        for part in self.parts_list:
            if(type(part.color)==type(None)):
                #if the color isn't set, let's set it now!
                part.color = pind
            if(type(part.color2)==type(None)):
                if(isinstance(part,AttachmentSite) and part.site_type in ["attL","attR"]):
                    #this is the only scenario in which we need a color2
                    #TODO make this more general
                    pind+=1
                    part.color2 = pind
            pind+=1
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
    #def get_species(self):
    #    outspec = Construct.get_species(self)
        #for spec in outspec:
        #    spec.material_type = "rna"
        #outspec.material_type="rna"
    #    return outspec
    def __repr__(self):
        """the name of an RNA should be different from DNA, right?"""
        return "RNA_construct = "+self.name