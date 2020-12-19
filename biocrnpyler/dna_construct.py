################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 7/31/2020
#
################################################################


# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from .component import Component
from .components_basic import DNA, RNA
from .dna_part import DNA_part
from .dna_part_cds import CDS
from .dna_part_misc import AttachmentSite
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .dna_part_terminator import Terminator
import biocrnpyler.component_enumerator as ce
from .component_enumerator import ComponentEnumerator,LocalComponentEnumerator, GlobalComponentEnumerator
from .species import (ComplexSpecies, OrderedMonomer, OrderedPolymer,
                      OrderedPolymerSpecies)
from .utils import all_comb, remove_bindloc, rev_dir
import logging


class ConstructExplorer(ComponentEnumerator):
    def __init__(self, name, direction="forward",possible_directions=("forward","reverse"), max_loop_count = 2):
        ComponentEnumerator.__init__(self=self,name=name)
        self.direction=direction
        self.max_loop_count = max_loop_count
        self.possible_directions = possible_directions

    def enumerate_components(self, component):
        #Only works on Constructs!
        if isinstance(component, Construct):
            #iterate through the possible directions
            for direction in self.possible_directions:
                #Set the loop count
                self.current_loop_count = 0 #due to while loop, current_loop_count is instantly increment to 1
                self.initialize_loop(direction = direction)

                #Deep Copy Component
                parts_list = list(copy.deepcopy(component.parts_list))
                #Set the direction
                if direction == "reverse":
                    parts_list.reverse()

                #Main iteration loop
                while self.check_loop():
                    #Check each part
                    for part_ind, part in enumerate(parts_list):
                        #Do something to the part
                        self.iterate_part(part, direction)

            return self.return_components()
        else:
            return []

    def check_loop(self):
        """if we already went around the plasmid, then what we're checking for is continuing
        transcripts or proteins. We don't want to start making new transcripts
        because we already checked this area for promoters"""
        logging.debug(f"loop count = {self.current_loop_count}")
        if self.current_loop_count > self.max_loop_count:
            self.terminate_loop()
            return False
        else:
            self.current_loop_count += 1
            return True

    def initialize_loop(self, direction):
        #Starts the loop
        #MUST SUBCLASS
        pass

    def iterate_part(self, part, direction):
        #Part is a DNA_part
        #direction is the direction we are looking in along the DNA_construct
        #runs on every part and does something!
        #MUST SUBCLASS
        pass

    def terminate_loop(self):
        #called when we are done enumerating
        pass
        #MUST SUBCLASS

    def return_components(self):
        #returns components at the end of the loop
        #MUST SUBCLASS
        return []

class TxExplorer(ConstructExplorer):
    def __init__(self, name = "TxExplorer", direction="forward",possible_directions=("forward","reverse"), max_loop_count = 2):
        ConstructExplorer.__init__(self, name, direction="forward",possible_directions=("forward","reverse"), max_loop_count = 2)

        self.all_rnas = {} #Stores all transcripts from the DNA_Construct


    def initialize_loop(self, direction):
        self.current_rnas = {} #Stores RNA's being examined in loop promoter --> [(DNA_part, direction) list]


    def iterate_part(self, part, reading_direction):
        #Part is a DNA_part
        #direction is the direction we are looking in along the DNA_construct
        
        """the explorer sees a dna_part.
        1. it calculates effective direction of part relative to transcription direction
        2. it adds parts to growing transcripts
        3. it makes new transcripts and terminates growing transcripts.
        """
        logging.debug("looking at "+str(part))
        logging.debug("current_rnas are "+str(self.current_rnas))

        #The orientation of the part absolutely to the DNA_construct
        part_direction = part.direction

        #The relative orientation of the part to the reading direction
        if part_direction == reading_direction:
            effective_direction = "forward"
        else:
            effective_direction = "reverse"

        #Grow Existing Transcripts
        #For each promoter being transcribed, we will add the part and its direction relative to that promoter
        for promoter in self.current_rnas:
            self.current_rnas[promoter]+= [(part, effective_direction)]

        #Case for Different Parts doing different things

        #Case 1: Promoter in the correct direction (and first loop around a circular plasmid)
        if effective_direction=="forward" and isinstance(part,Promoter) and self.current_loop_count == 1:
            #this if statement makes sure that the current part wants to transcribe, and
            #that is something that we are allowed to do
            self.current_rnas[part] = []

        #Case 2: Terminator
        elif effective_direction == "forward" and isinstance(part, Terminator): #Kill Sarah Day O'Connor
            self.terminate_loop()
    

    def terminate_loop(self):
        #Transfers current rnas into all_rnas

        #Create RNA_Constructs for each promoter

        for promoter in self.current_rnas:
            rna_parts_list = self.current_rnas[promoter]
            rna_construct = RNA_construct(rna_parts_list, promoter = promoter)
            self.all_rnas[promoter] = rna_construct
            promoter.transcript = rna_construct.get_species()

        self.current_rnas = {}

    #Returns a list of RNAconstructs
    def return_components(self):
        print("TxEnumerator returning,", [self.all_rnas[k] for k in self.all_rnas])
        return [self.all_rnas[k] for k in self.all_rnas]



class TxTlExplorer_CE(LocalComponentEnumerator):
    def __init__(self,name="TxTlExplorer",possible_rxns=("transcription","translation"),\
                    direction="forward",possible_directions=("forward","reverse"),\
                        return_rnas = True, return_proteins = True):
        """this class goes through a parts_list of a DNA_construct and decides what RNAs are made
        and what proteins are made based on orientation and location of parts"""
        LocalComponentEnumerator.__init__(self=self,name=name)
        self.return_rnas = return_rnas
        self.return_proteins = return_proteins
        self.current_rnas = {}
        self.current_proteins = {}
        self.made_rnas = {}
        self.all_rnas = {}
        self.made_proteins = {}
        self.current_rbs = None
        self.possible_rxns = possible_rxns
        self.direction=direction
        self.second_looping = False
        self.possible_directions = possible_directions

    def enumerate_components(self,component):
        assert(isinstance(component,Construct))
        logging.debug("enumerating "+str(component))
        proteins = {}
        rnas = {}
        #copy the components so 
        newlist = copy.deepcopy(component.parts_list)
        for direction in self.possible_directions:
            self.__init__(direction=direction,possible_rxns=self.possible_rxns,possible_directions=self.possible_directions)
            logging.debug("direction is "+str(direction))
            if("transcription" not in self.possible_rxns and component.material_type == "rna"):
                #if we are looking at an RNA, then track the fact that it is an RNA (and properly evaluate RNA_specific parts)
                    self.make_rna(component.my_promoter)
            startnum = 0
            if(direction == "reverse"):
                startnum = len(component.parts_list)-1
                #if we go backwards then also list the parts backwards
            elif(direction not in ["forward","reverse"]):
                raise ValueError("Got incoherent direction "+str(direction)+" when it should be either forward or reverse")
            part_index = startnum
            keep_going = 1
            second_looping = 0
            while keep_going: #we may have to loop through the parts twice because the construct is circular.
                                #this does not account for infinite loops, RCA, etc.
                part = newlist[part_index] #pull out the part
                keep_going = self.see(part) #the explorer keeps track of everything
                if(direction=="forward"):
                    part_index+=1 #keep going+
                elif(direction=="reverse"):
                    part_index-=1 #go backwards

                if(part_index==len(component.parts_list) or part_index<0):
                    part_index = startnum #this takes care of looping
                    if(component.circular and second_looping == 0):
                        second_looping = 1 #during the second loop, we don't care about promoters
                        self.second_loop() #tell the explorer that we don't care about promoters
                    else:
                        self.end() #this completes all "in progress" RNAs and proteins
                        break
            proteins.update(self.get_proteins())
            rnas.update(self.get_rnas())
        protein_components = []
        #the output of this function should be a list of new components.
        #but txtl explorer knows more than that! it knows who made what also
        #perhaps these reactions should be *inside* the components which we spit out.
        #
        #we will update the appropriate RBS or promoter parts with the species they are
        #making, then output those components as new components.
        rna_constructs = [rnas[a] for a in rnas]
        #this next part is configuring the output list with the qualities of the previous components
        for partind in range(len(newlist)):
            if(hasattr(newlist[partind],"transcript")):
                component[partind].transcript = newlist[partind].transcript
            if(hasattr(newlist[partind],"protein")):
                component[partind].protein = newlist[partind].protein
        
        for rna in proteins:
            for rbs in proteins[rna]:
                protein_components += proteins[rna][rbs]
        outlist = []
        ######################################
        if(self.return_rnas):
            outlist += rna_constructs
        if(self.return_proteins):
            outlist += protein_components
        ######################################
        return outlist
    def see(self,part):
        """the explorer sees a part. This does different stuff depending on
        1) if there is a current rna or current rbs
        2) what the part is
        it returns whether:
        1: not done, keep going
        0: done, stop
        """
        effective_direction = part.direction
        logging.debug("looking at "+str(part))
        logging.debug("possible reactions are "+str(self.possible_rxns))
        logging.debug("current_rnas are "+str(self.current_rnas))
        if(self.direction=="reverse"):
            effective_direction = rev_dir(effective_direction)
        if(self.current_rnas!={}):
            #this means we are making an RNA
            terminated = 0
            #TODO how do we know that promoters transcribe and RBSes translate?
            #maybe you can transcribe or translate in both directions?
            if(isinstance(part,RBS) and "translation" in self.possible_rxns and \
                                                    effective_direction == "forward"):
                
                self.current_proteins[part]= []
                for promoter in self.current_rnas:
                    #this tells us that whatever RNAs are currently being made contain this RBS
                    if(part not in self.current_rnas[promoter][1]):
                        self.current_rnas[promoter][1] += [part]
                if(self.current_rbs == None):
                    #this means nobody is translating yet
                    self.current_rbs = part
                else:
                    #another RBS is translating. this could mean a few things:
                    #1: translational coupling. the previous translation ends and turns
                    #   into this one, with an added coupling term when that gets implemented
                    #2: it terminates the previous translation. Chances are it would, right?
                    self.current_rbs = part
                    # TODO add coupling
            elif(effective_direction in part.no_stop_codons and self.current_rbs != None):
                #this means we can translate through the current part!
                self.current_proteins[self.current_rbs] += [[part,effective_direction]]
            elif(self.current_rbs != None):
                #this means we CAN'T translate through the current part, but we are translating
                self.current_proteins[self.current_rbs] += [[part,effective_direction]]
                self.current_rbs = None #current RBS has done its work
            for promoter in list(self.current_rnas.keys()):
                self.current_rnas[promoter][0]+=[[part,effective_direction]]
                #we add the current part
                if(isinstance(part,Terminator) and effective_direction == "forward"):
                    terminated = 1 #a terminator ends everything. This does not account for leakiness
                    self.terminate_transcription(promoter)
            if(terminated):
                assert(self.current_rnas=={}) #this should be true if terminate_transcription() worked
                # this means that all RNAs have been stopped. But, since the proteins didn't know
                # which RNAs they belonged to, they are still hanging around so we have to
                # clean them up
                self.current_proteins = {}
                self.current_rbs = None
                terminated = 0
        if((effective_direction=="forward") and \
                isinstance(part,Promoter) and "transcription" in self.possible_rxns and not self.second_looping):
            #this if statement makes sure that the current part wants to transcribe, and
            #that is something that we are allowed to do
            self.current_rnas.update({part:[[],[]]})
            return 1
        elif((("transcription" not in self.possible_rxns) or self.second_looping) and \
                                         self.current_rnas=={}):
            return 0 #this means we are pretty sure there is nothing else left to do.
        else:
            return 1 #this means keep going! As far as we know there could be more to see
    def second_loop(self):
        """if we already went around the plasmid, then what we're checking for is continuing
        transcripts or proteins. We don't want to start making new transcripts
        because we already checked this area for promoters"""
        logging.debug("second loop!")
        if(self.second_looping==False):
            #new_rxns = []
            #for rxn in self.possible_rxns:
            #    if(rxn!="transcription"):
            #        #remove transcription from the set of possible reactions!
            #        new_rxns += [rxn]
            #self.possible_rxns = new_rxns
            self.second_looping=True
        else:
            self.end()
    def make_rna(self,promname="external"):
        """we are making an rna! imposed by external force"""
        if(promname in self.current_rnas):
            warn("RNA already started; it looks like this: "+str(self.current_rnas[promname]))
        else:
            self.current_rnas.update({promname:[[],[]]})
    def terminate_transcription(self,promoter=None):
        """this terminates the transcription of a specific RNA, started by the promoter
        given in the argument. the RNA is deleted from the list of currently transcribing RNAs,
        and any proteins that it was making are catalogued"""
        logging.debug("transcription terminated! Promoter was "+str(promoter))
        if(promoter==None and len(self.current_rnas)>0):
            #if you don't specify then it terminates the first one
            promoter = list(self.current_rnas.keys())[0]
        elif(len(self.current_rnas)==0):
            #in this case you terminated transcription when nothing was transcribing
            return
        rna_partslist = self.current_rnas[promoter][0]
        copied_partslist = copy.deepcopy(rna_partslist)
        rna_construct = RNA_construct(copied_partslist,made_by = promoter)
        rna_construct.set_mixture(promoter.mixture)
        self.made_proteins[rna_construct]={}
        #compile the name, keeping track of which promoter made the rna and all the parts on it
        proteins_per_promoter = []
        for rbs in self.current_rnas[promoter][1]:
            #ending the RNA also ends all the proteins being generated here.
            proteins_per_rbs = []
            for protein_part in self.current_proteins[rbs]:
                if(protein_part[1] == "forward"):
                    #it is essential to be forwards
                    if(isinstance(protein_part[0],CDS) and protein_part[0].protein != None):
                        #this happens if we do in fact make a protein
                        proteins_per_rbs+=[protein_part[0]]
            if(self.return_rnas):
                #this means we are creating new RNAs, and so the correct RBS is the one in the new RNA.
                correct_rbs = rna_construct.parts_list[rna_partslist.index([rbs,"forward"])]
            else:
                #this means we are reading an existing RNA, so then the correct rbs is the one
                #in the original "DNA", which is actually an RNA.
                correct_rbs = rbs
            proteins_per_promoter+=proteins_per_rbs
            self.made_proteins[rna_construct].update({correct_rbs:proteins_per_rbs})
            correct_rbs.protein = [a.get_species() for a in proteins_per_rbs]
            # so this has to be done at the end of the RNA only because only then do we know
            # exactly what that full RNA is going to be.
            # we want to clear the current proteins, but there could be multiple RNAs that run over the same proteins,
            # so wait until we are done tallying all the RNAs before removing everything from current_proteins
        del self.current_rnas[promoter] # this removes the current RNA from the list, because it's
                                        # being terminated!!
        if(self.return_rnas):
            promoter.transcript = rna_construct.get_species() #set the promoter's transcript to be correct
        promoter_proteins = [a.get_species() for a in proteins_per_promoter]
        if(len(promoter_proteins)>0):
            promoter.protein = promoter_proteins
        self.all_rnas.update({promoter:rna_construct})
        return
    def end(self):
        """we've reached the end of the dna! End everything!"""
        logging.debug("end of the DNA!")
        for promoter in list(self.current_rnas.keys()):
            self.terminate_transcription(promoter)
        self.current_rnas = {}
        self.current_proteins = {}
        self.current_rbs = None
    def get_proteins(self):
        """return all the proteins made during the latest exploration"""
        return self.made_proteins
    def get_rnas(self):
        """return all RNAs made during the latest exploration"""
        return self.all_rnas

class TxExplorerOLD(TxTlExplorer_CE):
    def __init__(self,name="TxExplorer",possible_rxns=("transcription",),\
                        direction="forward",possible_directions=("forward","reverse"),debug=False):
        TxTlExplorer_CE.__init__(self,name=name,possible_rxns=possible_rxns,\
                        direction=direction,possible_directions=possible_directions,return_proteins=False)

class TlExplorer(TxTlExplorer_CE):
    def __init__(self,name="TlExplorer",possible_rxns=("translation",),\
                        direction="forward",possible_directions=("forward",),debug=False):
        TxTlExplorer_CE.__init__(self,name=name,possible_rxns=possible_rxns,\
                        direction=direction,possible_directions=possible_directions,return_rnas = False)

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
        myparts = []
        for part in parts_list:
            newpart = []
            if(type(part)==list or type(part)==tuple):
                if(copy_parts):
                    newpart = [copy.deepcopy(part[0]).remove(),part[1]]
                else:
                    newpart = [part[0].remove(),part[1]]
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

    def show(self):
        txt = self.name
        for part in self.parts_list:
            if(isinstance(part,Promoter)):
                if(part.transcript is not None):
                    txt += "\n\t"+repr(part)
                    txt += "\n\t" + repr(part.transcript)
                if(part.protein is not None):
                    for prot in part.protein:
                        txt += "\n\tprotein = "+repr(prot)
        return txt

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
                new_obj2 = copy.deepcopy(obj2).unclone()
                uncloned_list = [a.unclone() for a in copy.deepcopy(self.parts_list)]
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
        return hash(self.__repr__())
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
                component_enumerators = (TxExplorer(),),
                **keywords):

        self.material_type = "dna"
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
                component_enumerators = (TlExplorer(),),\
                **keywords):
        """an RNA_construct is a lot like a DNA_construct except it can only translate, and
        can only be linear"""
        self.material_type = "rna"
        self.my_promoter = promoter
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