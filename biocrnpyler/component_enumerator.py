
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.
from typing import List, Union
from warnings import warn
from .component import Component
from .dna_part_cds import CDS
from .dna_part_misc import AttachmentSite
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .dna_part_terminator import Terminator
from dna_construct import RNA_construct, DNA_construct
import copy
from .utils import rev_dir

class ComponentEnumerator:
    def __init__(self, name:str,enumerator_type=""):
        """a component enumerator's job is to create new components in a process similar to mechanisms."""
        self.name = name
        self.enumerator_type = enumerator_type
    def update_components(self,component=None) -> List:
        """this will create new components based on the input component somehow
        The child class should implement this
        :return: empty list"""
        warn(f"Default update_components called for ComponentEnumerator = {self.name}.")
        return []
    def __repr__(self):
        return self.name

class TxTlExplorer(ComponentEnumerator):
    def __init__(self,name="TxTlExplorer",enumerator_type="local",possible_rxns=("transcription","translation"),direction="forward"):
        """this class goes through a parts_list of a DNA_construct and decides what RNAs are made
        and what proteins are made based on orientation and location of parts"""
        self.current_rnas = {}
        self.current_proteins = {}
        self.made_rnas = {}
        self.all_rnas = {}
        self.made_proteins = {}
        self.current_rbs = None
        self.possible_rxns = possible_rxns
        self.direction=direction
        self.second_looping = False
    def update_components(self,component=None,possible_rxns=("transcription","translation")):
        #TODO the component we are given better be a dna_construct or something that can be txtled
        proteins = {}
        rnas = {}
        for direction in ["forward","reverse"]:
            self.__init__(possible_rxns=possible_rxns,direction=direction)
            startnum = 0
            newlist = component.parts_list
            if(direction == "reverse"):
                startnum = len(component.parts_list)-1
                #if we go backwards then also list the parts backwards
                #deepcopy so we don't mangle the list
                #newlist = copy.deepcopy(component.parts_list)[::-1]
            else:
                pass
                #newlist = copy.deepcopy(component.parts_list)
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
                        #print('restart from beginning')
                        second_looping = 1 #during the second loop, we don't care about promoters
                        self.second_loop() #tell the explorer that we don't care about promoters
                    else:
                        self.end() #this completes all "in progress" RNAs and proteins
                        break
            #proteins.update(explorer.get_proteins())
            rnas.update(self.get_rnas())
        proteins = {}
        for promoter in rnas:
            _,prots = rnas[promoter].explore_txtl()
            proteins.update(prots)
        #the output of this function should be a list of new components.
        #but txtl explorer knows more than that! it knows who made what also
        #perhaps these reactions should be *inside* the components which we spit out.
        #
        #we could update the appropriate RBS or promoter parts with the species they are
        #making, then output those components as new components.
        return rnas,proteins
    def see(self,part):
        """the explorer sees a part. This does different stuff depending on
        1) if there is a current rna or current rbs
        2) what the part is
        it returns whether:
        1: not done, keep going
        0: done, stop
        """
        effective_direction = part.direction
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
                (isinstance(part,Promoter) and "transcription" in self.possible_rxns)):
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
        if(self.second_looping==False):
            new_rxns = []
            for rxn in self.possible_rxns:
                if(rxn!="transcription"):
                    #remove transcription from the set of possible reactions!
                    new_rxns += [rxn]
            self.possible_rxns = new_rxns
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
        if(promoter==None and len(self.current_rnas)>0):
            #if you don't specify then it terminates the first one
            promoter = list(self.current_rnas.keys())[0]
        elif(len(self.current_rnas)==0):
            #in this case you terminated transcription when nothing was transcribing
            return
        rna_partslist = self.current_rnas[promoter][0]
        #print(self.current_rnas)
        #print(self.current_proteins)
        #TODO copy parts here and not in the constructor
        rna_construct = RNA_construct(copy.deepcopy(rna_partslist),made_by = promoter)
        rna_construct.set_mixture(promoter.mixture)
        #current_rna_name = str(promoter)+"-"+str(rna_construct)
        self.made_proteins[rna_construct]={}
        #compile the name, keeping track of which promoter made the rna and all the parts on it
        proteins_per_promoter = []
        for rbs in self.current_rnas[promoter][1]:
            #ending the RNA also ends all the proteins being generated here.
            proteins_per_rbs = []
            # TODO this will run multiple times for each RNA. make it so it runs once!
            for protein_part in self.current_proteins[rbs]:
                if(protein_part[1] == "forward"):
                    #it is essential to be forwards
                    if(isinstance(protein_part[0],CDS) and protein_part[0].protein != None):
                        #this happens if we do in fact make a protein
                        proteins_per_rbs+=[protein_part[0]]
            #print(rna_partslist)
            #print("currently we are translating from "+str(rbs)+ " which is located at " + str(rbs.pos))
            #TODO correct_rbs is possibly not needed
            #print(rbs)
            #print(rna_partslist)
            correct_rbs = rna_construct.parts_list[rna_partslist.index([rbs,"forward"])]
            proteins_per_promoter+=proteins_per_rbs
            self.made_proteins[rna_construct].update({correct_rbs:proteins_per_rbs})
            correct_rbs.protein = proteins_per_rbs
            # so this has to be done at the end of the RNA only because only then do we know
            # exactly what that full RNA is going to be.
            # we want to clear the current proteins, but there could be multiple RNAs that run over the same proteins,
            # so wait until we are done tallying all the RNAs before removing everything from current_proteins
        del self.current_rnas[promoter] # this removes the current RNA from the list, because it's
                                        # being terminated!!
        promoter.transcript = rna_construct.get_species() #set the promoter's transcript to be correct
        promoter.protein = proteins_per_promoter
        self.all_rnas.update({promoter:rna_construct})
        return
    def end(self):
        """we've reached the end of the dna! End everything!"""
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