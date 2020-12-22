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
from .dna_part_cds import CDS
from .dna_part_misc import AttachmentSite
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .dna_part_terminator import Terminator
from .dna_construct import Construct, RNA_construct
import biocrnpyler.component_enumerator as ce
from .component_enumerator import ComponentEnumerator,LocalComponentEnumerator, GlobalComponentEnumerator
from .species import (ComplexSpecies, OrderedMonomer, OrderedPolymer,
                      OrderedPolymerSpecies)
from .utils import all_comb, remove_bindloc, rev_dir
import logging


class ConstructExplorer(LocalComponentEnumerator):
    def __init__(self, name, direction="forward",possible_directions=("forward","reverse")):
        ComponentEnumerator.__init__(self=self,name=name)
        self.direction=direction
        self.possible_directions = possible_directions
        self.reset_enumerator()

    def reset_enumerator(self):
        """clear accumulator variables from the numerator. MUST BE SUBCLASSED"""
        pass

    def enumerate_components(self, component):
        #Only works on Constructs!
        self.reset_enumerator()
        if isinstance(component, Construct):

            #if we are circular, then we can go around twice
            #otherwise, once is the correct amount
            if(component.circular):
                self.max_loop_count = 2
            else:
                self.max_loop_count = 1
            #iterate through the possible directions
            for direction in self.possible_directions:
                #Set the loop count
                self.current_loop_count = 0 #due to while loop, current_loop_count is instantly increment to 1
                self.initialize_loop(direction = direction)

                #Deep Copy Component
                parts_list = list(component.parts_list)
                #Set the direction
                if direction == "reverse":
                    parts_list.reverse()
                #Main iteration loop
                while self.check_loop():
                    #Check each part
                    for part_ind, part in enumerate(parts_list):
                        #Do something to the part
                        self.iterate_part(part, direction)
                
                #must stop before switching directions
                self.terminate_loop()

            return self.return_components(component)
        else:
            return []

    def check_loop(self):
        """if we already went around the plasmid, then what we're checking for is continuing
        transcripts or proteins. We don't want to start making new transcripts
        because we already checked this area for promoters"""
        logging.debug(f"loop count = {self.current_loop_count}")
        if self.current_loop_count >= self.max_loop_count:
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

    def return_components(self,component):
        #returns components at the end of the loop
        #MUST SUBCLASS
        return []

class TxExplorer(ConstructExplorer):
    def __init__(self, name = "TxExplorer", direction="forward",possible_directions=("forward","reverse")):
        ConstructExplorer.__init__(self, name, direction=direction,possible_directions=possible_directions)


    def reset_enumerator(self):
        self.all_rnas = {} #must clear this before filling it up again

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
            if(len(rna_parts_list) > 0):
                rna_construct = RNA_construct(rna_parts_list, promoter = promoter)
                self.all_rnas[promoter] = rna_construct
                promoter.transcript = rna_construct.get_species()
            else:
                warn(f"{promoter} makes empty transcript! Protein binding to this promoter will not be included in the CRN")
        self.current_rnas = {}

    #Returns a list of RNAconstructs
    def return_components(self,component):
        rna_construct_list = [self.all_rnas[k] for k in self.all_rnas]
        return rna_construct_list

class TlExplorer(ConstructExplorer):
    def __init__(self, name = "TlExplorer", direction="forward",possible_directions=("forward",)):
        ConstructExplorer.__init__(self, name, direction=direction,possible_directions=possible_directions)
    def reset_enumerator(self):
        self.all_proteins = {} #Stores all proteins from the RNA_construct
    def initialize_loop(self, direction):
        self.current_proteins = {} #Stores proteins being examined in loop rbs --> [(DNA_part, direction) list]
    def iterate_part(self, part, reading_direction):
        #Part is a DNA_part
        #direction is the direction we are looking in along the DNA_construct
        
        """the explorer sees a dna_part.
        1. it calculates effective direction of part relative to transcription direction
        2. it adds parts to growing transcripts
        3. it makes new transcripts and terminates growing transcripts.
        """
        logging.debug("looking at "+str(part))
        logging.debug("current_proteins are "+str(self.current_proteins))

        #The orientation of the part absolutely to the DNA_construct
        part_direction = part.direction
        #print("my part is ",part)
        #print("my direction is ",part.direction)

        #The relative orientation of the part to the reading direction
        if part_direction == reading_direction:
            effective_direction = "forward"
        else:
            effective_direction = "reverse"
        
        #count the number of RBSes that have been seen
        #currently this can only be 0 or 1 because RBS coupling is not implemented
        rbses_translating = len(self.current_proteins)

        #Case for Different Parts doing different things

        #Case 1: RBSes in the correct orientation
        if(isinstance(part,RBS) and effective_direction == "forward"):    
            if(rbses_translating>0):
                #rbs stops translation
                self.terminate_loop()
            #start a new protein
            self.current_proteins[part]= []

        #Case 2: CDSes in the correct orientation
        elif(rbses_translating >0 and isinstance(part,CDS) and effective_direction=="forward"):
            
            for rbs in self.current_proteins:
                #add to what is being translated
                self.current_proteins[rbs] += [[part,effective_direction]]
            #CDSes that aren't specifically designed to allow read-through
            #will terminate by default
            if(effective_direction not in part.no_stop_codons):
                #if the part has stop codons, terminate
                self.terminate_loop()

        #Case 3: all other parts usually terminate, unless they are designed to allow read-through
        elif(rbses_translating >0):
            if(effective_direction not in part.no_stop_codons):
                #if the part has stop codons, terminate
                self.terminate_loop()
                
    def terminate_loop(self):
        for rbs in self.current_proteins:
            self.all_proteins[rbs] = self.current_proteins[rbs]
        self.current_proteins = {} #reset who is translating

    def return_components(self,component):
        returnable_proteins = []
        for rbs in self.all_proteins:
            proteins = [a[0] for a in self.all_proteins[rbs]]
            #print("protein directions")
            #print([a.direction for a in proteins])
            returnable_proteins += proteins
            protein_species = [a.get_species() for a in proteins]
            if(protein_species==[]):
                protein_species = None
            component[rbs.position].protein = protein_species
        return returnable_proteins #these are CDS components
