################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 4/30/2020
#
#
################################################################

from matplotlib import cm
import matplotlib.pyplot as plt
import random
class DNA_part:
    def __init__(self,name,biocrnpyler_class=None,attributes={},**keywords):
        """this represents a sequence. These get compiled into working components"""
        self.biocrnpyler_class = biocrnpyler_class
        self.sequence = None
        self.integrase = None
        self.dinucleotide = None
        self.part_type=None
        self.regulator = None
        self.name = name
        self.protein = None
        self.no_stop_codons = []
        self.attributes = attributes


        if("no_stop_codons" in attributes):
            self.no_stop_codons = attributes["no_stop_codons"]
        if("sequence" in attributes):
            self.sequence = attributes["sequence"]
        if("part_type" in attributes):
            #acceptable part 
            self.part_type = attributes["part_type"]
        if(self.part_type in ["attB","attP","attL","attR"]):
            if("integrase" in attributes):
                self.integrase = attributes["integrase"]
            else:
                self.integrase = "int1"
            if("dinucleotide" in attributes):
                self.dinucleotide = attributes["dinucleotide"]
            else:
                self.dinucleotide = 1
        if(self.part_type == "promoter"):
            if("regulator" in attributes):
                self.regulator = attributes["regulator"]
        if(self.part_type == "CDS"):
            if("protein" in attributes):
                self.protein = attributes["protein"]
    def __repr__(self):
        myname = self.name
        if(self.part_type in ["attB","attP","attL","attR"]):
            myname += "-" + self.integrase
            if(self.dinucleotide != 1):
                myname += "-"+str(self.dinucleotide) 
        return myname
    def dnaplotlib_design(self,direction=True,color=(1,4,2),color2=(3,2,4)):
        dnaplotlib_dict = {\
            "promoter":"Promoter",\
            "rbs":"RBS",\
            "CDS":"CDS",\
            "terminator":"Terminator",\
            "attP":"RecombinaseSite",\
            "attB":"RecombinaseSite",\
            "attL":"RecombinaseSite2",\
            "attR":"RecombinaseSite2"}
        dpl_type = dnaplotlib_dict[self.part_type]
        outdesign = [{'type':dpl_type,"name":self.name,"fwd":direction,'opts':{'color':color,'color2':color2}}]
        if(self.regulator!= None):
            outdesign += [{"type":"Operator","name":self.regulator,"fwd":direction,'opts':{'color':color,'color2':color2}}]
        if(not direction):
            outdesign = outdesign[::-1]
        return outdesign
from copy import deepcopy as dc
def rev_dir(dir):
    reversedict = {"forward":"reverse","reverse":"forward"}
    return reversedict[dir]
class DNA_construct:
    def __init__(self,parts_list,circular=False,**keywords):
        """this represents a bunch of parts in a row.
        A parts list has [[part,direction],[part,direction],...]"""
        self.parts_list = parts_list
        self.circular=circular
    def explore_txtl(self):
        """this function finds promoters and terminators and stuff in the construct"""
        transcribing = False
        translating = False
        txrxns = []
        tlrxns = []
        made_proteins = {}
        all_rnas = []
        all_proteins = []
        current_rnas = {}
        current_proteins = {}
        current_rbs = None
        for direction in ["forward","reverse"]:
            #go forwards and also backwards!
            if(direction == "reverse"):
                #if we go backwards then also list the parts backwards
                #deepcopy so we don't mangle the list
                newlist = dc(self.parts_list)[::-1]
            else:
                newlist = dc(self.parts_list)
            for part_index in range(len(newlist)):
                #
                # iterate through all parts!
                #
                part = newlist[part_index]
                '''
                #this part is for picking out the next and previous parts.
                #is it needed? unclear
                try:
                    #look ahead!
                    next_part = newlist[part_index+1]
                except IndexError:
                    #this means we fell off the end of the list!
                    if(self.circular == True):
                        #if we are circular then look back at the beginning
                        next_part = newlist[0]
                    else:
                        #otherwise, it's the end!
                        next_part = None
                previous_part = newlist[part_index-1]
                #unlike looking off the end of the list, -1 is a valid index.
                if(self.circular == False):
                    previous_part = None
                #'''
                effective_direction = part[1]

                if(direction=="reverse"):
                    #if we are reverse then everything is backwards
                    effective_direction = rev_dir(effective_direction)
                part_object = part[0]


                if(current_rnas!={}):
                    #this means we are compiling an RNA
                    terminated = 0
                    if(part_object.part_type=="rbs" and effective_direction == "forward"):
                        #rbs makes protein!! we add this to a list and then LATER decide which RNA gets it.
                        #because, there could be no protein after the rbs, and then we dont care!
                        #we don't track ORFs or do any kind of ribosome simulation so this could get tricky
                        #ideally you just look at one protein following the RBS, but...
                        #1- conditional RNA motifs could exist
                        #2- attachment sites inside of protein coding sequences
                        #3- fusion proteins
                        #4- translational coupling
                        #
                        #i think the solution is just to track cds'es which have stop codons.
                        #let's just assume that valid CDSes have stop codons.
                        #logically though 
                        current_proteins[part_object]= []
                        if(current_rbs == None):
                            #this means nobody is translating yet
                            current_rbs = part_object
                        else:
                            #another RBS is translating. this could mean a few things:
                            #1: translational coupling. the previous translation ends and turns
                            #   into this one, with an added coupling term when that gets implemented
                            #2: it just terminates the previous translation. Chances are it would, right?
                            current_rbs = part_object
                            # TODO add coupling
                    elif(effective_direction in part_object.no_stop_codons and current_rbs != None):
                        #this means we can translate through the current part!
                        current_proteins[current_rbs] = current_proteins[current_rbs]+[[part_object,effective_direction]]
                    elif(current_rbs != None):
                        #this means we CAN'T translate through the current part, but we are translating
                        current_proteins[current_rbs] = current_proteins[current_rbs] + [[part_object,effective_direction]]
                        #this also means we perhaps have a chance to evaluate the latest 'protein' and see if it's good.
                        current_rbs = None
                    
                    for promoter in current_rnas:
                        current_rnas[promoter]+=[[part_object,effective_direction]]
                        #we add the current part
                        if(part_object.part_type=="terminator" and effective_direction == "forward"):
                            terminated = 1 #a terminator ends everything. This does not account for leakiness
                            #this stops the RNA! but only in the forward direction
                            current_rna_name = str(promoter)+"="+'_'.join([str(a[0])+"-r"*(a[1]=="reverse") for a in current_rnas[promoter]])
                            made_proteins[current_rna_name]={}
                            #compile the name, keeping track of which promoter made the rna and all the parts on it
                            for rbs in current_proteins:
                                #ending the RNA also ends all the proteins being generated here.
                                proteins_per_rbs = []
                                # TODO this will run multiple times for each RNA. make it so it runs once!
                                for protein_part in current_proteins[rbs]:
                                    if(protein_part[1] == "forward"):
                                        #it is essential to be forwards
                                        if(protein_part[0].protein != None):
                                            #this happens if we do in fact make a protein
                                            proteins_per_rbs+=[protein_part[0]]
                                made_proteins[current_rna_name].update({rbs:proteins_per_rbs})
                                #we want to clear the current proteins, but there could be multiple RNAs that run over the same proteins,
                                #so wait until we are done tallying all the RNAs before removing everything from current_proteins
                            all_rnas +=[[promoter,current_rnas[promoter]]]
                    
                        
                    if(terminated):
                        current_rnas = {}
                        current_proteins = {}
                        current_rbs = None
                        terminated = 0
                if(part_object.part_type=="promoter" and effective_direction=="forward"):
                    #this part is a promoter, so it transcribes!
                    current_rnas.update({part_object:[]})
        return all_rnas,made_proteins
    def __repr__(self):
        output = ""
        output = '_'.join([str(a[0])+"-r"*(a[1]=="reverse") for a in self.parts_list])
        if(self.circular):
            output+="-o"
        return output
    def dnaplotlib_design(self):
        outdesign = []
        cmap = cm.Set1(range(len(self.parts_list)))
        pind = 0
        for part in self.parts_list:

            outdesign+=part[0].dnaplotlib_design(part[1]=="forward",color=cmap[pind][:-1],color2 = random.choice(cmap)[:-1])
            pind+=1
        return outdesign