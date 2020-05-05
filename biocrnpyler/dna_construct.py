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
from warnings import warn
integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]
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
        self.pos = None
        self.parent_dna = None
        self.direction=None

        if("no_stop_codons" in attributes):
            self.no_stop_codons = attributes["no_stop_codons"]
        if("sequence" in attributes):
            self.sequence = attributes["sequence"]
        if("part_type" in attributes):
            #acceptable part 
            self.part_type = attributes["part_type"]
        if(self.part_type in integrase_sites):
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
        if(self.part_type in integrase_sites):
            myname += "-" + self.integrase
            if(self.dinucleotide != 1):
                myname += "-"+str(self.dinucleotide) 
        if(self.direction =="reverse"):
            myname += "-r"
        return myname
    def clone(self,position,direction,parent_dna):
        """this defines where the part is in what piece of DNA"""
        self.pos = position
        self.direction = direction
        self.parent_dna = parent_dna
        return self
    def unclone(self):
        """removes the current part from anything"""
        self.pos = None
        self.direction = None
        self.parent_dna = None
        return self
    def reverse(self):
        if(self.direction=="forward"):
            self.direction = "reverse"
        elif(self.direction=="reverse"):
            self.direction = "forward"
        elif(self.direction==None):
            warn(self.name+" has no direction. Perhaps it wasn't cloned?")
        else:
            raise ValueError("direction is not forward or reverse! It's "+self.direction)
    def dnaplotlib_design(self,direction=None,color=(1,4,2),color2=(3,2,4)):
        if(direction==None and self.direction != None):
            direction == self.direction=="forward"
        elif(direction==None):
            direction = True
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

class IntegraseMechanism:
    def __init__(self,name,reactions={("attB","attP"):("attL","attR")}):
        self.reactions = reactions
        self.attsites = []
        for reactants in reactions:
            self.attsites+=list(reactants)+list(reactions[reactants])
    def binds_to(self):
        return self.attsites
    def integrate(site1,site2,chassis):
        """perform an integration reaction between the chosen sites and alter the chassis"""
        reactants1 = chassis.extract_construct(site1)

        reactants2 = chassis.extract_construct(site2)

        for reactant1 in reactants1:

            for reactant2 in reactants2:
                pass

        

from copy import deepcopy as dc
def rev_dir(dir):
    reversedict = {"forward":"reverse","reverse":"forward"}
    return reversedict[dir]
class DNA_construct:
    def __init__(self,parts_list,circular=False,**keywords):
        """this represents a bunch of parts in a row.
        A parts list has [[part,direction],[part,direction],...]"""
        myparts = []
        curpos = 0
        for part in parts_list:
            if(isinstance(part,list)):
                #if you give it a list, then that means the second element of the
                #list is the direction
                myparts += [dc(part[0]).clone(curpos,part[1],self)]
            elif(isinstance(part,DNA_part)):
                #in this case the direction wasn't provided, but we assume it's forward!
                myparts += [dc(part[0]).clone(curpos,"forward",self)]
            curpos += 1
        self.parts_list = myparts
        self.circular=circular
    def explore_txtl(self):
        """this function finds promoters and terminators and stuff in the construct"""
        made_proteins = {}
        all_rnas = {}
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
                effective_direction = part.direction

                if(direction=="reverse"):
                    #if we are reverse then everything is backwards
                    effective_direction = rev_dir(effective_direction)
                part_object = part


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
                            rna_partslist = current_rnas[promoter]
                            rna_construct = DNA_construct(rna_partslist,circular=False)
                            current_rna_name = str(promoter)+"="+str(rna_construct)
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
                            all_rnas.update({promoter:current_rnas[promoter]})
                    
                        
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
        output = '_'.join([str(a) for a in self.parts_list])
        if(self.circular):
            output+="-o"
        return output
    def list_integrase(self,int_dict={}):
        """lists all the parts that can be acted on by integrases"""
        for part in self.parts_list:
            if(part.integrase != None):
                if(part.integrase in int_dict):
                    int_dict.update({part.integrase:int_dict[part.integrase]+[part]})
                else:
                    int_dict[part.integrase]=[part]
        return int_dict
                    
    def __contains__(self,obj2):
        """we should be able to check if this construct has a certain part"""
        if(isinstance(obj2,DNA_part)):
            #if we got a DNA part it could mean one of two things:
            #1 we want to know if a dna part is anywhere
            #2 we want to know if a specific DNA part is in here
            #this is complicated by the fact that we want to have the same DNA part be reusable
            #in many locations
            if(obj2.parent_dna==self):
                #the object should already know if it's a part of me
                return True
            elif(obj2.parent_dna==None):
                #this object has been orphaned. 
                #that means we are looking for matching objects in any position
                new_obj2 = dc(obj2).unclone()
                uncloned_list = [a.unclone() for a in dc(self.parts_list)]
                return new_obj2 in uncloned_list
            else:
                return False
        elif(isinstance(obj2,str)):
            #if we get a string, that means we want to know if the name exists anywhere
            return obj2 in str(self)

    def dnaplotlib_design(self):
        outdesign = []
        cmap = cm.Set1(range(len(self.parts_list)))
        pind = 0
        for part in self.parts_list:
            outdesign+=part.dnaplotlib_design(part.direction=="forward",color=cmap[pind][:-1],color2 = random.choice(cmap)[:-1])
            pind+=1
        return outdesign
class Chassis():
    def __init__(self,dna_constructs,name="myChassis"):
        """a chassis is a combination of DNA parts. This is either your TXTL tube which can contain
        multiple different DNA molecules, or your cell, which also contains multiple different
        DNA molecules. This is useful because interactions between different DNA
        molecules that are part of the same test tube can occur"""
        self.dna_constructs = dna_constructs
        self.name = name
    def extract_construct(self,part):
        """extract a DNA_construct which contains part.
        part can be a DNA_part or a string which matches part of
        the DNA_construct's name"""
        conlist = []
        for construct in self.dna_constructs:
            if(part in construct):
                conlist+= [construct]
        return conlist

    def explore_integrases(self,intnames=None):
        """this explores all the possible integrase-motivated DNA configurations. If some
        integrases aren't present, then define intnames to be a list of names of the
        integrases which are present.

        An integrase can act in different ways. 
        * serine integrases recombine B and P sites that turn into L and R sites, 
                    and only sites with the same dinucleotide can be recombined.
        * serine integrases with directionality factors recombine L and R sites
                    with the same dinucleotide
        * Invertases only do flipping reactions
        * resolvases only do deletion reactions
        * FLP or CRE react with homotypic sites, so site1+site1 = site1+site1. But
                    there are still different types of sites which are orthogonal. For
                    example, a CRE type 1 or a CRE type 2 site. The sites can also be palindromic,
                    which means that they can react in either direction.
        """
        int_dict = {}
        for construct in self.dna_constructs:
            #list each integrase that exists and which sites they react with
            int_dict = construct.list_integrase(int_dict)
        for integrase in int_dict:
            #now, going through each one, generate the reactions and species that arise
            if(integrase in intnames or intnames==None):
                #this only reacts with integrases that we said exist, or all of them if we didn't
                #say anyhting
                attsites = int_dict[integrase]
                #but now we need to know what kind of integrase reactions are possible

