################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 4/30/2020
#
#
################################################################

import matplotlib.pyplot as plt
import random
from warnings import warn
from matplotlib import cm
from .component import DNA
#import copy.copy as 

integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]
class DNA_part:
    #TODO can this be a component?
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
        self.color = None
        self.color2 = None
        if("color" in attributes):
            self.color = attributes["color"]
        if("color2" in attributes):
            self.color2 = attributes["color2"]
        if("direction" in attributes):
            self.direction = attributes["direction"]
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
        #TODO make this a copy function, so we don't have circular references
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
        return self




class IntegraseMechanism:
    def __init__(self,name,reactions={("attB","attP"):"attL",("attP","attB"):"attR"}):
        self.reactions = reactions
        self.attsites = []
        for reactants in reactions:
            self.attsites+=list(reactants)+list(reactions[reactants])
    def binds_to(self):
        return self.attsites
    def generate_products(self,site1,site2):
        """generates DNA_part objects corresponding to the products of recombination"""
        #the sites should have the same integrase and dinucleotide, otherwise it won't work
        assert(site1.integrase == site2.integrase)
        assert(site1.dinucleotide == site2.dinucleotide)
        integrase = site1.integrase
        dinucleotide = site1.dinucleotide
        rxn1 = (site1.part_type,site2.part_type)
        rxn2 = (site2.part_type,site1.part_type)
        try:
            prod1 = self.reactions[rxn1]
        except KeyError:
            #this means the reaction is not possible!!
            raise KeyError("{} not found to react with {} in {}".format(site1,site2,self.reactions))
        try:
            prod2 = self.reactions[rxn2]
        except KeyError:
            raise KeyError("{} not found to react with {} in {}".format(site2,site1,self.reactions))
        
        part_prod1 = DNA_part(prod1,attributes={"part_type":prod1,"dinucleotide":dinucleotide,
                                                    "integrase":integrase,"direction":site1.direction,
                                                    "color":site1.color,"color2":site2.color})
        part_prod2 = DNA_part(prod2,attributes={"part_type":prod2,"dinucleotide":dinucleotide,
                                                    "integrase":integrase,"direction":site2.direction,
                                                    "color":site2.color,"color2":site1.color})
        if(site1.direction=="forward"):
            return (part_prod1,part_prod2)
        else:
            part_prod2.direction = site1.direction
            part_prod1.direction = site2.direction
            return(part_prod2,part_prod1)

    def integrate(self,site1,site2):
        """perform an integration reaction between the chosen sites and alter the chassis"""
        #if one of the sites is not part of a construct then raise an error!
        if(not(isinstance(site1.parent_dna,DNA_construct))):
            raise ValueError("{} not part of a construct".format(site1))
        elif(not(isinstance(site2.parent_dna,DNA_construct))):
            raise ValueError("{} not part of a construct".format(site2))
        site1_initial_facing = dc(site1.direction) #keep track of what we reverse so we can flip it back!
        site2_initial_facing = dc(site2.direction)
        
        #if(site1.direction=="reverse"):
        #    reversed_site1 = True
        #if(site2.direction=="reverse"):

        #if(site1.direction=="reverse"):
        #    site1.parent_dna.reverse()
        #    reversed_site1 = True
        #reversed_site2 = False
        #if(site2.direction=="reverse"):
        #    site2.parent_dna.reverse()
        #    reversed_site2 = True
        cutpos1 = site1.pos
        cutpos2 = site2.pos
        if(site1.parent_dna==site2.parent_dna):
            #these sites are part of the same piece of DNA, so they are going to do an intramolecular reaction
            if(site1.pos > site2.pos):
                #we reverse which position is where
                cutpos2 = site1.pos
                cutpos1 = site2.pos
                prod1,prod2 = self.generate_products(site2,site1)
            else:
                prod1,prod2 = self.generate_products(site1,site2)
            
            dna = dc(site1.parent_dna)
            circularity = dna.circular
            
            if(site1.direction == site2.direction):
                #if the sites point in the same direction, then we are doing a deletion reaction
                cutdna_list = dna.parts_list[:cutpos1]+[prod1]+dna.parts_list[cutpos2+1:]
                cutdna = DNA_construct(cutdna_list,circular=circularity)
                if(site1_initial_facing != prod1.direction):
                    cutdna.reverse()
                newdna_list = [prod2]+dna.parts_list[1+cutpos1:cutpos2]
                newdna2 = DNA_construct(newdna_list,circular=True)
                if(site2_initial_facing != prod2.direction):
                    newdna2.reverse()
                newdna = [newdna2,cutdna]
            else:
                #this means we are dealing with an inversion
                inv_segment = dna.parts_list[cutpos1+1:cutpos2][::-1]
                inv_seg2 = []
                for ipart in inv_segment:
                    inv_seg2+=[ipart.reverse()]
                #print(inv_segment)
                #print(inv_seg2)
                #print([a.reverse() for a in dna.parts_list[cutpos1+1:cutpos2][::-1]])
                invertdna_list = dna.parts_list[:cutpos1]+ [prod1] + \
                                        inv_segment+ \
                                        [prod2]+dna.parts_list[cutpos2+1:]
                print(invertdna_list)
                inverteddna = DNA_construct(invertdna_list,circular=circularity)
                if(site1_initial_facing != prod1.direction):
                    inverteddna.reverse()
                newdna = [inverteddna]
        else:
            #otherwise these sites are on different pieces of DNA, so they are going to combine
            dna1 = dc(site1.parent_dna)
            dna2 = dc(site2.parent_dna)
            circ1 = dna1.circular
            circ2 = dna2.circular
            prod1,prod2 = self.generate_products(site1,site2)
            prod1.direction = "forward"
            prod2.direction = "forward"
            #now, cut the DNAs
            pcs1 = dna1.cut(site1.pos)
            print(pcs1)
            pcs2 = dna2.cut(site2.pos)
            #flip the DNAs around so that they are pointing forwards
            if(site1.direction=="reverse"):
                pcs1 = [a.reverse() for a in pcs1][::-1]
            if(site2.direction=="reverse"):
                pcs2 = [a.reverse() for a in pcs2][::-1]
            if(len(pcs1)== 1 and len(pcs2)==1):
                #in this case we are combining two circular plasmids to make one circular plasmid
                newdna = [DNA_construct(pcs1[0].parts_list+[prod1]+pcs2[0].parts_list+[prod2],circular=True)]
            elif(len(pcs1)==1 and len(pcs2)==2):
                #here we are combining a circular plasmid with a linear DNA
                newdna = [DNA_construct(pcs2[0].parts_list+[prod2]+pcs1[0].parts_list+[prod1]+pcs2[1].parts_list,circular=False)]
            elif(len(pcs1)==2 and len(pcs2)==1):
                #here we are combining a linear DNA with a circular plasmid
                newdna = [DNA_construct(pcs1[0].parts_list+[prod1]+pcs2[0].parts_list+[prod2]+pcs1[1].parts_list,circular=False)]
            elif(len(pcs1)==2 and len(pcs2)==2):
                #here we are combining two linear plasmids
                newdna_1 = DNA_construct(pcs1[0].parts_list+[prod1]+pcs2[1].parts_list,circular=False)
                newdna_2 = DNA_construct(pcs2[0].parts_list+[prod2]+pcs1[1].parts_list,circular=False)
                if(site1_initial_facing != prod1.direction):
                    newdna_1.reverse()
                if(site2_initial_facing != prod2.direction):
                    newdna_2.reverse()
                newdna = [newdna_1,newdna_2]
        return(newdna)   
        

from copy import deepcopy as dc
def rev_dir(dir):
    reversedict = {"forward":"reverse","reverse":"forward"}
    return reversedict[dir]




class DNA_construct(DNA):
    def __init__(self,
                parts_list,
                name=None,
                circular=False,
                mechanisms={},  # custom mechanisms
                parameters={},  # customized parameters
                attributes=[],
                initial_conc=None,
                parameter_warnings = True,
                **keywords):
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
                if(part.direction==None):
                    #in this case the direction wasn't provided, but we assume it's forward!
                    myparts += [dc(part).clone(curpos,"forward",self)]
                else:
                    #in this case the direction is provided in the part. This is a "recloned" part
                    #TODO make sure it doesnt screw up the parts
                    myparts += [dc(part).clone(curpos,part.direction,self)]
            curpos += 1
        self.parts_list = myparts
        cmap = cm.Set1(range(len(self.parts_list)+5))
        pind = 0
        for part in self.parts_list:
            if(type(part.color)==type(None)):
                #if the color isn't set, let's set it now!
                part.color = cmap[pind][:-1]
            if(type(part.color2)==type(None)):
                if(part.part_type in ["attL","attR"]):
                    #this is the only scenario in which we need a color2
                    #TODO make this more general
                    pind+=1
                    part.color2 = cmap[pind][:-1]
            pind+=1
        self.circular=circular
        if(name == None):
            name = str(self)
        DNA.__init__(self=self,name=name,length = len(parts_list),
                    mechanisms=mechanisms,parameters=parameters,
                    attributes=attributes,initial_conc = initial_conc,
                    parameter_warnings=parameter_warnings, **keywords)
    def reverse(self):
        """reverses everything, without actually changing the DNA"""
        newlist = self.parts_list[::-1]
        newpos = 0
        for part in newlist:
            part.pos = newpos
            part.direction = rev_dir(part.direction)
            newpos+=1
        self.parts_list = newlist
        return self
    def explore_txtl(self):
        """this function finds promoters and terminators and stuff in the construct"""
        # lets try to make this more modular shall we?\
        proteins = {}
        rnas = {}
        for direction in ["forward","reverse"]:
            explorer = TxTl_Explorer()
            explorer.direction=direction
            if(direction == "reverse"):
                #if we go backwards then also list the parts backwards
                #deepcopy so we don't mangle the list
                newlist = dc(self.parts_list)[::-1]
            else:
                newlist = dc(self.parts_list)
            #explorer.make_rna(self.name)
            part_index = 0
            keep_going = 1
            second_looping = 0
            while keep_going:
                part = newlist[part_index]
                keep_going = explorer.see(part)
                part_index+=1
                if(part_index==len(self.parts_list)):
                    part_index = 0
                    if(self.circular and second_looping == 0):
                        second_looping = 1
                        explorer.second_loop()
                    else:
                        explorer.end()
                        break
            proteins.update(explorer.get_proteins())
            rnas.update(explorer.get_rnas())
        return rnas,proteins
    def __repr__(self):
        """the name of a DNA has to be unique and usable as a component name in biocrnpyler"""
        output = ""
        output = "dna_"+'_'.join([str(a) for a in self.parts_list])
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
        """checks if this construct contains a certain part, or a copy of a certain part"""
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
    
    def cut(self,position,keep_part=False):
        """cuts the construct and returns the resulting piece or pieces"""
        left = self.parts_list[:position]
        if(keep_part):
            rightpos = position
        else:
            rightpos = position+1
        right = self.parts_list[rightpos:]
        newDNA = []
        if(self.circular):
            #this means we return one linear piece
            newDNA += [DNA_construct(right+left,circular=False)]
        else:
            #this means we return two linear pieces
            newDNA += [DNA_construct(left,circular=False),DNA_construct(right,circular=False)]
        return newDNA
    def append(self,part,direction="forward"):
        """stick a part on at the end"""
        self.parts_list+=[part]
        part.clone(len(self.parts_list)-1,direction,self)
    def prepend(self,part,direction="forward"):
        """stick a part on at the end"""
        self.parts_list=[part]+self.parts_list
        part.clone(0,direction,self)
        for part_ind,part in zip(range(len(self.parts_list)),self.parts_list):
            part.pos = part_ind
    def update_species(self):
        #TODO should we have this get all the species that our DNA is responsible for?
        species = [self.get_species()]
        rnas,proteins = self.explore_txtl()
        components = self.update_components(rnas,proteins)

        for component in components:
            species += component.update_species()
        return species

class TxTl_Explorer():
    def __init__(self,possible_rxns=("transcription","translation"),\
                        part_defs={"promoter":"transcription", "rbs":"translation"},direction="forward"):
        self.current_rnas = {}
        self.current_proteins = {}
        self.made_rnas = {}
        self.all_rnas = {}
        self.made_proteins = {}
        self.current_rbs = None
        self.possible_rxns = possible_rxns
        self.part_defs = part_defs
        self.direction=direction
        self.second_looping = False
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
            terminated = 0
            if(part.part_type in self.part_defs and \
                self.part_defs[part.part_type]=="translation" and \
                            effective_direction == "forward"):
                self.current_proteins[part]= []
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
                self.current_proteins[self.current_rbs] = self.current_proteins[self.current_rbs]+\
                                                                        [[part,effective_direction]]
            elif(self.current_rbs != None):
                #this means we CAN'T translate through the current part, but we are translating
                self.current_proteins[self.current_rbs] = self.current_proteins[self.current_rbs] + \
                                                    [[part,effective_direction]]
                self.current_rbs = None #current RBS has done its work
            for promoter in list(self.current_rnas.keys()):
                self.current_rnas[promoter]+=[[part,effective_direction]]
                #we add the current part
                if(part.part_type=="terminator" and effective_direction == "forward"):
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
        #print(part.part_type)
        #if(part.part_type in self.part_defs):
            #print(self.part_defs[part.part_type])
        if((effective_direction=="forward") and \
                (part.part_type in self.part_defs) and \
                (self.part_defs[part.part_type] in self.possible_rxns) and \
                (self.part_defs[part.part_type] == "transcription")):
            #this if statement makes sure that the current part wants to transcribe, and
            #that is something that we are allowed to do
            self.current_rnas.update({part:[]})
            return 1
        elif((("transcription" not in self.possible_rxns) or self.second_looping) and \
                                         len(self.current_rnas.keys())==0):
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
            self.current_rnas.update({promname:[]})
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
        rna_partslist = self.current_rnas[promoter]
        rna_construct = RNA_construct(rna_partslist)
        del self.current_rnas[promoter] # this removes the current RNA from the list, because it's
                                        # being terminated!!
        #current_rna_name = str(promoter)+"-"+str(rna_construct)
        self.made_proteins[rna_construct]={}
        #compile the name, keeping track of which promoter made the rna and all the parts on it
        for rbs in self.current_proteins:
            #ending the RNA also ends all the proteins being generated here.
            proteins_per_rbs = []
            # TODO this will run multiple times for each RNA. make it so it runs once!
            for protein_part in self.current_proteins[rbs]:
                if(protein_part[1] == "forward"):
                    #it is essential to be forwards
                    if(protein_part[0].protein != None):
                        #this happens if we do in fact make a protein
                        proteins_per_rbs+=[protein_part[0]]
            self.made_proteins[rna_construct].update({rbs:proteins_per_rbs})
            # so this has to be done at the end of the RNA only because only then do we know
            # exactly what that full RNA is going to be.
            # we want to clear the current proteins, but there could be multiple RNAs that run over the same proteins,
            # so wait until we are done tallying all the RNAs before removing everything from current_proteins
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


class RNA_construct(DNA_construct):
    #TODO this must work as an RNA! How is that different from DNA?
    # 1: there are no promoters to detect
    # 2: it must be linear
    # 3: parts know that they are part of an RNA
    # 4: the representation is different
    # 5: maybe a DNA object contains RNAs??
    def __init__(self,parts_list,name=None,**keywords):
        
        DNA_construct.__init__(self=self,parts_list=parts_list,circular=False,name=name,**keywords)
        if(name == None):
            name = str(self)
        self.name = name
        self.material_type = "rna"
    def explore_txtl(self):
        """an RNA has no tx, only TL! central dogma exists, right?"""
        # lets try to make this more modular shall we?
        explorer = TxTl_Explorer(possible_rxns = ("translation",))
        explorer.make_rna(self.name)
        part_index = 0
        keep_going = 1
        while keep_going:
            part = self.parts_list[part_index]
            keep_going = explorer.see(part)
            part_index+=1
            if(part_index==len(self.parts_list)):
                part_index = 0
                if(self.circular):
                    explorer.second_loop()
                else:
                    explorer.end()
                    break
        proteins = explorer.get_proteins()
        return proteins[list(proteins.keys())[0]]
    def __repr__(self):
        """the name of an RNA should be different from DNA, right?"""
        output = ""
        output = 'rna_'+'_'.join([str(a) for a in self.parts_list])
        return output

class Chassis():
    def __init__(self,dna_constructs,name="myChassis"):
        """a chassis is a combination of DNA parts. This is either your TXTL tube which can contain
        multiple different DNA molecules, or your cell, which also contains multiple different
        DNA molecules. This is useful because interactions between different DNA
        molecules that are part of the same test tube can occur"""
        self.dna_constructs = dna_constructs
        self.name = name
    def extract_constructs_containing(self,part):
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
            if(integrase in intnames or type(intnames)==type(None)):
                #this only reacts with integrases that we said exist, or all of them if we didn't
                #say anyhting
                attsites = int_dict[integrase]
                #but now we need to know what kind of integrase reactions are possible

