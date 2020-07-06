################################################################
#       DNA_construct: a higher level for construct compilation
#       Author: Andrey Shur
#       Latest update: 6/4/2020
#
#
################################################################

import matplotlib.pyplot as plt
import random

from .chemical_reaction_network import ComplexSpecies, Species, OrderedComplexSpecies
from .mechanisms_binding import One_Step_Cooperative_Binding, Combinatorial_Cooperative_Binding
from matplotlib import cm
from .components_basic import DNA, Protein, RNA

from .dna_part import DNA_part
from .dna_part_misc import AttachmentSite
from .dna_part_rbs import RBS
from .dna_part_cds import CDS
from .dna_part_terminator import Terminator
from .dna_part_promoter import Promoter
#from .dna_part_promoter import 
import itertools as it
import copy

from warnings import warn

#integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]

def all_comb(input_list):
    out_list = []
    for i in range(1,len(input_list)+1):
        out_list += it.combinations(input_list,i)
    return out_list


def rev_dir(dir):
    reversedict = {"forward":"reverse","reverse":"forward"}
    return reversedict[dir]

class dummyAssembly:
    """this class exists because of how promoters work. Is it necessary? I don't know"""
    def __init__(self,species):
        self._dna = copy.deepcopy(species)
    def get_dna(self):
        return copy.deepcopy(self._dna)
    def set_dna(self,species):
        self._dna = copy.deepcopy(species)
    def del_dna(self):
        del self._dna
    dna = property(get_dna,set_dna,del_dna)

        
        

class DNA_construct(DNA):
    def __init__(self,
                parts_list,
                name=None,
                circular=False,
                mechanisms={},  # custom mechanisms
                parameters={},  # customized parameters
                attributes=[],
                initial_conc=None,
                parameter_warnings = True, copy_parts=True,
                **keywords):
        """this represents a bunch of parts in a row.
        A parts list has [[part,direction],[part,direction],...]"""
        myparts = []
        curpos = 0
        if(copy_parts):
            parts_list = [[copy.deepcopy(a[0]).unclone(),a[1]] for a in parts_list]
        for part in parts_list:
            if(isinstance(part,list)):
                #if you give it a list, then that means the second element of the
                #list is the direction
                myparts += [part[0].clone(curpos,part[1],self)]
            elif(isinstance(part,DNA_part)):
                if(part.direction==None):
                    #in this case the direction wasn't provided, but we assume it's forward!
                    myparts += [part.clone(curpos,"forward",self)]
                else:
                    #in this case the direction is provided in the part. This is a "recloned" part
                    myparts += [part.clone(curpos,part.direction,self)]
            curpos += 1
        self.parts_list = myparts
        cmap = cm.Set1(range(len(self.parts_list)+5))
        pind = 0
        for part in self.parts_list:
            if(type(part.color)==type(None)):
                #if the color isn't set, let's set it now!
                part.color = cmap[pind][:-1]
            if(type(part.color2)==type(None)):
                if(isinstance(part,AttachmentSite) and part.site_type in ["attL","attR"]):
                    #this is the only scenario in which we need a color2
                    #TODO make this more general
                    pind+=1
                    part.color2 = cmap[pind][:-1]
            pind+=1
        self.circular=circular
        if(name is None):
            name = self.make_name() #automatic naming
        self.update_dna(name) #this creates the dna Species object
        DNA.__init__(self=self,name=name,length = len(parts_list),
                    mechanisms=mechanisms,parameters=parameters,
                    attributes=attributes,initial_conc = initial_conc,
                    parameter_warnings=parameter_warnings, **keywords)
        self.transcripts = []
        self.set_parameter_warnings(parameter_warnings)
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
    def update_dna(self, dna, attributes = None):
        if dna is None:
            self.dna = self.set_species(self.name, material_type = "dna", attributes = attributes)
        else:
            self.dna = self.set_species(dna, material_type = "dna", attributes = attributes)
    def set_parameter_warnings(self, parameter_warnings):
        """updates all components with the proper parameter warnings"""
        self.parameter_warnings = parameter_warnings
        if(self.parameter_warnings is not None):
            for part in self.parts_list:
                part.set_parameter_warnings(parameter_warnings)
    def update_parameters(self, mixture_parameters = {}, parameters = {},
                          overwrite_custom_parameters = True):
        """update parameters of all parts in the construct"""
        DNA.update_parameters(self = self,
                              mixture_parameters = mixture_parameters,
                              parameters = parameters)
        for part in self.parts_list:
            part.update_parameters(mixture_parameters = mixture_parameters,
                                    parameters = parameters,
                                    overwrite_custom_parameters = False)
    def update_mechanisms(self, mixture_mechanisms = {}, mechanisms = {},
                          overwrite_custom_parameters = True,overwrite_custom_mechanisms = False):
        DNA.update_mechanisms(self = self,
                              mixture_mechanisms = mixture_mechanisms,
                              mechanisms = mechanisms,overwrite_custom_mechanisms = overwrite_custom_mechanisms)
        for part in self.parts_list:
            part.update_mechanisms(mechanisms = mechanisms,
                                    mixture_mechanisms = mixture_mechanisms,
                                    overwrite_custom_mechanisms = overwrite_custom_mechanisms)
    def explore_txtl(self):
        """this function finds promoters and terminators and stuff in the construct"""
        # lets try to make this more modular shall we?\
        proteins = {}
        rnas = {}
        for direction in ["forward","reverse"]:
            explorer = TxTl_Explorer(parameter_warnings=self.parameter_warnings)
            explorer.direction=direction
            if(direction == "reverse"):
                #if we go backwards then also list the parts backwards
                #deepcopy so we don't mangle the list
                newlist = copy.deepcopy(self.parts_list)[::-1]
            else:
                newlist = copy.deepcopy(self.parts_list)
            part_index = 0
            keep_going = 1
            second_looping = 0
            while keep_going: #we may have to loop through the parts twice because the construct is circular.
                                #this does not account for infinite loops, RCA, etc.
                part = newlist[part_index] #pull out the part
                #print("part is " + str(part))
                keep_going = explorer.see(part) #the explorer keeps track of everything
                part_index+=1 #keep going+
                if(part_index==len(self.parts_list)):
                    part_index = 0 #this takes care of looping
                    if(self.circular and second_looping == 0):
                        #print('restart from beginning')
                        second_looping = 1 #during the second loop, we don't care about promoters
                        explorer.second_loop() #tell the explorer that we don't care about promoters
                    else:
                        explorer.end() #this completes all "in progress" RNAs and proteins
                        break
            #proteins.update(explorer.get_proteins())
            rnas.update(explorer.get_rnas())
        proteins = {}
        for promoter in rnas:
            _,prots = rnas[promoter].explore_txtl()
            proteins.update(prots)
        return rnas,proteins
    def __repr__(self):
        """this is just for display purposes"""
        txt = "dna = "+ self.name
        rnas,proteins = self.explore_txtl()
        if(len(rnas)>0):
            for promoter in rnas:
                txt += "\n\t"+repr(promoter)
                txt += "\n\t" + repr(rnas[promoter])
                if(rnas[promoter] in proteins):
                    for rbs in proteins[rnas[promoter]]:
                        txt += "\n\t"+repr(rbs)
                        #print("promoter is "+str(promoter))
                        #print("rbs is "+ str(rbs))
                        #print("rnas is "+ str(rnas))
                        #print("proteins is "+ str(proteins))
                        for protein in proteins[rnas[promoter]][rbs]:
                            txt += "\n\tprotein = "+repr(protein)
        return txt
    
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
                new_obj2 = copy.deepcopy(obj2).unclone()
                uncloned_list = [a.unclone() for a in copy.deepcopy(self.parts_list)]
                return new_obj2 in uncloned_list
            else:
                return False
        elif(isinstance(obj2,str)):
            #if we get a string, that means we want to know if the name exists anywhere
            return obj2 in str(self)
    def get_species(self):
        ocomplx = [Species(a.name,material_type="dna") for a in self.parts_list]+[Species(self.name,material_type="dna")]
        return OrderedComplexSpecies(ocomplx,name=self.name,material_type="dna")
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
        """stick a part on at the beginning"""
        self.parts_list=[part]+self.parts_list
        part.clone(0,direction,self)
        for part_ind,part in zip(range(len(self.parts_list)),self.parts_list):
            part.pos = part_ind
    def update_attachment_sites(self,int_sites,combinatorial_species=None):
        """this function updates integrase sites ONLY WITH integrase binding"""
        if(combinatorial_species is None):
            combinatorial_species = [self.get_species()]
        sites = []
        for comb_specie in combinatorial_species:
                if(isinstance(comb_specie,OrderedComplexSpecies) and comb_specie.species[-1].name == self.name):
                    for site in int_sites:
                        if(isinstance(comb_specie.species[site.pos],ComplexSpecies)):
                            #if we already made a complex, then don't do anything
                            continue
                        new_site = copy.deepcopy(site)
                        new_mv_self = copy.copy(comb_specie)
                        new_mv_self.attributes = ["bindloc_"+str(site.pos)]
                        new_assy = dummyAssembly(new_mv_self)
                        new_site.assembly = new_assy
                        sites += [new_site]
        return sites
                        
    def update_promoters(self,rnas,proteins,combinatorial_species=None):
        """this function properly populates all promoters
        returns a combinatorial list of transcribable states"""
        if(combinatorial_species is None):
            combinatorial_species = [self.get_species()]
        promoters = []
        for comb_specie in combinatorial_species:
            if(isinstance(comb_specie,OrderedComplexSpecies) and comb_specie.species[-1].name == self.name):
                for promoter in rnas:
                    #print(promoter)
                    if(promoter == "unknown"):
                        continue
                    new_promoter = copy.deepcopy(promoter) #copy the promoter because we could have the same one many times
                    if(isinstance(comb_specie.species[promoter.pos],ComplexSpecies)):
                        #for each promoter, we generate all complexes which
                        #result in binding at a specific site. If there's something already there, then
                        #we skip those
                        continue
                    new_mv_self = copy.deepcopy(comb_specie)
                    #TODO make this attribute thing better
                    new_mv_self.attributes = ["bindloc_"+str(promoter.pos)]
                    new_assy = dummyAssembly(new_mv_self) #this allows the promoter to get the right DNA
                    #a promoter can make multiple RNAs, however...
                    new_promoter.transcript = rnas[promoter].get_species() 
                    rbslist = proteins[rnas[promoter]]
                    new_promoter.protein = [rbslist[a] for a in rbslist]
                    #currently, it's only possible for one promoter to make one RNA
                    new_promoter.assembly = new_assy 
                    promoters += [new_promoter]
        
        return promoters
    def update_rbses(self,proteins,combinatorial_species=None):
        my_rbses = []
        for rna in proteins:
            my_rbses += rna.update_rbses()
        return my_rbses
    def remove_bindloc(self,spec_list):
        """go through every species on a list and remove any "bindloc" attributes"""
        spec_list2 = copy.copy(spec_list)
        out_sp_list = []
        for specie in spec_list2:
            #go through the species and remove the "bindloc" attribute
            #I don't care about the binding now that I am done generating species
            new_attrib = []
            for attrib in specie.attributes:
                if("bindloc" in attrib):
                    pass
                else:
                    new_attrib+=[attrib]
            specie.attributes = new_attrib
            out_sp_list += [specie]
        return out_sp_list
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
        spec_indexes = [a[1] for a in spec_list] #extract all indexes
        compacted_indexes = sorted(list(set(spec_indexes)))
        prototype_list = [None]*len(compacted_indexes)
        for spec in spec_list:
            #go through every element and put it in the right place
            proto_ind = compacted_indexes.index(spec[1]) #where to put it?
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
        #print("comb_list is {}".format(str(comb_list)))
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
    def update_combinatorial_complexes(self,active_components):
        species = [self.get_species()]
        for part in active_components:
            #first we make binary complexes
            sp_list =  self.remove_bindloc(part.update_species())
            species+=sp_list
        unique_complexes = {}
        possible_backbones = {self.name:self.get_species()}
        #possible_backbones = {a.name:a.get_species() for a in [self]+[a for a in proteins]}
        #species need to be uniqueified
        unique_species = list(set(species)) 
        for specie in unique_species:
            #in this list we extract all the variants of the complexes from possible_backbones
            #that exist in our species list.
            if(isinstance(specie,OrderedComplexSpecies) and specie.species[-1].name in possible_backbones):
                #we only care about OrderedComplexSpecies made from this construct
                if(specie.species[-1].name in unique_complexes):
                    unique_complexes[specie.species[-1].name] += [specie]
                else:
                    unique_complexes[specie.species[-1].name] = [specie]
        #unique_complexes now has a list of all the non-combinatorial complexes we can make
        combinatorial_complexes = []
        for bb_name in unique_complexes:
            #for each backbone, make all combinatorial combinations.
            comp_binders = []
            for unit in unique_complexes[bb_name]:
                #for each complex for this backbone, find out what is bound at which location
                comp_bound = None
                pos_i = 0
                for pos in unit.species:
                    #find the position that has a ComplexSpecies in it
                    if(isinstance(pos,ComplexSpecies) and comp_bound is None):
                        comp_bound = [copy.deepcopy(pos),pos_i] #this is how we record the location
                    elif(isinstance(pos,ComplexSpecies)):
                        raise ValueError("complex exists in two locations!")
                    pos_i+=1
                if(comp_bound is not None):
                    comp_binders += [comp_bound] #record what is bound, and at what position
            allcomb = self.located_allcomb(comp_binders) #all possible combinations of binders are made here
            allcomb += [[]] #unbound dna should also be used
            for combo in allcomb:
                backbone = possible_backbones[bb_name]
                new_backbone_species = copy.deepcopy(backbone.species)
                new_material = backbone.material_type
                for spec in combo:
                    new_material = "ordered_complex"
                    new_backbone_species[spec[1]] = spec[0]
                newname = None
                if(new_backbone_species==self.get_species().species):
                    newname = self.get_species().name
                combinatorial_complexes += [OrderedComplexSpecies(new_backbone_species,\
                    material_type=new_material,name=newname)] #we make a new OrderedComplexSpecies
        return combinatorial_complexes
    def update_components(self,rnas=None,proteins=None):
        multivalent_self = self.get_species()
        species = [multivalent_self]
        if((rnas is None) or (proteins is None)):
            rnas,proteins = self.explore_txtl()
        #problem is they reach beyond just this DNA_construct so we can't just run
        #explore_integrases here
        #above will be solved by Integrase_construct class in integrases.py
        
        active_attsites = []
        for a in self.parts_list:
            #TODO can we update all components in this loop? Why are promoters and RBSes special?
            #TODO make AttachmentSite a ProteinBindingSite or something
            #pull out all integrase sites.
            if(isinstance(a,AttachmentSite)):
                active_attsites+=[a]
        attsites = self.update_attachment_sites(active_attsites)
        promoters = self.update_promoters(rnas,proteins)
        #rbses = self.update_rbses(proteins)

        active_components = promoters+attsites
        combinatorial_complexes = self.update_combinatorial_complexes(active_components)
        #print("combinatorial complexes are "+str(combinatorial_complexes))
        attsites = self.update_attachment_sites(active_attsites,combinatorial_complexes)
        promoters= self.update_promoters(rnas,proteins,combinatorial_complexes)
        #rbses= self.update_rbses(proteins,combinatorial_complexes)
        out_components = promoters+attsites
        return out_components
    def __hash__(self):
        return hash(self.__repr__())
    def __eq__(self,construct2):
        """equality means comparing the parts list in a way that is not too deep"""
        if(self.__repr__()==construct2.__repr__()):
            return True
        else:
            return False
    
    def update_species(self):
        #TODO make a seperate function for RNA_construct
        species = [self.get_species()]
        rnas = None
        proteins = None
        rnas,proteins = self.explore_txtl()
        #rnas:
        #this is a dictionary of the form:
        #{promoter:rna_construct,promoter2:rna_construct2}
        #proteins:
        #this is a dictionary of the form:
        #{rna_construct:{RBS:[Product1,Product2],RBS2:[Product3]}}
        out_components = self.update_components(rnas,proteins)
        for part in out_components:
            sp_list =  self.remove_bindloc(part.update_species())
            species+=sp_list
        for rna in proteins:
            if(not rna == self):
                species += rna.update_species()
        return species
    def update_reactions(self,norna=False):
        reactions = []
        rnas = None
        proteins = None
        rnas,proteins = self.explore_txtl()
        out_components = self.update_components(rnas,proteins)
        for part in out_components:
            rx_list = []
            for rxn in part.update_reactions():
                new_rxn = copy.copy(rxn)
                new_rxn.inputs = self.remove_bindloc(rxn.inputs)
                new_rxn.outputs = self.remove_bindloc(rxn.outputs)
                rx_list+=[new_rxn]
            reactions+= rx_list
        for rna in proteins:
            if(not rna == self):
                reactions += rna.update_reactions()
        return reactions

class DNAassembly_inprog(DNA_construct):
    def __init__(self, 
                name: str, 
                promoter = None, 
                transcript = None,
                rbs = None, 
                protein = None, 
                length = None,
                attributes = [], 
                mechanisms = {}, 
                parameters = {}, 
                initial_conc = None,
                parameter_warnings = True, 
                **keywords):
        if(isinstance(promoter,str)):
            #if the promoter is a string that means make
            #a default constitutitve promoter
            part_promoter = Promoter(promoter)
        elif(promoter is None):
            part_promoter = Promoter(name)
        if(isinstance(transcript,str)):
            #the name of the transcript must be settable but currently it is not
            self.transcriptName = transcript
        elif(isinstance(transcript,RNA)):
            self.transcriptName = transcript.name
        elif(transcript is None):
            self.transcriptName = name
        
        if(isinstance(rbs,str)):
            part_rbs = RBS(rbs)
        elif(rbs is None):
            part_rbs = RBS(name)
        
        if(isinstance(protein,Protein)):
            part_cds = CDS(protein.name,protein = protein)
        elif(isinstance(protein,str)):
            part_cds = CDS(protein,protein = Protein(protein))
        
        parts_list = [[part_promoter,"forward"],[part_rbs,"forward"],[part_cds,"forward"]]
        DNA_construct(self=self,
                    parts_list=parts_list,
                    name=name,
                    circular=False,
                    mechanisms=mechanisms,  # custom mechanisms
                    parameters=parameters,  # customized parameters
                    attributes=attributes,
                    initial_conc=initial_conc,
                    parameter_warnings = parameter_warnings,
                    **keywords
                    )


class TxTl_Explorer():
    def __init__(self,possible_rxns=("transcription","translation"),direction="forward",parameter_warnings=True):
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
        self.parameter_warnings =parameter_warnings
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

        rna_construct = RNA_construct(copy.deepcopy(rna_partslist),made_by = promoter,parameter_warnings=self.parameter_warnings)
        
        #current_rna_name = str(promoter)+"-"+str(rna_construct)
        self.made_proteins[rna_construct]={}
        #compile the name, keeping track of which promoter made the rna and all the parts on it
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
            self.made_proteins[rna_construct].update({correct_rbs:proteins_per_rbs})

            # so this has to be done at the end of the RNA only because only then do we know
            # exactly what that full RNA is going to be.
            # we want to clear the current proteins, but there could be multiple RNAs that run over the same proteins,
            # so wait until we are done tallying all the RNAs before removing everything from current_proteins
        del self.current_rnas[promoter] # this removes the current RNA from the list, because it's
                                        # being terminated!!
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
    def __init__(self,parts_list,name=None,made_by="unknown",**keywords):
        """an RNA_construct is a lot like a DNA_construct except it can only translate, and
        can only be linear"""
        #TODO make sure we are subclassing RNA and not DNA. I am not sure how to do this,
        #since DNA_construct subclasses DNA
        self.my_promoter = made_by
        DNA_construct.__init__(self=self,parts_list=parts_list,circular=False,name=name,material_type="rna",**keywords)
        #if(name == None):
        #    name = self.make_name()
        #self.name = super().name
        #self.material_type = "rna"
    
    def update_rbses(self,in_proteins,combinatorial_species=None):
        my_rbses = [] #output list of RBSes
        #print(proteins)
        if(len(in_proteins)>1):
            warn("expected one RNA, but we detected more! proteins looks like "+str(in_proteins))
        proteins = in_proteins[self]
        rbses = list(proteins.keys())
        
        #proteins looks like this:
        # {RBS:[protein],RBS1:[protein1,protein2]}
        if(combinatorial_species is None):
            combinatorial_species = [self.get_species()] #list of RNAs
        #print("update_rbs got: combinatorial_species is "+str(combinatorial_species))
        for comb_specie in combinatorial_species: #combinatorial_species will contain RNAs
            
            if(isinstance(comb_specie,OrderedComplexSpecies) and comb_specie.species[-1].name == self.name):
                #TODO below is very similar to what happens in update_promoters. Can we consolidate?
                for rbs in rbses:
                    new_rbs = copy.deepcopy(rbs)
                    if(isinstance(comb_specie.species[rbs.pos],ComplexSpecies)):
                        #if something is already bound then forget it
                        continue
                    new_mv_self = copy.deepcopy(comb_specie)
                    new_mv_self.attributes = ["bindloc_"+str(rbs.pos)] #assign binding location
                    made_proteins = proteins[rbs] #this is a list
                    new_rbs.protein = [a.get_species() for a in made_proteins]
                    new_rbs.transcript = new_mv_self
                    my_rbses+=[new_rbs]
        return my_rbses

    def explore_txtl(self):
        """an RNA has no tx, only TL! central dogma exists, right?"""
        # lets try to make this more modular shall we?
        explorer = TxTl_Explorer(possible_rxns = ("translation",),parameter_warnings=self.parameter_warnings)
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
        rnadict = {self.my_promoter:self}
        return rnadict, proteins
    def get_species(self):
        return OrderedComplexSpecies([a.name for a in self.parts_list]+[self.name],material_type="rna",name=self.name)
    def update_components(self,rnas=None,proteins=None):
        if(proteins is None):
            proteins = self.explore_txtl()
        rbses = self.update_rbses(proteins)
        active_components = rbses
        combinatorial_complexes = self.update_combinatorial_complexes(active_components)
        my_rbses = self.update_rbses(proteins,combinatorial_complexes)
        out_components = my_rbses
        return out_components
    def __repr__(self):
        """the name of an RNA should be different from DNA, right?"""
        output = "rna = "+self.name
        return output



#DEPRECATED circular matching code below
#apparently below won't work because you can't hash the parts_list
    #this means the equality won't work for circular constructs. oops!
    '''
        if(len(self.parts_list) != len(construct2.parts_list)):
            return False
        else:
            if(self.circular == construct2.circular):
                if(self.circular == False):
                    for p1,p2 in zip(self.parts_list,construct2.parts_list):
                        if((p1.name != p2.name) or type(p1)!=type(p2) or p1.pos!=p2.pos):
                            #name and type are the crucial elements
                            #the pos should always match up. If it doesn't, then one of the two
                            #DNA_constructs is super wrong somehow...
                            
                            return False
                    #if you look through the whole list and don't find a mismatch then it's a match
                    return True
                elif(self.circular == True):
                    #in this case we need to rotate them until they align
                    #this opens up a can of worms though. How can we refer to positions of things
                    #if the sequence can be rotated???
                    poslist = [a for a in range(len(construct2.parts_list))]
                    for part_id in range(len(construct2.parts_list)):
                        #try every possible rotation of the thing we're comparing against
                        permuted_id_list = poslist[part_id:]+poslist[:part_id]
                        #we're just making a new list with the parts_list rotated
                        permuted_p2list = [construct2.parts_list[a] for a in permuted_id_list]
                        match = True #if you make it through the loop, then the match is true
                        for p1,p2 in zip(self.parts_list,permuted_p2list):
                            #compare each element to each other element
                            if((p1.name != p2.name) or type(p1)!=type(p2)):
                                match = False
                                #as soon as you find something wrong, quit
                                break
                        if(match==False):
                            #try the next permutation
                            continue
                        else:
                            #if we made it all the way here then it matches!
                            return True
                    #if we are here then we've gone through all permutations and hit "continue" every time.
                    #that means there is no orientation in which it matches
                    return False
    #'''