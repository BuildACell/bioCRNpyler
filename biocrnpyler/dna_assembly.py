#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from .component import Component, DNA
from .chemical_reaction_network import Species
from .mechanism import One_Step_Cooperative_Binding
from warnings import warn as pywarn
import itertools as it
import numpy as np

def warn(txt):
    pywarn(txt)


class Promoter(Component):
    def __init__(self, name, assembly=None,
                 transcript=None, length=0,
                 mechanisms={}, parameters={}, **keywords):
        self.assembly = assembly
        self.length = length
        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type="rna")
        else:
            self.transcript = self.set_species(transcript, material_type = 'rna')

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        species = []
        species += mech_tx.update_species(dna = self.assembly.dna, \
            transcript = self.transcript, protein = self.assembly.protein)
        return species

    def update_reactions(self):
        mech_tx = self.mechanisms["transcription"]
        reactions = []


        reactions += mech_tx.update_reactions(dna = self.assembly.dna, \
                        component = self, part_id = self.name, complex = None,
                        transcript = self.transcript, protein = self.assembly.protein)
        return reactions

class RegulatedPromoter(Promoter):
    def __init__(self, name, regulators, leak = True, assembly = None,
                 transcript = None, length = 0, mechanisms = {},
                 parameters = {}, **keywords):

        if not isinstance(regulators, list):
            regulators = [regulators]

        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]

        self.leak = leak

        self.default_mechanisms = {"binding": One_Step_Cooperative_Binding()}

        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters,
                          **keywords)

    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']
        species = []
        self.complexes = []
        if self.leak is not False:
            species += mech_tx.update_species(dna = self.assembly.dna)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            species_b = mech_b.update_species(regulator, self.assembly.dna,\
                                            part_id = self.name,component=self)
            species += species_b
            complex_ = species_b[0]
            self.complexes += [complex_]
            species += mech_tx.update_species(dna = complex_, \
                    transcript = self.transcript, protein = self.assembly.protein)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak != False:
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, component = self, part_id = self.name, \
                                                            transcript = self.transcript, protein = self.assembly.protein)

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]

            reactions += mech_b.update_reactions(regulator, self.assembly.dna, component = self, \
                                                                    part_id = self.name+"_"+regulator.name)
            reactions += mech_tx.update_reactions(dna = complex_, component = self, part_id = self.name+"_"+regulator.name, \
                                            transcript = self.transcript, protein = self.assembly.protein)

        return reactions

class CombinatorialPromoter(Promoter):
    def __init__(self, name, regulators, leak = True, assembly = None,
                 transcript = None, length = 0, mechanisms = {},
                 parameters = {},tx_capable_list = None, **keywords):
        """
        tx_capable_list = [[1,2,3],[0,2,3]] this means having regulator 1, 2, and 3 will transcribe
                                           but also 3, 2, and 0.
                                           #TODO make it force sorted list
        """
        if not isinstance(regulators, list):
            regulators = [regulators]
        if(tx_capable_list == None):
            tx_capable_list = [[a for a in range(len(regulators))]]
        self.tx_capable_list = tx_capable_list #TODO make this better
        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]
        self.regulators = sorted(self.regulators)
        self.leak = leak

        self.default_mechanisms = {"binding": One_Step_Cooperative_Binding()}

        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters,
                          **keywords)
        self.complex_combinations = {}
    def dp_complex_combination(self,key,dp_dict,core_dna=None,\
                                            mech_b=None):
        """retrieves a key from a list, otherwise, makes an element
        and stores it"""

        if(core_dna==None):
            core_dna = self.assembly.dna
        if(mech_b == None):
            mech_b = self.mechanisms['binding']
        if(key in dp_dict):
            return dp_dict[key]
        else:
            made_complex = core_dna
            for element in key:
                made_complex = mech_b.update_species(element,made_complex,\
                                        part_id = self.name,component=self)[0]
            dp_dict[key]=made_complex
            return made_complex
    def update_species(self):
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']
        
        species = []

        
        self.complexes = []
        if self.leak is not False:
            species += mech_tx.update_species(dna = self.assembly.dna)
        complex_combinations = {}
        for i in range(1, len(self.regulators)+1):
            # Get all unique complexes of len i
            #species += [self.regulators[i]]
            for combination in it.combinations(self.regulators, i):

                temp_complex = self.dp_complex_combination(combination,self.complex_combinations)
                species+=[temp_complex]
                if([self.regulators.index(a) for a in sorted(combination)] in self.tx_capable_list):
                    species += mech_tx.update_species(dna = temp_complex, transcript = self.transcript, protein = self.assembly.protein)
        print(species)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.mechanisms["transcription"]
        mech_b = self.mechanisms['binding']

        if self.leak != False:
            reactions += mech_tx.update_reactions(dna = self.assembly.dna, component = self, part_id = self.name, \
                                                            transcript = self.transcript, protein = self.assembly.protein)

        for i in range(1, len(self.regulators)+1):
            # Get all unique complexes of len i
            for bigger in it.combinations(self.regulators, i):
                # Get all unique complexes len i-1
                for smaller in it.combinations(self.regulators, i-1):
                    # If smaller is only 1 activator away
                    if (bigger > smaller):
                        #this is the regulator that is different
                        regulator_to_add = np.setdiff1d(bigger,smaller)[0]
                        #so now we have to come up with the "bigger" complex
                        big_complex = self.dp_complex_combination(bigger,self.complex_combinations)
                        #now we do the same thing for "smaller"
                        small_complex = self.dp_complex_combination(smaller,self.complex_combinations)
                        reactions+=mech_b.update_reactions(regulator_to_add,small_complex,\
                                  complex=big_complex,part_id = self.name, component = self,)
                        if([self.regulators.index(a) for a in sorted(bigger)] in self.tx_capable_list):
                            mech_tx.update_reactions(dna=big_complex,component=self,\
                                            transcript=self.transcript,protein=self.assembly.protein)
        
        return reactions


class RBS(Component):
    def __init__(self, name, assembly = None,
                 transcript = None, protein = None, length = 0,
                 mechanisms = {}, parameters = {}, **keywords):
        self.assembly = assembly
        self.length = length

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type = "rna")
        else:
            self.transcript = self.set_species(transcript, material_type = "rna")
 
        if protein is None:
            self.protein = Species(assembly.name, material_type = "protein")
        else:
            self.protein = self.set_species(protein, material_type = "protein")

        Component.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)

    def update_species(self):
        mech_tl = self.mechanisms['translation']
        species = []
        species += mech_tl.update_species(transcript = self.transcript, protein = self.protein)
        return species

    def update_reactions(self):
        mech_tl = self.mechanisms['translation']
        reactions = []

        reactions += mech_tl.update_reactions(transcript = self.transcript,
                                              protein = self.protein, component = self, part_id = self.name)
        return reactions


class DNAassembly(DNA):
    def __init__(self, name, dna = None, promoter = None, transcript = None,
                 rbs = None, protein = None, length = None,
                 attributes = [], mechanisms = {}, parameters = {}, initial_conc = None,
                 parameter_warnings = True, **keywords):
        self.promoter = None
        self.rbs = None
        self.transcript = None
        self.initial_concentration = initial_conc
        self.name = name

        DNA.__init__(self, name, length = length, mechanisms = mechanisms,
                     parameters = parameters, initial_conc = initial_conc,
                     parameter_warnings = parameter_warnings,
                     attributes = list(attributes), **keywords)

        self.update_dna(dna, attributes = list(attributes))
        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_promoter(promoter, transcript = self.transcript)
        self.update_rbs(rbs, transcript = self.transcript,
                        protein = self.protein)

        self.set_parameter_warnings(parameter_warnings)

            

    def set_parameter_warnings(self, parameter_warnings):
        self.parameter_warnings = parameter_warnings

        if self.parameter_warnings != None:
            if self.promoter != None:
                self.promoter.set_parameter_warnings(parameter_warnings)
            if self.rbs != None:
                self.rbs.set_parameter_warnings(parameter_warnings)


    def update_dna(self, dna, attributes = None):
        if dna is None:
            self.dna = self.set_species(self.name, material_type = "dna", attributes = attributes)
        else:
            self.dna = self.set_species(dna, material_type = "dna", attributes = attributes)
        

    def update_transcript(self, transcript, attributes = None):
        if transcript is None:
            self.transcript = self.set_species(self.name, material_type = "rna", attributes = attributes)
        else:
            self.transcript = self.set_species(transcript, material_type = "rna", attributes = attributes)

        if self.promoter is not None:
            self.promoter.transcript = self.transcript
        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(self, protein, attributes = None):
        if protein is None:
            self._protein = self.set_species(self.name, material_type = "protein", attributes = attributes)
        else:
            self._protein = self.set_species(protein, material_type = "protein", attributes = attributes)

        if self.rbs is not None:
            self.rbs.transcript = self.protein

    def update_promoter(self, promoter, transcript=None):
        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(promoter, str):
            self.promoter = Promoter(assembly = self, name = promoter,
                                     transcript = self.transcript,
                                     parameters = self.parameters)
        elif isinstance(promoter, Promoter):
            self.promoter = promoter
            self.promoter.assembly = self
            self.promoter.transcript = self.transcript
        elif promoter is not None:
            raise ValueError("Improper promoter type recieved by DNAassembly. "
                             "Expected string or promoter object. "
                             f"Recieved {repr(promoter)}.")
        if promoter is not None:
            self.promoter.update_parameters(
                                        mixture_parameters = self.parameters,
                                        overwrite_custom_parameters = False)

    def update_rbs(self, rbs, transcript = None, protein = None):
        if protein is not None:
            self.update_protein(protein)

        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(rbs, str):
            self.rbs = RBS(assembly = self, name = rbs, protein = self._protein,
                           transcript = self.transcript,
                           parameters = self.parameters)
        elif isinstance(rbs, RBS):
            self.rbs = rbs
            self.rbs.assembly = self
            self.rbs.transcript = self.transcript
            self.rbs.protein = self._protein
        elif rbs is not None:
            raise ValueError("Improper rbs type recieved by DNAassembly. "
                             "Expected string or RBS object. Recieved "
                            f"{repr(rbs)}.")

        if rbs is not None:
            self.rbs.update_parameters(mixture_parameters = self.parameters,
                                       overwrite_custom_parameters = False)

    @property
    def protein(self):
        return self._protein

    def update_species(self):
        species = []
        species.append(self.dna)
        if self.promoter is not None and self.rbs is not None:
            species += self.promoter.update_species()
            species += self.rbs.update_species()

        elif self.promoter is not None and self.rbs is None:
            species += self.promoter.update_species()

        if "rna_degredation" in self.mechanisms and self.promoter is not None:
            deg_mech = self.mechanisms["rna_degredation"]
            species += deg_mech.update_species(rna = self.transcript)

        # TODO raise a warning if there were duplicate species
        return list(set(species))

    def update_reactions(self):
        reactions = []
        if self.promoter is not None:
            self.promoter.parameter_warnings = self.parameter_warnings
            reactions += self.promoter.update_reactions()

        if self.rbs is not None:
            self.rbs.parameter_warnings = self.parameter_warnings
            reactions += self.rbs.update_reactions()

        if "rna_degredation" in self.mechanisms and self.promoter is not None:
            deg_mech = self.mechanisms["rna_degredation"]


            reactions += deg_mech.update_reactions(rna = self.transcript, component = self.promoter)
        # TODO check that the reaction list is unique
        return reactions

    def update_parameters(self, mixture_parameters = {}, parameters = {},
                          overwrite_custom_parameters = True):
        DNA.update_parameters(self = self,
                              mixture_parameters = mixture_parameters,
                              parameters = parameters)

        if self.promoter is not None:
            self.promoter.update_parameters(
                                    mixture_parameters = mixture_parameters,
                                    parameters = parameters,
                                    overwrite_custom_parameters = False)
        if self.rbs is not None:
            self.rbs.update_parameters(mixture_parameters = mixture_parameters,
                                       parameters = parameters,
                                       overwrite_custom_parameters = False)

    def update_mechanisms(self, mixture_mechanisms = {}, mechanisms = {},
                          overwrite_custom_parameters = True):
        DNA.update_mechanisms(self = self,
                              mixture_mechanisms = mixture_mechanisms,
                              mechanisms = mechanisms)

        if self.promoter is not None and "transcription" in self.mechanisms:
            mech_tx = self.mechanisms["transcription"]
            mechs = {"transcription": mech_tx}
            self.promoter.update_mechanisms(mechanisms = mechs,
                                            overwrite_custom_mechanisms = False)
        if self.rbs is not None and "translation" in self.mechanisms:
            mech_tl = self.mechanisms["translation"]
            mechs = {"translation": mech_tl}
            self.rbs.update_mechanisms(mechanisms = mechs,
                                       overwrite_custom_mechanisms = False)

    def __str__(self):
        return type(self).__name__ + ": " + self.name

    def __repr__(self):
        txt = str(self)
        if self.promoter is not None:
            txt += "\n\t" + repr(self.promoter)
            txt += "\n\ttranscript = " + repr(self.transcript)
        if self.rbs is not None:
            txt += "\n\t" + repr(self.rbs)
            txt += "\n\tprotein = " + repr(self._protein)
        return txt





#ANDREY'S STUFF BELOW!!! here be dagrons
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
