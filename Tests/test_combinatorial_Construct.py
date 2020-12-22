#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


import pytest
#from unittest import TestCase
import copy

from itertools import permutations
from biocrnpyler import *

def combinatorial_DNAconstruct_RNAconstruct():

    P = Promoter("pconst") #constitutive promoter
    U = RBS("rbs")
    C = CDS("gfp")
    T = Terminator("term")

    correct_order_list = [P,U,C,T]

    directions_list = ["forward","reverse"]

    index_list = range(4)

    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein")),"translation":Translation_MM(Species("Ribo"))}



    iter_list = []

    for order in permutations(correct_order_list):

        for i in range(4):
            
            for direction in directions_list:
                new_order = list(order)
                new_order[i] = [order[i],direction]
                iter_list += [new_order]
    try:
        for combination in iter_list:
            x = DNA_construct(combination,mechanisms = mechs,parameters=parameters)
            print(x)
            y = x.enumerate_components()
            z = []
            #promoter at the beginning facing reverse doesn't make a transcript
            if((isinstance(x[0],Promoter) and x[0].direction=="reverse")):
                assert(y == [x[0]])
            #promoter at the end facing forward doesn't make a transcript
            elif( (isinstance(x[-1],Promoter) and x[-1].direction=="forward")):
                assert(y == [x[-1]])
            #everything else should make one transcript
            else:
                assert(isinstance(y[1],RNA_construct))
            for element in y:
                #find the RNA
                if(isinstance(element,RNA_construct)):
                    rna_enumeration = element.enumerate_components()
                    #reverse RBS does not function in RNA
                    if any([p.direction == 'reverse' for p in element if isinstance(p, RBS)]):
                        assert(len(rna_enumeration)==0)
                    elif(U in element):
                        #if there is an RBS and it is not reverse, then it should be active
                        assert(len(rna_enumeration)>=1)
                    else:
                        #if there is no RBS in the RNA then the RNA cannot do anything
                        assert(len(rna_enumeration)==0)
                            
                    if(len(rna_enumeration)==2):
                        #if two things are returned, that means that a protein is made, and
                        #this is only possible if there is an RBS and CDS in the forward
                        #orientation
                        assert(any([p.direction == 'forward' for p in element if isinstance(p, RBS)]))
                        assert(any([p.direction == 'forward' for p in element if isinstance(p, CDS)]))
    except Exception as e:
        error_txt = f"Combinatorial testing failed when parts list was {combination}. \n Unexpected Error: {str(e)}"
        raise Exception(error_txt)