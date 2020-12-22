#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


import pytest
#from unittest import TestCase
import copy

from itertools import permutations
from biocrnpyler import Promoter, RBS, CDS, Terminator, Transcription_MM, Translation_MM, RNA_construct, DNA_construct, Species

def test_combinatorial_DNAconstruct_RNAconstruct():

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
        error_txt = f"Combinatorial enumerate_components failed when parts list was {combination}. \n Unexpected Error: {str(e)}"
        raise Exception(error_txt)



from biocrnpyler import Promoter, ActivatablePromoter, RepressiblePromoter, RegulatedPromoter, CombinatorialPromoter

#list of possible promoter types
promoters = [
    Promoter("promoter"),
    ActivatablePromoter("activatable_promoter", activator = Species("A")),
    RepressiblePromoter("repressible_promoter", repressor = Species("R")),
    RegulatedPromoter("regulated_promoter", regulators = [Species("S1"), Species("S2")]),
    CombinatorialPromoter("combinatorial_promoter", regulators = [Species("S1"), Species("S2")])
    ]

from biocrnpyler import Terminator
#list of possible terminators
terminators = [
    Terminator("terminator")
    ]

from biocrnpyler import SimpleTxTlExtract, TxTlExtract, SimpleTxTlDilutionMixture, TxTlDilutionMixture
#list of possible mixtures
parameters = {"kb":1.0, "ku":1.0, "ktx":1.0, "ktl":1.0, "kdeg":1.0, "kdil":1.0, "kexpress":1.0,"kcat":1.0, "K":10, 
    "cooperativity":2, "n":2, "k":1.0, "kleak":1.0}

#Mixtures have to be instantiated each time, to reset their internal variables
mixture_classes = [
    (SimpleTxTlExtract, {"name":"simple_tx_tl_extract"}),
    (TxTlExtract, {"name":"tx_tl_extract"}),
    (SimpleTxTlDilutionMixture, {'name':"simple_tx_tl_dilution_mixture"}),
    (TxTlDilutionMixture, {"name":"tx_tl_dilution_mixture"})
    ]

#Helper function to do sanity checks on CRN reaction inputs and outputs
def all_reaction_inputs_and_outputs(CRN):
    inputs = []
    outputs = []
    for r in CRN.reactions:
        inputs += [w.species for w in r.inputs]
        outputs += [w.species for w in r.outputs]

    return inputs, outputs

def test_combinatorial_DNAconstruct_in_Mixtures():
    #Tests the construct [Promoter, Terminator] in many different mixtures
    try:
        for mclass, args in mixture_classes:
            #create the Mixture with parameters
            args["parameters"] = dict(parameters)
            m = mclass(**args)

            #Promoters and terminators are instances
            for p in promoters:
                for t in terminators:
                    construct = DNA_construct([p, t])
                    m.add_component(construct)
                    #Get the construct out of the Mixture
                    construct = m.get_component(construct)
                    crn = m.compile_crn()
                    inputs, outputs = all_reaction_inputs_and_outputs(crn)

                    #the DNA construct should be the input of a reaction
                    assert construct.get_species() in inputs

                    #the DNA_construct's promoter should produce a transcript
                    assert construct[0].transcript in outputs

    except Exception as e:
        error_txt = f"Combinatorial compilation failed when DNA_constrct parts list was {[p, t]} in mixture {m}. \n Unexpected Error: {str(e)}"
        raise Exception(error_txt)


#list of CDSs
CDSs = [CDS("cds")]

#list of RBSs
RBSs = [RBS("rbs")]

def test_combinatorial_RNAconstruct_in_Mixtures():
    #Tests the RNA_construct [RBS, CDS] in many different mixtures
    try:
        for mclass, args in mixture_classes:
            #create the Mixture with parameters
            args["parameters"] = dict(parameters)
            m = mclass(**args)

            #RBSs and CDSs are instances
            for rbs in RBSs:
                for cds in CDSs:
                    construct = RNA_construct([rbs, cds])
                    m.add_component(construct)
                    #Get the construct out of the Mixture
                    construct = m.get_component(construct)
                    crn = m.compile_crn()
                    inputs, outputs = all_reaction_inputs_and_outputs(crn)

                    #the RNA construct should be the input of a reaction
                    assert construct.get_species() in inputs

                    #the RNA_construct's CDS should produce a transcript
                    assert construct[1].protein in outputs

    except Exception as e:
        error_txt = f"Combinatorial compilation failed when RNA_construct parts list was {[rbs, cds]} in mixture {m}. \n Unexpected Error: {str(e)}"
        raise Exception(error_txt)

def test_combinatorial_DNAconstruct_RNAconstruct_in_Mixtures():
    #Tests the DNA_construct [Promoter, RBS, CDS, Terminator] in many different mixtures
    #This will also produce an RNA_construct during component enumeration and can be seen
    #as a combination of the above two tests
    try:
        for mclass, args in mixture_classes:
            #create the Mixture with parameters
            args["parameters"] = dict(parameters)
            m = mclass(**args)

            #Promoters and terminators are instances
            for p in promoters:
                for t in terminators:
                    for rbs in RBSs:
                        for cds in CDSs:
                            construct = DNA_construct([p, rbs, cds, t])

                            m.add_component(construct)
                            #Get the construct out of the Mixture
                            construct = m.get_component(construct)
                            crn = m.compile_crn()
                            inputs, outputs = all_reaction_inputs_and_outputs(crn)

                            #the DNA construct should be the input of a reaction
                            assert construct.get_species() in inputs

                            #the DNA_construct's promoter should produce a transcript
                            assert construct[0].transcript in outputs

                            #the RNA construct transcript should be the input of a reaction
                            assert construct[0].transcript in inputs

                            #the RNA_construct's CDS should produce a transcript
                            assert construct[2].protein in outputs

    except Exception as e:
        error_txt = f"Combinatorial compilation failed when DNA_construct parts list was {[p, rbs, cds, t]} in mixture {m}. \n Unexpected Error: {str(e)}"
        raise Exception(error_txt)