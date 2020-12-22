#two promoters, rbs, cds, etc
from biocrnpyler import *
import pytest
import copy


def test_two_promoters():
    P = Promoter("pconst") #constitutive promoter
    U = RBS("rbs")
    C = CDS("gfp")
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein")),"translation":Translation_MM(Species("Ribo"))}

    #two promoters
    x = DNA_construct([P,P],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(y)==5) #enumerated components should have five things
    assert(x[0] in y) #first promoter is active
    assert(x[1] in y) #second promoter is active
    assert(len(rnas)==1) #one RNAconstruct is made
    assert(x[0].transcript==rnas[0].get_species()) #promoter makes the correct RNA
    assert(x[1].transcript is None) #second promoter makes no transcript

    #circular construct
    x = DNA_construct([P,P],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(y)==6) #enumerated components should have six things now because the second RNA goes around the plasmid
    assert(x[0] in y) #first promoter is active
    assert(x[1] in y) #second promoter is active
    assert(len(rnas)==2) #twp RNAconstructs are made
    assert(x[0].transcript in [a.get_species() for a in rnas]) #promoter makes a transcript
    assert(x[1].transcript in [a.get_species() for a in rnas]) #second promoter makes a transcript

def test_two_RBSes():
    P = Promoter("pconst") #constitutive promoter
    U = RBS("rbs")
    C = CDS("gfp")
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein")),"translation":Translation_MM(Species("Ribo"))}

    #two RBSes
    x = DNA_construct([P,U,U],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(rnas)==1) #one RNA is made
    z = rnas[0].enumerate_components()
    print(z)
    assert(len(z)==4) #four RBSes are made, but no proteins
    assert(rnas[0][0] in z) #RBS from RNA is active
    assert(rnas[0][1] in z) #second RBS is also active

    #two RBSes with a protein
    x = DNA_construct([P,U,U,C],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(len(y)==2) #two components are made
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(rnas)==1) #one RNA is made
    z = rnas[0].enumerate_components()
    assert(sum([isinstance(a,CDS) for a in z])==1) #only one protein is made

    #circularly permuted
    x = DNA_construct([U,U,C,P],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    assert(len(y)==2) #two components are made
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(rnas)==1) #one RNA is made
    z = rnas[0].enumerate_components()
    assert(len(z)==5) #four RBS parts and one CDS
    assert(sum([isinstance(a,CDS) for a in z])==1) #only one protein is made

    #two RBSes with two proteins

    x = DNA_construct([P,U,C,U,C],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(len(y)==2) #two components are made
    rnas = []
    for a in y:
        if(isinstance(a,RNA_construct)):
            rnas += [a]
    assert(len(rnas)==1) #one RNA is made
    z = rnas[0].enumerate_components()
    assert(sum([isinstance(a,CDS) for a in z])==2) #two proteins are made
    assert(len(z)==6) #four RBS parts and two CDS