
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler.components import Promoter, DNA_construct, Terminator,\
    Species, RBS, CDS, Complex
from biocrnpyler.mechanisms import Transcription_MM
import copy

def test_promoter_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    x = DNA_construct([P],mechanisms = mechs,parameters=parameters)
    assert(not(P in x.parts_list)) #make sure P is copied
    assert(P in x) #testing "contains" function of Construct
    y = x.enumerate_components()
    assert([x[0]] == y)
    assert(x[0].transcript is None)

    #works in reverse as well
    x = DNA_construct([[P,"reverse"]],mechanisms = mechs,parameters=parameters)
    assert(not(P in x.parts_list)) #make sure P is copied
    assert(P in x) #testing "contains" function of Construct
    y = x.enumerate_components()
    assert([x[0]] == y)
    assert(x[0].transcript is None)

def test_promoter_terminator_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    #minimal RNA transcription
    x = DNA_construct([P,T],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

    #"incoherent" transcription
    x = DNA_construct([[P,"reverse"],T],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y == [x[0]]) #correct promoter is returned

    #transcription in reverse works as well
    x = DNA_construct([T,[P,"reverse"]],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y[0] == x[1]) #correct promoter is returned
    assert(x[1].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[1]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

def test_get_part_DNAconstruct():
    P = Promoter("pconst")
    U = RBS("utr1")
    C = CDS("coolprot1")
    T = Terminator("T11")
    T2 = Terminator("iamsoterminator")

    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    x = DNA_construct([P,U,C,T,T2],mechanisms = mechs,parameters=parameters)
    with pytest.raises(ValueError, match=r"get_component requires a single keyword. Recieved component=None, name=None, index=None."):
        y = x.get_part()
    with pytest.raises(ValueError):
        y = x.get_part(part = P,index=3)
    with pytest.raises(ValueError):
        y = x.get_part(part = "oogabooga")
    with pytest.raises(ValueError):
        y = x.get_part(part_type = "oogabooga")
    with pytest.raises(ValueError):
        y = x.get_part(name = 5)
    with pytest.raises(ValueError):
        y = x.get_part(index = "oogabooga")
    assert(x.get_part(part_type=Promoter)==x[0])
    assert(x.get_part(name="pconst")==x[0])
    assert(x.get_part(index=2)==x[2])
    assert(x.get_part(part=P)==x[0])
    assert(len(x.get_part(part_type=Terminator))==2)

def test_circular_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    #circular construct
    x = DNA_construct([P,T],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs


    x = DNA_construct([[P,"reverse"],T],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()
    assert(y[0] == x[0]) #correct promoter is returned
    assert(x[0].transcript == y[1].get_species()) #promoter in the DNA_construct has the correct transcript species
    assert(y[1].promoter == x[0]) #rna is linked back to the DNA_construct's promoter
    assert(y[1].promoter == y[0]) #correct promoter is returned by enumerate_constructs

    x = DNA_construct([[P,"reverse"],[T,"reverse"],],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.enumerate_components()

    z = DNA_construct([T,P],mechanisms = mechs,parameters=parameters,circular=True)
    e = z.enumerate_components()
    assert(y[1]==e[1]) #different DNA constructs lead to the same RNA construct
    assert(y[0]==x[0]) #correct promoter is working
    assert(e[0]==z[1]) #correct promoter is working

def test_combinatorial_enumeration_DNAconstruct():
    P = Promoter("pconst") #constitutive promoter
    T = Terminator("term")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"transcription":Transcription_MM(Species("RNAP",material_type="protein"))}

    #circular construct
    x = DNA_construct([P,P,T],mechanisms = mechs,parameters=parameters,circular=True)
    y = x.combinatorial_enumeration()
    prom_positions = []
    dna_species = {}
    for prom_comp in y:
        if(prom_comp.position not in dna_species):
            dna_species[prom_comp.position] = [prom_comp.dna_to_bind.parent]
        else:
            dna_species[prom_comp.position] += [prom_comp.dna_to_bind.parent]
    
    
    p1_bound = Complex([x.get_species()[0],Species("RNAP",material_type="protein")]).parent
    p2_bound = Complex([x.get_species()[1],Species("RNAP",material_type="protein")]).parent
    both_bound = Complex([p1_bound[1],Species("RNAP",material_type="protein")]).parent
    assert(x.get_species() in dna_species[0]) #unbound polymer found in combinatorial enumeration
    assert(p2_bound in dna_species[0]) #proper bound polymer found in combinatorial enumeration
    assert(x.get_species() in dna_species[1]) #unbound polymer found in combinatorial enumeration
    assert(p1_bound in dna_species[1]) #proper bound polymer found in combinatorial enumeration
    assert(both_bound not in dna_species[0]) #dna with both promoters bound is not in combinatorial enumeration
    assert(both_bound not in dna_species[1]) #dna with both promoters bound is not in combinatorial enumeration
    assert(len(y)==4)