
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler import Promoter, DNA_construct, Terminator, Transcription_MM,\
                         Species, RBS, CDS, Polymer_transformation, NamedPolymer,\
                         OrderedPolymer, IntegraseSite, IntegraseRule
import copy

def test_polymer_transformation():
    prom = Promoter("p1")
    utr = RBS("utr")
    cds = CDS("GFP")
    cds2 = CDS("RFP")
    term = Terminator("t16")

    construct = DNA_construct([prom,utr,cds,term])
    construct2 = DNA_construct([cds2])
    dummy = Polymer_transformation.dummify(construct,"input1")
    assert(dummy.name == "input1") #dummy polymer has the right name
    assert(len(dummy.polymer)==len(construct)) #dummy polymer has the right length
    assert([a.direction for a in dummy]==[a.direction for a in construct]) #dummy polymer monomers have the right direction
    assert([a.position for a in dummy]==[a.position for a in construct]) #dummy polymer monomers have the right position
    assert(isinstance(dummy,NamedPolymer)) #make sure the output of dummify has the right type
    assert(isinstance(dummy,OrderedPolymer)) #make sure the output of dummify has the right type
    dummy2 = Polymer_transformation.dummify(construct,"input2")
    transformation = Polymer_transformation([dummy[0],[prom,"forward"],dummy2[0]],circular = True)
    assert(transformation.circular == True) #make sure circular is set properly
    assert(len(transformation.parentsdict)==2) #two parents!
    assert(transformation.create_polymer([construct,construct2])==DNA_construct([prom,prom,cds2],circular=True)) #transformations are performed correctly
    assert(transformation.create_polymer([construct2,construct])==DNA_construct([cds2,prom,prom],circular=True))
    assert("(input1-0)(p1)(input2-0)" in str(transformation)) #contains vital information in the repr

def test_integrase_rule():
    prom = Promoter("p1")
    utr = RBS("utr")
    cds = CDS("GFP")
    cds2 = CDS("RFP")
    term = Terminator("t16")
    ap = IntegraseSite("attP","attP")
    ab = IntegraseSite("attB","attB")
    aflp = IntegraseSite("FLP","FLP")
    delete = DNA_construct([ab,cds,ap])
    flip = DNA_construct([ab,cds,[ap,"reverse"]])


    bxb1_rule = IntegraseRule(name="Bxb1",reactions={("attB","attP"):"attL",("attP","attB"):"attR"})
    assert(set(["attP","attB","attL","attR"])==set(bxb1_rule.attsites))
    assert(set(["attP","attB","attL","attR"])==set(bxb1_rule.binds_to()))
    assert(bxb1_rule.integrase_species == Species("Bxb1",material_type="protein"))

    
    flp_rule = IntegraseRule(name="FLP",reactions={("FLP","FLP"):"FLP"})
    bad_rule = IntegraseRule(name="Bxb1",reactions={("attB","attP"):"attL"})
    with pytest.raises(AssertionError):
        bxb1_rule.generate_products(delete[0],delete[2]) #must have the right integrase
    delete[0].update_integrase("Bxb1")
    delete[2].update_integrase("Bxb1")
    with pytest.raises(KeyError):
        bxb1_rule.generate_products(delete[0],delete[0]) #attP+attP is not possible
    with pytest.raises(KeyError):
        bad_rule.generate_products(delete[0],delete[2]) #attP->attB reaction is not possible with a bad rule
    
    aL = IntegraseSite("attL","attL",integrase="Bxb1")
    aR = IntegraseSite("attR","attR",integrase="Bxb1")
    productsites = bxb1_rule.generate_products(delete[0],delete[2])
    assert(productsites[0]==aL.set_dir("forward"))
    assert(productsites[1]==aR.set_dir("forward"))
    #
        
