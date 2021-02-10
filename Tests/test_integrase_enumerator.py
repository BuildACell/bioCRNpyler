
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler import Promoter, DNA_construct, Terminator, Transcription_MM,\
                         Species, RBS, CDS, Polymer_transformation, NamedPolymer,\
                         OrderedPolymer, AttachmentSite, IntegraseRule, Integrase_Enumerator, TxTlExtract
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
    assert(len(dummy._polymer)==len(construct)) #dummy polymer has the right length
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
    ap = AttachmentSite("attP","attP")
    ab = AttachmentSite("attB","attB")
    aflp = AttachmentSite("FLP","FLP")
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
    
    aL = AttachmentSite("attL","attL",integrase="Bxb1")
    aR = AttachmentSite("attR","attR",integrase="Bxb1")
    productsites = bxb1_rule.generate_products(delete[0],delete[2])
    assert(productsites[0]==aL.set_dir("forward"))
    assert(productsites[1]==aR.set_dir("forward"))
    #
        
def test_compilation():
    #Create an infinite polymer system and compile it at different recursion depths

    pconst = Promoter("pconst") #constitutive promoter
    gfp = CDS("GFP")
    rfp = CDS("RFP")
    t16 = Terminator("t16") #a terminator stops transcription

    #some parameters are useful also. note the "kint" parameter which determines the rate of recombination
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    #here is where we define the integrase attachment sites
    attP = AttachmentSite("attP","attP",integrase="Bxb1") #the first argument is the name of this attachment site, the second argument is the type of attachment site it is ("attP", "attB", "attL" or "attR") and the integrase denotes which integrase binds to this site (default is "int1")
    attB = AttachmentSite("attB","attB",integrase="Bxb1") #we define two attachment sites, as one site doesn't do anything on its own besides bind integrases

    #Create an integrase enumerator
    bxb1_mechanism = IntegraseRule("Bxb1", reactions={("attB","attP"):"attL",("attP","attB"):"attR"})
    bxb1 = Integrase_Enumerator("Bxb1", int_mechanisms={"Bxb1":bxb1_mechanism}) #we must also define an integrase enumerator. The default integrase is always "int1", but here we are specifying "Bxb1" as the integrase. The Integrase enumerator gets a name, and the mechanism we defined above.

    #create DNA_constructs with integrase sites
    plasmid1_construct = DNA_construct([t16,attP,attB,gfp],circular=True)
    genome_construct = DNA_construct([pconst,attB,rfp])

    #Create a Mixture
    myMixture = TxTlExtract(name = "txtl", parameters = parameters, components = [plasmid1_construct,genome_construct],global_component_enumerators=[bxb1]) 

    #Test Recursion Depth = 0
    #Note compile directives are used to speed up the test
    myCRN0, comps0 = myMixture.compile_crn(recursion_depth = 0, 
                                         return_enumerated_components = True,
                                         initial_concentrations_at_end = True,
                                         copy_objects = False,
                                         add_reaction_species = False)
    assert len(comps0) == 14
    assert len(myCRN0.species) == 14
    assert len(myCRN0.reactions) == 12

    #Test recursion depth = 1
    myCRN1, comps1 = myMixture.compile_crn(recursion_depth = 1, 
                                         return_enumerated_components = True,
                                         initial_concentrations_at_end = True,
                                         copy_objects = False,
                                         add_reaction_species = False)

    assert len(comps1) == 52
    assert len(myCRN1.species) == 36
    assert len(myCRN1.reactions) == 61

    #Recompiling at a different length doesn't change anything
    myCRN0b, comps0b = myMixture.compile_crn(recursion_depth = 0, 
                                         return_enumerated_components = True,
                                         initial_concentrations_at_end = True,
                                         copy_objects = False,
                                         add_reaction_species = False)
    assert len(comps0) == len(comps0b)
    assert all(myCRN0.species == myCRN0b.species)
    assert all(myCRN0.reactions == myCRN0b.reactions)