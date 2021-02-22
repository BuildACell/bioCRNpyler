
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler import Promoter, DNA_construct, Terminator, Transcription_MM,\
                         Species, RBS, CDS, Polymer_transformation, NamedPolymer,\
                         OrderedPolymer, IntegraseSite, IntegraseRule,Integrase_Enumerator, TxTlExtract
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
    plasp = DNA_construct([cds,ap],circular=True)
    plasb = DNA_construct([cds,ab],circular=True)
    genp = DNA_construct([cds2,ap])
    genb = DNA_construct([cds2,ab])


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
    

    ap = IntegraseSite("attP","attP",integrase="Bxb1")
    ab = IntegraseSite("attB","attB",integrase="Bxb1")
    delete = DNA_construct([ab,cds,ap])
    flip = DNA_construct([ab,cds,[ap,"reverse"]],circular=True)
    plasp = DNA_construct([cds,ap],circular=True)
    plasb = DNA_construct([cds,ab],circular=True)
    genp = DNA_construct([cds2,ap])
    genb = DNA_construct([cds2,ab])
    #testing the actual integration reactions
    #inversion
    integ_funcs = bxb1_rule.integrate(flip[0],flip[2]) #sites forwards
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[None,1,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["forward","reverse","reverse"])
    integ_funcs = bxb1_rule.integrate(flip[2],flip[0]) #sites reverse
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[None,1,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["reverse","reverse","forward"])

    integ_funcs = bxb1_rule.integrate(flip[0],flip[2],also_inter=True) #allow intermolecular reactions
    assert(len(integ_funcs)==2)#now, two products are made
    assert(sum([a.circular for a in integ_funcs])==2) #all of them is circular
    for prod in integ_funcs:
        if(len(prod.partslist)==3):
            assert([a[0].position for a in prod.partslist]==[None,1,None])
            assert([a[1] for a in prod.partslist]==["forward","reverse","reverse"])
        else:
            assert([a[0].position for a in prod.partslist]==[None,1,0,None,1,2])
            assert([a[1] for a in prod.partslist]==["forward","reverse","reverse","forward","forward","reverse"])
    #deletion
    integ_funcs = bxb1_rule.integrate(delete[0],delete[2]) #sites forwards
    assert(len(integ_funcs)==2)#two product is made
    assert(sum([a.circular for a in integ_funcs])==1) #one of them is circular
    for prod in integ_funcs:
        if(prod.circular):
            assert([a[0].position for a in prod.partslist]==[None,1])
            assert([a[1] for a in prod.partslist]==["forward","forward"])
        else:
            assert([a[0].position for a in prod.partslist]==[None])
            assert([a[1] for a in prod.partslist]==["forward"])
    #now, what if we also calculate the intermolecular reactions?
    integ_funcs = bxb1_rule.integrate(delete[2],delete[0],also_inter=True) #sites reverse
    print(integ_funcs)
    assert(len(integ_funcs)==4)#four product is made
    assert(sum([a.circular for a in integ_funcs])==2) #two of them is circular
    for prod in integ_funcs:
        if(prod.circular):
            if(len(prod.partslist)==2):
                assert([a[0].position for a in prod.partslist]==[None,1])
                assert([a[1] for a in prod.partslist]==["forward","forward"])
                assert(prod.number_of_inputs == 1)
            else:
                assert([a[0].position for a in prod.partslist]==[None,1,None,1])
                assert([a[1] for a in prod.partslist]==["forward","forward","reverse","reverse"])
                assert(prod.number_of_inputs == 2)
        else:
            assert([a[0].position for a in prod.partslist]==[None])
            assert([a[1] for a in prod.partslist]==["forward"])
            assert(prod.number_of_inputs == 1)


    #integration
    integ_funcs = bxb1_rule.integrate(plasp[1],plasb[1]) #two plasmids -> one plasmid p then b
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[0,None,0,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["forward","forward","forward","forward"])
    assert(integ_funcs[0].circular)
    assert(integ_funcs[0].number_of_inputs == 2)

    integ_funcs = bxb1_rule.integrate(plasb[1],plasp[1]) #two plasmids -> one plasmid b then p
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[0,None,0,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["forward","forward","forward","forward"])
    assert(integ_funcs[0].circular)
    assert(integ_funcs[0].number_of_inputs == 2)

    integ_funcs = bxb1_rule.integrate(genb[1],plasp[1]) #linear with circular
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[0,None,0,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["forward","forward","forward","forward"])
    assert(not integ_funcs[0].circular)
    assert(integ_funcs[0].number_of_inputs == 2)
    #reverse integration
    integ_funcs = bxb1_rule.integrate(plasp[1],genb[1]) #circular with linear
    assert(len(integ_funcs)==1)#one product is made
    assert([a[0].position for a in integ_funcs[0].partslist]==[0,None,0,None])
    assert([a[1] for a in integ_funcs[0].partslist]==["forward","forward","forward","forward"])
    assert(not integ_funcs[0].circular)
    assert(integ_funcs[0].number_of_inputs == 2)
    #recombination
    integ_funcs = bxb1_rule.integrate(genp[1],genb[1]) #linear with linear
    assert(len(integ_funcs)==2) #two products
    assert(sum([a.circular for a in integ_funcs])==0) #all products are linear
    for prod in integ_funcs:
        assert([a[0].position for a in prod.partslist]==[0,None])
        assert([a[1] for a in prod.partslist]==["forward","forward"])
        assert(prod.number_of_inputs == 2)

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
    attP = IntegraseSite("attP","attP",integrase="Bxb1") #the first argument is the name of this attachment site, the second argument is the type of attachment site it is ("attP", "attB", "attL" or "attR") and the integrase denotes which integrase binds to this site (default is "int1")
    attB = IntegraseSite("attB","attB",integrase="Bxb1") #we define two attachment sites, as one site doesn't do anything on its own besides bind integrases

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
    assert myCRN0.species==myCRN0b.species
    assert all([r in myCRN0.reactions for r in myCRN0b.reactions]) and all([r in myCRN0b.reactions for r in myCRN0.reactions])