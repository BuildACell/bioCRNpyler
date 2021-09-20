
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
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

    bxb1_enumerator = Integrase_Enumerator("Bxb1", int_mechanisms={"Bxb1":bxb1_rule}) 

    
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
    aL.direction= "forward"
    aR = IntegraseSite("attR","attR",integrase="Bxb1")
    aR.direction = "forward"
    productsites = bxb1_rule.generate_products(delete[0],delete[2])

    testaL = copy.copy(aL)
    testaL.direction = "forward"
    testaL.position = 0
    testaL.parent = delete
    testaR = copy.copy(aR)
    testaR.direction = "forward"
    testaR.position = 2
    testaR.parent = delete
    assert(productsites[0]==testaL) #generated product sites have the right parent, position, and direction
    assert(productsites[1]==testaR) #generated product sites have the right parent, position, and direction
    

    ap = IntegraseSite("attP","attP",integrase="Bxb1")
    ab = IntegraseSite("attB","attB",integrase="Bxb1")
    delete = DNA_construct([ab,cds,ap])
    flip = DNA_construct([ab,cds,[ap,"reverse"]],circular=True)
    flipr = flip.get_reversed()
    plasp = DNA_construct([cds,ap],circular=True)
    plasb = DNA_construct([cds,ab],circular=True)
    genp = DNA_construct([cds2,ap])
    genb = DNA_construct([cds2,ab])
    #testing the actual integration reactions
    #inversion
    new_constructs = bxb1_rule.integrate(flip[0],flip[2],also_inter=False) #sites forwards
    assert(len(new_constructs)==1)#one product is made
    assert([part.name for part in new_constructs[0]]==["attL","GFP","attR"])
    assert([part.direction for part in new_constructs[0]]==["forward","reverse","reverse"])

    new_constructs = bxb1_rule.integrate(flip[2],flip[0],also_inter=False) #sites reverse
    assert(len(new_constructs)==1)#one product is made
    assert([part.name for part in new_constructs[0]]==["attL","GFP","attR"])
    assert([part.direction for part in new_constructs[0]]==["forward","reverse","reverse"])

    new_constructs = bxb1_rule.integrate(flipr[2],flipr[0],also_inter=False) #sites forwards
    assert(len(new_constructs)==1)#one product is made
    assert([part.name for part in new_constructs[0]]==["attR","GFP","attL"])
    assert([part.direction for part in new_constructs[0]]==["forward","forward","reverse"])
    
    new_constructs = bxb1_rule.integrate(flip[0],flip[2],also_inter=True) #allow intermolecular reactions
    assert(len(new_constructs)==2)#now, two products are made
    assert(sum([a.circular for a in new_constructs])==2) #all of them is circular
    for prod in new_constructs:
        if(len(prod)==3):
            assert([part.name for part in prod]==["attL","GFP","attR"])
            assert([part.direction for part in prod]==["forward","reverse","reverse"])
        else:
            assert([part.name for part in prod]==['attL', 'GFP', 'attB', 'attR', 'GFP', 'attP'])
            assert([part.direction for part in prod]==["forward","reverse","reverse","forward","forward","reverse"])
    #deletion
    new_constructs = bxb1_rule.integrate(delete[0],delete[2],also_inter=False) #sites forwards
    assert(len(new_constructs)==2)#two product is made
    assert(sum([a.circular for a in new_constructs])==1) #one of them is circular
    for prod in new_constructs:
        if(prod.circular):
            assert([part.name for part in prod]==["attR","GFP"])
            assert([part.direction for part in prod]==["forward","forward"])
        else:
            assert([part.name for part in prod]==["attL"])
            assert([part.direction for part in prod]==["forward"])
    #now, what if we also calculate the intermolecular reactions?
    new_constructs = bxb1_rule.integrate(delete[2],delete[0],also_inter=True) #sites reverse
    assert(len(new_constructs)==3)#four product is made
    assert(sum([a.circular for a in new_constructs])==1) #one of them is circular
    for prod in new_constructs:
        if(prod.circular):
            if(len(prod)==2):
                print(prod)
                assert([part.name for part in prod]==["attR","GFP"])
                assert([part.direction for part in prod]==["forward","forward"])
        elif(len(prod)==1):
            assert([part.name for part in prod]==["attL"])
            assert([part.direction for part in prod]==["forward"])
        elif(len(prod)==5):
            assert([part.name for part in prod]==["attB","GFP","attR","GFP","attP"])
            assert([part.direction for part in prod]==["forward"]*5)


    #integration
    new_constructs = bxb1_rule.integrate(plasp[1],plasb[1]) #two plasmids -> one plasmid p then b
    assert(len(new_constructs)==1)#one product is made
    prod = new_constructs[0]
    assert([part.name for part in prod]==["GFP","attR","GFP","attL"])
    assert([part.direction for part in prod]==["forward","forward","forward","forward"])
    assert(new_constructs[0].circular)

    new_constructs = bxb1_rule.integrate(plasb[1],plasp[1]) #two plasmids -> one plasmid b then p
    assert(len(new_constructs)==1)#one product is made
    prod = new_constructs[0]
    assert([part.name for part in prod]==["GFP","attL","GFP","attR"])
    assert([part.direction for part in prod]==["forward","forward","forward","forward"])
    assert(new_constructs[0].circular)


    bxb1_enumerator.reset([plasb,plasp])
    assert(len(plasb[1].linked_sites)==0) #linked sites are nonexistant if it is reset
    assert(len(plasp[1].linked_sites)==0) #linked sites are nonexistant if it is reset

    #include result construct
    result_construct = DNA_construct([cds,aL,cds,aR],circular=True)
    new_constructs = bxb1_rule.integrate(plasb[1],plasp[1], existing_dna_constructs=[result_construct]) #two plasmids -> one plasmid b then p
    assert(len(new_constructs)==0)#no new products are made, only the existing one
    assert((plasp[1],True) in plasb[1].linked_sites) #proper reaction found inside integrase site
    assert((plasb[1],True) in plasp[1].linked_sites ) #proper reaction found inside integrase site

    #including circularly permuted existing construct
    #reaction is:
    #[cds,ab] + [cds,ap] =[cds, aL, cds, aR]
    bxb1_enumerator.reset([plasb,plasp])
    permuted_construct = DNA_construct([cds,aR,cds,aL],circular=True)
    new_constructs = bxb1_rule.integrate(plasb[1],plasp[1], existing_dna_constructs=[permuted_construct]) #two plasmids -> one plasmid b then p
    assert(len(new_constructs)==0)#no new products are made, only the existing one
    assert((plasp[1],True) in plasb[1].linked_sites ) #proper reaction found inside integrase site
    assert((plasb[1],True) in plasp[1].linked_sites ) #proper reaction found inside integrase site


    bxb1_enumerator.reset([plasb,plasp])
    new_constructs = bxb1_rule.integrate(genb[1],plasp[1]) #linear with circular
    assert(len(new_constructs)==1)#one product is made
    prod = new_constructs[0]
    assert([part.name for part in prod]==['RFP', 'attL', 'GFP', 'attR'])
    assert([part.direction for part in prod]==["forward","forward","forward","forward"])
    assert(not new_constructs[0].circular)


    #reverse integration
    bxb1_enumerator.reset([plasb,plasp])
    new_constructs = bxb1_rule.integrate(plasp[1],genb[1]) #circular with linear
    assert(len(new_constructs)==1)#one product is made
    prod = new_constructs[0]
    assert([part.name for part in prod]==['RFP', 'attL', 'GFP', 'attR'])
    assert([part.direction for part in prod]==["forward","forward","forward","forward"])
    assert(not new_constructs[0].circular)
    #recombination
    bxb1_enumerator.reset([plasb,plasp])
    new_constructs = bxb1_rule.integrate(genp[1],genb[1]) #linear with linear
    assert(len(new_constructs)==2) #two products
    assert(sum([a.circular for a in new_constructs])==0) #all products are linear
    for prod in new_constructs:
        assert([part.name for part in prod] in [['RFP', 'attL'],['RFP', 'attR']])
        assert([part.direction for part in prod]==["forward","forward"])

    #

def test_integrase_enumerator():
    cds = CDS("GFP")
    cds2 = CDS("RFP")
    ap = IntegraseSite("attP","attP",integrase="Bxb1")
    ab = IntegraseSite("attB","attB",integrase="Bxb1")
    al = IntegraseSite("attL","attL",integrase="Bxb1")
    ar = IntegraseSite("attR","attR",integrase="Bxb1")


    flip = DNA_construct([ab,cds,[ap,"reverse"]])
    plasp = DNA_construct([cds,ap],circular=True)
    plasb = DNA_construct([cds,ab],circular=True)



    genp = DNA_construct([cds2,ap])
    genb = DNA_construct([cds2,ab])
    


    bxb1_rule = IntegraseRule(name="Bxb1",reactions={("attB","attP"):"attL",("attP","attB"):"attR"})

    bxb1_enumerator = Integrase_Enumerator("Bxb1", int_mechanisms={"Bxb1":bxb1_rule}) 

    int_dict = bxb1_enumerator.list_integrase(flip)
    assert("Bxb1" in int_dict) #the correct integrase is listed
    assert(flip[0] in int_dict["Bxb1"]) #both sites are found in the dictionary
    assert(flip[2] in int_dict["Bxb1"])  #both sites are found in the dictionary

    x = bxb1_rule.integrate(flip[0],flip[2])
    assert((flip[2],False) in flip[0].linked_sites) #sites are properly integrated

    bxb1_enumerator.reset([flip])
    assert(len(flip[0].linked_sites)==0)
    assert(len(flip[2].linked_sites)==0)

    permuted_plasp = DNA_construct([ap,cds],circular=True)
    assert(plasp != permuted_plasp) #normal equality doesn't cut it
    found,perm_func = bxb1_enumerator.find_dna_construct(plasp,[plasp])
    assert([perm_func(a) for a in range(len(plasp))] == [(lambda a: (a,"f"))(a) for a in range(len(plasp))]) #1:1 find matching
    assert(found == plasp) #found the same construct we put in

    found,perm_func = bxb1_enumerator.find_dna_construct(plasp,[permuted_plasp])
    assert([perm_func(a) for a in range(len(plasp))] == [(1,"f"),(0,"f")]) #find permuted mapping
    assert(found == permuted_plasp) #found the same construct we put in
    assert(bxb1_enumerator.find_dna_construct(plasp,[plasb]) is None) #we don't find anything

    reversed_plasp = DNA_construct([(ap,"reverse"),(cds,"reverse")],circular=True)
    reversed_permuted_plasp = DNA_construct([(cds,"reverse"),(ap,"reverse")],circular=True)
    found,perm_func = bxb1_enumerator.find_dna_construct(plasp,[reversed_plasp])
    assert(found == reversed_plasp) #found the same construct we put in
    assert([perm_func(a) for a in range(len(plasp))] == [(1,"r"),(0,"r")]) #find permuted mapping

    found,perm_func = bxb1_enumerator.find_dna_construct(plasp,[reversed_permuted_plasp])
    assert(found == reversed_permuted_plasp) #found the same construct we put in
    assert([perm_func(a) for a in range(len(plasp))] == [(0,"r"),(1,"r")]) #find permuted mapping

    reversed_genp = DNA_construct([(ap,"reverse"),(cds2,"reverse")])
    found,perm_func = bxb1_enumerator.find_dna_construct(genp,[reversed_genp])
    assert(found == reversed_genp) #found the same construct we put in
    assert([perm_func(a) for a in range(len(genp))] == [(1,"r"),(0,"r")]) #find permuted mapping

    plaspb = DNA_construct([cds,ar,cds,al],circular=True)
    x = bxb1_enumerator.enumerate_components([plasp,plasb])
    assert(plaspb in x) #the proper construct is made

    y = bxb1_enumerator.enumerate_components([plasp,plasb],previously_enumerated=[plaspb])
    assert (y == []) #nothing new is made since the proper construct has been "previously enumerated"

    plaspb = DNA_construct([ar,cds,al,cds,ar,cds,ar,cds],circular=True)
    for i in range(1,len(plaspb)):
        #test all circular permutations of this repetitive plasmid. They should all  yield the same hash!
        plasperm = plaspb.get_circularly_permuted(i)
        h1 = DNA_construct.omnihash(plaspb)
        h2 = DNA_construct.omnihash(plasperm)
        assert(h1[0]==h2[0])



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
    assert len(comps1) == 85
    assert len(myCRN1.species) == 52
    assert len(myCRN1.reactions) == 97

    #Recompiling at a different length doesn't change anything
    myCRN0b, comps0b = myMixture.compile_crn(recursion_depth = 0, 
                                         return_enumerated_components = True,
                                         initial_concentrations_at_end = True,
                                         copy_objects = False,
                                         add_reaction_species = False)
    assert len(comps0) == len(comps0b)
    assert myCRN0.species==myCRN0b.species
    assert all([r in myCRN0.reactions for r in myCRN0b.reactions]) and all([r in myCRN0b.reactions for r in myCRN0.reactions])