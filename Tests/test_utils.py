#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Promoter, IntegraseSite, CDS
from biocrnpyler.utils import member_dictionary_search, process_initial_concentration_dict

def test_process_initial_concentration_dict():
    A = Species('A')
    B = Species('B')
    C = Species('C')
    initial_concentration = {A: 10, B: 25.34, str(C): 77.24}
    processed = process_initial_concentration_dict(initial_concentration)
    for (key_old, key_new) in zip(initial_concentration.keys(), processed.keys()):
        assert str(key_old) == key_new
def test_dictionary_search():
    colordict = {"p1":"blue","Promoter":"red","Bxb1":"orange",\
                    ("Bxb1","attP"):"mauve","attP":"green",\
                    "protein":"purple",("phosphorylated",):"yellow"}


    prom = Promoter("p1")
    prom2 = Promoter("p2")
    c1 = CDS("ooga")
    prot = Species("myprot",material_type="protein")
    rna1 = Species("myrna",material_type="rna")
    rna1.add_attribute("phosphorylated")
    prot2 = Species("myprot2",material_type="protein")
    prot2.add_attribute("phosphorylated")
    int1 = Species("Bxb1","protein")
    int2 = IntegraseSite("attP","attP","Int1")
    int3 = IntegraseSite("attP","attP","Bxb1")
    int4 = IntegraseSite("attP2","attP","Bxb1")
    int5 = IntegraseSite("attP3e","attP")
    print(int4.name)
    print(repr(int4))
    print((int4.integrase,int4.site_type))
    assert(member_dictionary_search(prom,colordict)=="blue") #DNA_part name
    assert(member_dictionary_search(prom2,colordict)=="red") #DNA_part type
    assert(member_dictionary_search(c1,colordict)==None) #this one is not found
    assert(member_dictionary_search(prot,colordict)=="purple") #material type
    assert(member_dictionary_search(rna1,colordict)=="yellow") #attribute is colored
    assert(member_dictionary_search(prot2,colordict)=="purple") #attribute is colored
    colordict[("protein",("phosphorylated",))]="darkblue"
    assert(member_dictionary_search(prot2,colordict)=="darkblue") #(material_type,attribute) is colored differntly
    assert(member_dictionary_search(int1,colordict)=="orange") #integrase protein is colored
    assert(member_dictionary_search(int5,colordict)=="green") #site type in the dictionary
    assert(member_dictionary_search(int3,colordict)=="green") #name of attP is more important than tuple of integrase, attP type
    assert(member_dictionary_search(int4,colordict)=="mauve") #now integrase and type are supreme, since name doesn't match