#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import re
import warnings

import pytest

import biocrnpyler as bcp
from biocrnpyler.utils.plotting import CRNPlotter
from biocrnpyler.components import RNA_construct, DNA_construct,\
    Promoter, IntegraseSite, RBS, CDS, Terminator, Complex,\
    Species, DNA_part, RegulatedPromoter


def test_CRNPlotter():
    class dummy_renderer:
        class dummy_parts:
            def __init__(self):
                pass
        def __init__(self,scale=1,linewidth=2):
            self.linecolor = (0,0,0)
            pass
        def SBOL_part_renderers(self):
            return self.dummy_parts()
        def std_reg_renderers(self):
            return self.dummy_parts()
        def renderDNA(self,ax,designs,renderers,**keywords):
            ax.plot([0],[1])

            return(0,100)

    d_r = dummy_renderer()
    test_plotter = CRNPlotter(dna_renderer=dummy_renderer(),rna_renderer=dummy_renderer(),cmap=[0,1,2,3,4,5,6])

    test_construct = DNA_construct([Promoter("p1"),IntegraseSite("attP","attP"),Terminator("t1")])
    test_construct2 = DNA_construct([RBS("utr1"),IntegraseSite("attB","attB"),Terminator("t1")])
    test_construct3 = DNA_construct([IntegraseSite("attL","attL"),IntegraseSite("attR","attR"),CDS("o1")])
    ts1 = test_construct.get_species()
    test_species = Complex([ts1[1],Species("integrase",material_type="protein")]).parent


    assert(test_plotter.part_dpl_dict == {})
    assert(test_plotter.construct_dpl_dict == {})
    assert(test_plotter.species_dpl_dict == {})
    assert(test_plotter.species_image_dict == {})

    test_plotter.make_dpls_from_construct(test_construct)
    assert(len(test_plotter.part_dpl_dict)==len(test_construct)) #new parts are added to the dictionary
    assert(len(test_plotter.construct_dpl_dict) ==1) #the new construct is added
    #part added to part dict
    assert(test_plotter.part_dpl_dict[Promoter("p1")].get_dpl()[0]['type']=="Promoter")
    #construct added to construct dict
    condpl = test_plotter.construct_dpl_dict[test_construct].get_dpl()
    #test that everything has the right type
    assert(condpl[0]['type']=="Promoter")
    assert(condpl[1]['type']=="RecombinaseSite")
    assert(condpl[2]['type']=="Terminator")


    test_plotter.make_dpls_from_construct(test_construct2)
    assert(len(test_plotter.part_dpl_dict)==len(test_construct)+len(test_construct2)-1)
    #new parts from the new construct are added but the terminator is the same so don't add that one
    assert(len(test_plotter.construct_dpl_dict) ==2) #now there are two constructs
    test_plotter.make_dpls_from_construct(test_construct3)
    #almost all new parts are added
    assert(len(test_plotter.part_dpl_dict)==len(test_construct)+len(test_construct2)+len(test_construct3)-1)
    assert(len(test_plotter.construct_dpl_dict) ==3) #there are now 3 constructs!


    test_plotter.make_dpl_from_species(test_species) #making dpl from single species

    test_plotter.make_dpl_from_species(Species("integrase",material_type="protein"))
    assert(Species("integrase",material_type="protein") in test_plotter.species_dpl_dict)

    test_plotter.make_dpl_from_part(RegulatedPromoter("p1",regulators=["r1","r2"]))
    assert(RegulatedPromoter("p1",regulators=["r1","r2"]) in test_plotter.part_dpl_dict) #make a multipart promoter
    #normal promoter is a multipart with a promoter being first
    assert(test_plotter.part_dpl_dict[RegulatedPromoter("p1",regulators=["r1","r2"])].get_dpl()[0]["type"]=="Promoter")
    #try returning a reverse version of the same thing
    boundspec = test_plotter.make_dpl_from_species(Species("ooga",material_type="booga")) #make a species to be bound to the promoter
    revpart = test_plotter.part_dpl_dict[RegulatedPromoter("p1",regulators=["r1","r2"])].get_directed("reverse",bound=[boundspec]) #return reverse promoter with things bound
    assert(revpart.get_dpl()[0]["type"]=="Operator") #reverse promoter has operators on the left


    test_plotter.colordict = {"Promoter":"purple"}
    prom_dpl = test_plotter.make_dpl_from_part(Promoter("pconst")) #make a promoter part
    assert(prom_dpl.color == "purple") #found the proper color
    prom_dpl2 = test_plotter.make_dpl_from_part(Promoter("funkbeast")) #make a different promoter part
    assert(prom_dpl2.color=="purple")

    rnadpl = test_plotter.make_dpl_from_species(Species("testrna",material_type="rna"))
    assert(rnadpl.dpl_type == "NCRNA")
    #something is bound to the operator closest to the reversed promoter
    assert(revpart.parts_list[1].bound is not None)

    #construct having a multi-part promoter and a "user defined" part
    test_construct4 = DNA_construct([RegulatedPromoter("p1",regulators=["r1","r2"]),DNA_part("test")])
    test_plotter.make_dpl_from_species(test_construct4.get_species())
    assert(test_construct4.get_species() in test_plotter.species_dpl_dict) #construct is in the dictionary

    rna_construct = RNA_construct([RBS("ooga"),CDS("booga")])
    condpl1 = test_plotter.make_dpls_from_construct(rna_construct)
    assert(condpl1.material_type == "rna") #material type is set correctly

    rna_construct2 = RNA_construct([RBS("ooga"),CDS("booga")])
    condpl2 = test_plotter.make_dpls_from_construct(rna_construct2)
    assert(condpl2==condpl1) #this construct is already in the dictionary! Just return it



    bound_construct = Complex([test_construct4.get_species()[1],Species("integrase",material_type="protein")]).parent
    #making a construct that is bound to something
    dpl_binders = test_plotter.make_dpl_from_species(bound_construct).get_dpl_binders()
    #return the things that are bound
    assert(dpl_binders[0][0]["type"]=="Binding")


    test_plotter.make_dpl_from_species(Complex([Species("mydna",material_type="dna"),
                                                Species("integrase",material_type="protein")]))
    #dpl made from a complex
    #make sure the part inside the construct made it into the species dict
    assert(Species("mydna",material_type="dna") in test_plotter.species_dpl_dict)


def test_render_network_bokeh():
    # Create a simple reaction
    parameters = {
        (None, None, "ktx"): 0.05,
        (None, None, "kb"): 100,
        (None, None, "ku"): 10,
        (None, None, "ktl"): 0.05,
        ("rna_degredation_mm", None, "kdeg"): 0.001,
        (None, None, "cooperativity"): 2}
    txtl = bcp.TxTlExtract('mixture1', parameters=parameters)
    dna = bcp.DNAassembly(
        'mydna', promoter=bcp.RegulatedPromoter('plac', ['laci']),
        rbs='UTR1', protein='GFP', initial_concentration=10)
    txtl.add_component(dna)
    crn1 = txtl.compile_crn()
    crn1.add_reactions([
        bcp.Reaction.from_massaction(
            [bcp.Species('mydna', material_type='rna')], [], k_forward=0.1)])

    # Generate a graph and make sure it is created
    with warnings.catch_warnings(record=True) as records:
        plot = bcp.utils.plotting.render_network_bokeh(crn1)

    for w in records:
        if re.search("plotting disabled", w.message.args[0]):
            pytest.fail("plotting library import error")
        else:
            warnings.showwarning(
                w.message.args[0], w.category, w.filename, w.lineno)
