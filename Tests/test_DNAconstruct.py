
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from unittest import TestCase
from biocrnpyler import *
import copy

from biocrnpyler import Transcription_MM, Translation_MM, Degredation_mRNA_MM

class TestDNAConstruct(TestCase):
    def test_initialization(self):
        myprom = Promoter("prom")
        myrbs = RBS("rbs")
        mycds = CDS("mycds","GFP")
        myt = Terminator("t1")
        myz = AttachmentSite("attB","attB",direction="forward",color="blue",color2="skyblue",sequence="agcgcga")
        self.assertTrue("attB" in str(myz))

        myconst = DNA_construct([[myprom,"forward"],\
                                 [myrbs,"forward"],\
                                 [mycds,"forward"],\
                                 [myt,"forward"]])

        bspec = Species("prom_rbs_mycds_t1",material_type="dna")
        #base species naming
        self.assertEqual(myconst.base_species,bspec)
        #Components are contained
        self.assertTrue(myprom in myconst)
        self.assertTrue(myrbs in myconst)
        self.assertTrue(mycds in myconst)
        self.assertTrue(myt in myconst)
        self.assertTrue("prom" in myconst)

        spec = myconst.get_species()
        p = copy.deepcopy(spec[0])
        p.remove()
        u = copy.deepcopy(spec[1])
        u.remove()
        c = copy.deepcopy(spec[2])
        c.remove()
        t = copy.deepcopy(spec[3])
        t.remove()
        self.assertTrue(p == Species(myprom.name,material_type="dna"))
        self.assertTrue(u == Species(myrbs.name,material_type="dna"))
        self.assertTrue(c == Species(mycds.name,material_type="dna"))
        self.assertTrue(t == Species(myt.name,material_type="dna"))
        for a in spec:
            self.assertTrue(a.direction=="forward")
        myconst2 = DNA_construct([[myt,"reverse"],\
                                 [mycds,"reverse"],\
                                 [myrbs,"reverse"],\
                                 [myprom,"reverse"]],circular=True)
        self.assertTrue(myconst2.parts_list[0].direction=="reverse")
        self.assertTrue(myconst2.parts_list[1].direction=="reverse")
        self.assertTrue(myconst2.parts_list[2].direction=="reverse")
        self.assertTrue(myconst2.parts_list[3].direction=="reverse")
        myconst3 = copy.deepcopy(myconst)
        myconst3.circular = True
        myconst3.reverse()
        self.assertTrue(myconst3==myconst2)

        #testing automatic name creation where we insert _r for reverse and _o for circular
        self.assertTrue("_r_" in myconst2.name)
        self.assertTrue("_o" in myconst2.name)
    def test_txtl(self):
        myprom = Promoter("prom")
        myrbs = RBS("rbs")
        myrbs2 = RBS("rbs2")
        mycds = CDS("mycds","GFP")
        mycds2 = CDS("mycds2","RFP")
        myt = Terminator("t1")
        mycdsNOSTOP = CDS("mycds3","CFP",no_stop_codons=["forward"])

        myconst = DNA_construct([[myrbs,"forward"],\
                                 [myrbs,"forward"],\
                                 [mycdsNOSTOP, "forward"],\
                                 [mycds,"forward"],\
                                 [myt,"forward"], \
                                 [myprom,"forward"],\
                                 [myrbs2,"forward"],\
                                 [mycds2,"forward"]],circular=True)
        rnas,proteins = myconst.explore_txtl()
        #circular construct properly makes RNA
        self.assertTrue(myconst[5] in rnas)

        myTranscript = rnas[myconst[5]]
        myRBS = myTranscript[0]
        myRBSnotl = myTranscript[2]
        myRBS2 = myTranscript[3]
        myCDS = myTranscript[1]
        myCDS2 = myTranscript[4]
        myCDS3 = myTranscript[5]
        #transcript is made
        self.assertTrue(myTranscript in proteins)
        #two RBSes found in the transcript
        self.assertTrue(myRBS in proteins[myTranscript])
        self.assertTrue(myRBSnotl in proteins[myTranscript])
        self.assertTrue(myRBS2 in proteins[myTranscript])
        #proper protein is made by the right RBS
        self.assertTrue(myCDS in proteins[myTranscript][myRBS])
        self.assertTrue(myCDS2 in proteins[myTranscript][myRBS2])
        self.assertTrue(myCDS3 in proteins[myTranscript][myRBS2])
        self.assertTrue(proteins[myTranscript][myRBSnotl]==[])
        #CDSes have the right name
        self.assertTrue(proteins[myTranscript][myRBS][0].name==mycds2.name)
        self.assertTrue(proteins[myTranscript][myRBS2][1].name==mycds.name)
    def test_dna_construct_species(self):
        myprom = Promoter("prom")
        myrbs = RBS("rbs")
        myrbs2 = RBS("rbs2")
        mycds = CDS("mycds","GFP")
        mycds2 = CDS("mycds2","RFP")
        myt = Terminator("t1")
        myconst = DNA_construct([[myrbs,"forward"],\
                                 [mycds,"forward"],\
                                 [myt,"forward"], \
                                 [myprom,"forward"],\
                                 [myrbs2,"forward"],\
                                 [mycds2,"forward"]], \
                                 circular=True)
        parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
        myMixture = TxTlExtract(name = "txtl", parameters = parameters, \
                            components = [myconst])
        myCRN = myMixture.compile_crn()
        myspec = myCRN.species
        myrxn = myCRN.reactions
        #generating the right number of species and reactions
        self.assertTrue(len(myspec)==11)
        self.assertTrue(len(myrxn)==10)
        #one promoter makes one rna
        self.assertTrue(sum([spec.material_type=="rna" for spec in myspec])==1)

        newplist = list(myconst.parts_list[:3])+[[myprom,"reverse"]]+list(myconst.parts_list[3:])
        myconst2 = DNA_construct(newplist,circular=True)
        myMixture2 = TxTlExtract(name = "txtl", parameters = parameters, \
                            components = [myconst2])
        myCRN2 = myMixture2.compile_crn()
        #we add a promoter and more reactions and species are made
        self.assertTrue(len(myCRN2.species)==14)
        self.assertTrue(len(myCRN2.reactions)==16)
        #now there are two promoters
        self.assertTrue(sum([spec.material_type=="rna" for spec in myCRN2.species])==2)
        





