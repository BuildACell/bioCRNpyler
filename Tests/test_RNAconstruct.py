#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
#from unittest import TestCase
from biocrnpyler import *
import copy


def basic_RNAconstruct():
    R = RBS("rbs")
    C = CDS("gfp")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"translation":Translation_MM(Species("Ribo"))}

    #minimal RNA translation
    x = RNA_construct([R],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(not (R in y)) #RBS is copied correctly
    assert([x[0]] == y) #RBS with no protein is still active
    assert(y[0].protein is None) #protein is nothing

    x = RNA_construct([[R,"reverse"]],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y == []) #reverse RBS is not an RBS
    assert(x[0].protein is None) #protein is nothing

    #just a cds
    x = RNA_construct([C],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y==[]) #just a CDS does nothing

def basic_RBS_CDS_RNAconstruct():
    R = RBS("rbs")
    C = CDS("gfp")
    parameters={"cooperativity":2,"kb":100, "ku":10, "ktx":.05, "ktl":.2, "kdeg":2,"kint":.05}
    mechs = {"translation":Translation_MM(Species("Ribo"))}

    #minimal RNA translation
    x = RNA_construct([R,C],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(x[0] ==y[0])
    assert(y[1].protein==Species("gfp",material_type="protein"))
    assert(y[0].transcript.parent==x.get_species())

    #minimal RNA transcription
    x = RNA_construct([C,R],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y==[x[1]])
    assert(y[0].transcript.parent==x.get_species())

    #minimal RNA transcription
    x = RNA_construct([[R,"reverse"],C],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y==[])

    #minimal RNA transcription
    x = RNA_construct([R,[C,"reverse"]],mechanisms = mechs,parameters=parameters)
    y = x.enumerate_components()
    assert(y==[x[0]])
    assert(y[0].transcript.parent==x.get_species())