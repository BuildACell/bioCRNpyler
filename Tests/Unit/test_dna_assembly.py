import pytest
from biocrnpyler import DNAassembly, Promoter, RBS, Mechanism, Species


def test_instantiation():
    #Instantiate from strings
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom")
    assert type(ass.promoter) is Promoter and ass.promoter.name == "prom"
    assert type(ass.rbs) is RBS and ass.rbs.name == "rbs"
    assert type(ass.transcript) is Species and ass.transcript.name == ass.name
    assert type(ass.protein) is Species and ass.protein.name == ass.name
    assert type(ass.dna) is Species and ass.dna.name == ass.name

    #Instantiate with rbs = None
    ass = DNAassembly(name = "ass", rbs = None, promoter = "prom")
    assert ass.rbs is None

    #Instantiate with prom = None
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = None)
    assert ass.promoter is None
    assert type(ass.rbs) is RBS and ass.rbs.name == "rbs"
    assert type(ass.protein) is Species and ass.protein.name == ass.name

    #instantiate with transcript = Species
    T = Species("T")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom", transcript = T)
    assert ass.transcript == T

    #instantiate with dna = Species
    D = Species("D")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom", dna = D)
    assert ass.dna == D

    #instantiate with protein = Species
    P = Species("P")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom", protein = P)
    assert ass.protein == P


def test_update_promoter():
    T = Species("T")
    D = Species("D")
    P = Species("P")
    prom = Promoter("prom")
    mech = Mechanism(name = "mech", mechanism_type = "assembly")

    ass = DNAassembly("ass", rbs = None, promoter = None, transcript = T, dna = D, protein = P, mechanisms = [mech])
    ass.update_promoter(prom)

    #test copying
    assert ass.promoter != prom
    assert type(ass.promoter) == type(prom)
    assert ass.promoter.name == prom.name

    #test resetting transcript
    assert ass.promoter.transcript == T

    #test mechanism inheritance
    assert mech.mechanism_type in ass.promoter.mechanisms

def test_update_rbs():
    T = Species("T")
    D = Species("D")
    P = Species("P")
    rbs = RBS("rbs")
    mech = Mechanism(name = "mech", mechanism_type = "assembly")

    ass = DNAassembly("ass", rbs = None, promoter = None, transcript = T, dna = D, protein = P, mechanisms = [mech])
    ass.update_rbs(rbs)

    #test copying
    assert ass.rbs != rbs
    assert type(ass.rbs) == type(rbs)
    assert ass.rbs.name == rbs.name

    #test resetting transcript
    assert ass.rbs.protein == P

    #test mechanism inheritance
    assert mech.mechanism_type in ass.rbs.mechanisms

def test_update_transcript():
    T = Species("T")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom")
    ass.update_transcript(T)
    assert ass.transcript == T
    assert ass.promoter.transcript == T
    assert ass.rbs.transcript == T

def test_update_dna():
    D = Species("D")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom")
    ass.update_dna(D)
    assert ass.dna == D

def test_update_protein():
    P = Species("P")

    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom")
    ass.update_protein(P)
    assert ass.protein == P
    assert ass.promoter.protein == P
    assert ass.rbs.protein == P

def test_add_mechanism():
    mech = Mechanism(name = "mech", mechanism_type = "assembly")
    ass = DNAassembly(name = "ass", rbs = "rbs", promoter = "prom")
    ass.add_mechanism(mech)
    #Mechanism should be added to promoter and rbs
    assert mech.mechanism_type in ass.mechanisms
    assert mech.mechanism_type in ass.promoter.mechanisms
    assert mech.mechanism_type in ass.rbs.mechanisms

    #Mechaisms which are already in the Promoter or RBS should NOT be overwritten
    mech_rbs = Mechanism(name = "rbs", mechanism_type = "mech")
    mech_promoter = Mechanism(name = "promoter", mechanism_type = "mech")
    mech_assembly = Mechanism(name = "assembly", mechanism_type = "mech")

    #Create promoter and assembly with their own mechanisms
    prom = Promoter("prom", mechanisms = [mech_promoter])
    assert mech_promoter.mechanism_type in prom.mechanisms
    rbs = RBS("rbs", mechanisms = [mech_rbs])
    assert mech_rbs.mechanism_type in rbs.mechanisms
    ass = DNAassembly(name = "ass", rbs = rbs, promoter = prom)
    assert ass.rbs.mechanisms["mech"].name == mech_rbs.name
    assert ass.promoter.mechanisms["mech"].name == mech_promoter.name

    #Overwrite the assembly mechanism
    #this should not propogate down to rbs or promoter!
    ass.add_mechanism(mech_assembly)
    assert mech_assembly.mechanism_type in ass.mechanisms
    assert ass.rbs.mechanisms["mech"].name == mech_rbs.name
    assert ass.promoter.mechanisms["mech"].name == mech_promoter.name

