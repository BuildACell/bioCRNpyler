#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Enzyme, DNA, RNA, Protein, Metabolite, ChemicalComplex, Species, Complex, LengthComponent
import pytest


def test_enzyme():
    substrates = 'S1'
    products = 'P1'
    enzyme = 'E1'

    e = Enzyme(enzyme=enzyme, substrates=substrates, products=products)
    assert any([s.name == 'S1' for s in e.substrates])
    assert any([p.name == 'P1' for p in e.products])
    assert enzyme == e.enzyme.name

    assert e.get_species().name == 'E1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        e.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        e.update_reactions()

def test_length_component():
    c0 = LengthComponent(name = "C")
    assert c0.length == 0
    c1 = LengthComponent(name = "C", length = 10)
    assert c1.length == 10

    #Some lengths are not allowed
    with pytest.raises(ValueError):
        LengthComponent(name = "C",length = -1)

    with pytest.raises(ValueError):
        LengthComponent(name = "C",length = 0.1)

    with pytest.raises(ValueError):
        LengthComponent(name = "C",length = "10")

def test_DNA():
    dna = DNA(name = "G", length = 10)

    assert dna.get_species() == Species("G", material_type = "dna")
    assert len(dna.update_reactions()) == 0
    assert len(dna.update_species()) == 1
    assert dna.length == 10

def test_RNA():
    rna = RNA(name = "T", length = 10)

    assert rna.get_species() == Species("T", material_type = "rna")
    assert len(rna.update_reactions()) == 0
    assert len(rna.update_species()) == 1
    assert rna.length == 10

def test_Protein():
    protein = Protein(name = "X", length = 10)

    assert protein.get_species() == Species("X", material_type = "protein")
    assert len(protein.update_reactions()) == 0
    assert len(protein.update_species()) == 1
    assert protein.length == 10

def test_Metabolite():
    M = Metabolite(name = "m")
    assert M.get_species() == Species("m", material_type = "metabolite")
    assert len(M.update_reactions()) == 0
    assert len(M.update_species()) == 1

def test_ChemicalComplex():
    #Default Name
    C1 = ChemicalComplex([Species("S1"), Species("S2")])
    assert C1.get_species() == Complex([Species("S1"), Species("S2")])

    #Custom Name
    C2 = ChemicalComplex([Species("S1"), Species("S2")], name = "C2")
    assert C2.get_species() == Complex([Species("S1"), Species("S2")], name = "C2")


    with pytest.raises(KeyError, match='Unable to find mechanism of type binding in Component'):
        C1.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type binding in Component'):
        C1.update_reactions()