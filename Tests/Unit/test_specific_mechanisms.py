import pytest
"""
A catch all of tests to ensure different mechanisms behave as desired
"""

from biocrnpyler import Species, Complex

#Test Metabolite Mechanisms
from biocrnpyler import OneStepPathway

def test_OneStepPathway():
    osp = OneStepPathway()
    precursor = [Species("S")]
    product = [Species("P")]

    #Test Update Species
    assert len(osp.update_species(precursor, product)) == 2

    #This mechanism works with None as well
    assert len(osp.update_species(precursor, None)) == 1
    assert len(osp.update_species(None, product)) == 1

    #Test update reactions
    assert len(osp.update_reactions(precursor, product, k = 1)) == 1
    assert len(osp.update_reactions(precursor, None, k = 1)) == 1
    assert len(osp.update_reactions(None, product, k = 1)) == 1

#Test Binding Mechanisms
from biocrnpyler import One_Step_Cooperative_Binding, Two_Step_Cooperative_Binding, One_Step_Binding

def test_One_Step_Binding():
    osb = One_Step_Binding()
    binder = Species("S1")
    bindee = Species("S2")
    c_fake = Species("C")
    c1 = Complex([binder, bindee])


    #Test Update Species
    assert len(osb.update_species(binder, bindee, cooperativity = 1)) == 3
    assert c1 in osb.update_species(binder, bindee, cooperativity = 1)
    assert c_fake in osb.update_species(binder, bindee, complex_species = c_fake)


    #Test Update Reactions
    assert len(osb.update_reactions(binder, bindee, kb = 1.0, ku = 1.0)) == 1
    assert len(osb.update_reactions(binder, bindee, kb = 1.0, ku = 1.0, complex_species = c_fake,)) == 1

    #Also works with lists
    assert len(osb.update_species([binder], [bindee])) == 3
    assert len(osb.update_reactions([binder], [bindee], kb = 1.0, ku = 1.0)) == 1


def test_One_Step_Cooperative_Binding():
    oscb = One_Step_Cooperative_Binding()
    binder = Species("S1")
    bindee = Species("S2")
    c_fake = Species("C")
    c1 = Complex([binder, bindee])
    c2 = Complex([binder, binder, bindee])

    #Test Update Species
    #Try Cooperativity 1
    assert len(oscb.update_species(binder, bindee, cooperativity = 1)) == 3
    assert c1 in oscb.update_species(binder, bindee, cooperativity = 1)
    assert c_fake in oscb.update_species(binder, bindee, cooperativity = 1, complex_species = c_fake)
    #Try Cooperativity 2
    assert len(oscb.update_species(binder, bindee, cooperativity = 2)) == 3
    assert c2 in oscb.update_species(binder, bindee, cooperativity = 2)
    assert c_fake in oscb.update_species(binder, bindee, cooperativity = 2, complex_species = c_fake)


    #Test Update Reactions
    #Try Cooperativity 1
    assert len(oscb.update_reactions(binder, bindee, cooperativity = 1, kb = 1.0, ku = 1.0)) == 1
    assert len(oscb.update_reactions(binder, bindee, cooperativity = 2, kb = 1.0, ku = 1.0, complex_species = c_fake,)) == 1

    #Try Cooperativity 2
    assert len(oscb.update_reactions(binder, bindee, cooperativity = 2, kb = 1.0, ku = 1.0)) == 1
    assert len(oscb.update_reactions(binder, bindee, cooperativity = 2, kb = 1.0, ku = 1.0, complex_species = c_fake,)) == 1

def test_Two_Step_Cooperative_Binding():
    tscb = Two_Step_Cooperative_Binding()
    binder = Species("S1")
    bindee = Species("S2")
    c_fake = Species("C")
    nmer_fake = Species("n")
    c1 = Complex([nmer_fake, bindee])
    c2 = Complex([Complex([binder, binder]), bindee])

    #Test Update Species
    #Try Cooperativity 1

    #The following will fail
    with pytest.raises(ValueError):
        tscb.update_species(binder, bindee, cooperativity = 1)
    with pytest.raises(ValueError):
        tscb.update_species(binder, bindee, cooperativity = 1, complex_species = c_fake)

    #These work due to passing in the nmer
    assert nmer_fake in tscb.update_species(binder, bindee, cooperativity = 1, n_mer_species = nmer_fake)
    assert c1 in tscb.update_species(binder, bindee, cooperativity = 1, n_mer_species = nmer_fake)
    assert len(tscb.update_species(binder, bindee, cooperativity = 1, n_mer_species = nmer_fake)) == 4

    #Try Cooperativity 2
    assert len(tscb.update_species(binder, bindee, cooperativity = 2)) == 4
    assert c2 in tscb.update_species(binder, bindee, cooperativity = 2)
    assert c_fake in tscb.update_species(binder, bindee, cooperativity = 2, complex_species = c_fake)
    assert nmer_fake in tscb.update_species(binder, bindee, cooperativity = 2, n_mer_species = nmer_fake)


    #Test Update Reactions
    assert len(tscb.update_reactions(binder, bindee, cooperativity = 1, kb = (1.0, 1.0), ku = (1.0,1.0), n_mer_species = nmer_fake, complex_species = c_fake)) == 2
    assert len(tscb.update_reactions(binder, bindee, cooperativity = 2, kb = (1.0, 1.0), ku = (1.0, 1.0))) == 2

    #2 rates are required for kb and ku, each
    with pytest.raises(TypeError):
        assert tscb.update_reactions(binder, bindee, cooperativity = 2, kb = 1.0, ku = 1.0)
    