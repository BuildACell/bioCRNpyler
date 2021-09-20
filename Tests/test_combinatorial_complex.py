#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import CombinatorialComplex, One_Step_Binding, Reaction, CombinatorialConformation
from biocrnpyler import Species, OrderedPolymerSpecies, Complex, PolymerConformation
import pytest

def test_CombinatorialComplex_init_and_properties():
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    
    #Single final_state case
    C1 = Complex([Y, Z, Z])
    CC1 = CombinatorialComplex(final_states = [C1])

    #tests getters and setters
    assert set(CC1.final_states) == set([C1])
    assert set(CC1.sub_species) == set([Y, Z])
    assert set(CC1.initial_states) == set([Y, Z])
    assert CC1.intermediate_states is None

    #multiple final_state case
    C2 = Complex([X, X, Z, Z])
    CC2 = CombinatorialComplex(final_states = [C1, C2])

    #tests getters and setters
    assert set(CC2.final_states) == set([C1, C2])
    assert set(CC2.sub_species) == set([X, Y, Z])
    assert set(CC2.initial_states) == set([X, Y, Z])
    assert CC2.intermediate_states is None

    #test with initial states
    CC3 = CombinatorialComplex(final_states = [C1, C2], initial_states = [Complex([Z, Z]), X])
    assert set(CC3.final_states) == set([C1, C2])
    assert set(CC3.sub_species) == set([X, Y, Z])
    assert set(CC3.initial_states) == set([Complex([Z, Z]), X])
    assert CC3.intermediate_states is None

    #Test with intermediate states
    CC4 = CombinatorialComplex(final_states = [C1, C2], intermediate_states = [Complex([Z, Z])])
    assert set(CC4.final_states) == set([C1, C2])
    assert set(CC4.sub_species) == set([X, Y, Z])
    assert set(CC4.initial_states) == set([X, Y, Z])
    assert set(CC4.intermediate_states) == set([Complex([Z, Z])])

    #Test with excluded states
    CC5 = CombinatorialComplex(final_states = [C2], excluded_states = [C1])

    #Test cases that should produce errors

    #final_states must be ComplexSpecies
    with pytest.raises(ValueError):
        CC = CombinatorialComplex(final_states = [Species("C1")])

    #initial_states must be Species inside of final states
    with pytest.raises(ValueError):
        CC = CombinatorialComplex(final_states = [C1], initial_states = [Species("S")])

    #initial states which are Complexes can only contain things in the final_states
    with pytest.raises(ValueError):
        CC = CombinatorialComplex(final_states = [C1], initial_states = [Complex([Species("S"), Species("S2")])])

    #intermiade_states must be Complexes
    with pytest.raises(ValueError):
        CC = CombinatorialComplex(final_states = [C1], intermediate_states = [Species("S")])

    #intermediate_states must contain only things in final states
    with pytest.raises(ValueError):
        CC = CombinatorialComplex(final_states = [C1], intermediate_states = [Complex([Species("S"), Species("S2")])])

def test_CombinatorialComplex_compute_species_to_add():
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C = Complex([X, X, Y, Z])
    CC = CombinatorialComplex(final_states = [C])

    #A single X, Y, and Z should be added
    species_to_add = CC.compute_species_to_add(X, C)
    assert set(species_to_add) == set([X, Y, Z])
    assert species_to_add.count(X) == 1 and species_to_add.count(Y) == 1 and species_to_add.count(Z) == 1

    #Two X's, and Z should be added
    species_to_add = CC.compute_species_to_add(Y, C)
    assert set(species_to_add) == set([X, Z])
    assert species_to_add.count(X) == 2 and species_to_add.count(Z) == 1

    #Try with s0 as a ComplexSpecies
    species_to_add = CC.compute_species_to_add(Complex([X, X]), C)
    assert set(species_to_add) == set([Y, Z])
    assert species_to_add.count(Y) == 1 and species_to_add.count(Z) == 1

    #If sf is not a ComplexSpecies, there is an error
    with pytest.raises(ValueError):
        species_to_add = CC.compute_species_to_add(C, X)

    #In the following cases, sf cannot be created by adding species to s0, so None is returned

    #C contains more species than Complex([X, X])
    species_to_add = CC.compute_species_to_add(C, Complex([X, Y]))
    assert species_to_add is None
    #C contains more species than Complex([X, X])
    species_to_add = CC.compute_species_to_add(C, Complex([X, X]))
    assert species_to_add is None
    #S is not in C
    species_to_add = CC.compute_species_to_add(Species("S"), C)
    assert species_to_add is None
    species_to_add = CC.compute_species_to_add(Complex([Species("S"), X]), C)
    assert species_to_add is None

def test_CombinatorialComplex_get_combinations_between():
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C = Complex([X, Y, Z])
    CC = CombinatorialComplex(final_states = [C])

    combos = CC.get_combinations_between(X, C)

    #Below conditions uniquely specify all the correct returns
    assert len(combos) == 4
    #Assert first steps are in combos
    assert (Y, X, Complex([X, Y])) in combos and (Z, X, Complex([X, Z])) in combos
    #Second steps are in combos
    assert (Z, Complex([X, Y]), Complex([X, Y, Z])) in combos and (Y, Complex([X, Z]), Complex([X, Y, Z])) in combos

    #test again from a complex
    combos = CC.get_combinations_between(Complex([X, Y]), C)
    assert len(combos) == 1
    assert (Z, Complex([X, Y]), Complex([X, Y, Z])) in combos


    #Test with duplicate species
    C2 = Complex([X, X, Z])
    CC2 = CombinatorialComplex(final_states = [C2])
    combos = CC2.get_combinations_between(X, C2)
    assert len(combos) == 4
    assert (X, X, Complex([X, X])) in combos and (Z, X, Complex([X, Z])) in combos
    assert (Z, Complex([X, X]), C2) in combos and (X, Complex([X, Z]), C2) in combos

    #test again from a Complex
    combos = CC2.get_combinations_between(Complex([X, X]), C2)
    assert len(combos) == 1
    assert (Z, Complex([X, X]), C2) in combos
    combos = CC2.get_combinations_between(Complex([X, Z]), C2)
    assert len(combos) == 1
    assert (X, Complex([X, Z]), C2) in combos

    #Test with excluded states that are ComplexSpecies
    C3a = Complex([X, Y])
    C3b = Complex([X, Z])
    CC3 = CombinatorialComplex(final_states = [C], excluded_states = [C3a, C3b])
    combos = CC3.get_combinations_between(X, C)
    combos_list = [c[0] for c in combos]+[c[1] for c in combos] + [c[2] for c in combos]
    assert C3a not in combos and C3b not in combos

    #Test with excluded states that are Species
    CC4 = CombinatorialComplex(final_states = [C], excluded_states = [Y, Z])
    combos = CC4.get_combinations_between(X, C)
    combos_list = [c[0] for c in combos]+[c[1] for c in combos] + [c[2] for c in combos]
    assert Y not in combos and Z not in combos 

def test_CombinatorialComplex_update_species():
    mech_b = One_Step_Binding()

    #Basic Case
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C1 = Complex([X, Y, Z])
    CC1 = CombinatorialComplex(final_states = [C1], mechanisms = [mech_b])

    s1 = CC1.update_species()
    assert X in s1 and Y in s1 and Z in s1 and C1 in s1
    assert Complex([X, Y]) in s1 and Complex([X, Z]) in s1 and Complex([Y, Z]) in s1

    #Set initial states
    C2i1 = Complex([X, X, Y])
    C2i2 = Complex([Y, Z])
    C2 = Complex([X, X, Y, Z])
    CC2 = CombinatorialComplex(final_states = [C2], initial_states =[C2i1, C2i2], mechanisms = [mech_b])
    s2 = CC2.update_species()
    assert X in s2 and Y not in s2 and Z in s2
    assert C2i1 in s2 and C2i2 in s2 and C2 in s2
    assert Complex([X, Y, Z]) in s2
    assert Complex([X, X]) not in s2 and Complex([X, X, Z]) not in s2 
    assert Complex([X, X]) not in s2 and Complex([X, Y]) not in s2 and Complex([X, Z]) not in s2

    #set intermediate states
    CC3 = CombinatorialComplex(final_states = [C2], intermediate_states =[C2i1, C2i2], mechanisms = [mech_b])
    s3 = CC3.update_species()
    assert X in s3 and Y in s3 and Z in s3
    assert C2i1 in s3 and C2i2 in s3 and C2 in s3
    assert Complex([X, Y, Z]) in s3 and Complex([X, X]) in s3 and Complex([X, Y]) in s3
    assert Complex([X, X, Z]) not in s3 and Complex([X, Z]) not in s3

    #set initial and intermediate states
    C4 = Complex([X, Y, Y, Z, Z])
    CC4 = CombinatorialComplex(final_states = [C4], intermediate_states =[Complex([X, Y, Y]), Complex([X, Z, Z])], initial_states = [Complex([X, Y]), Complex([X, Z])], mechanisms = [mech_b])
    s4 = CC4.update_species()
    assert X not in s4 and Y in s4 and Z in s4
    assert Complex([X, Y, Y]) in s4 and Complex([X, Z, Z]) in s4 and Complex([X, Y]) in s4 and Complex([X, Z]) in s4
    assert Complex([Y, Y]) not in s4 and Complex([Z, Z]) not in s4 and Complex([X, Y, Z]) not in s4
    assert Complex([X, Y, Y, Z]) in s4 and Complex([X, Y, Z, Z]) in s4

    #multiple final states
    C5a = Complex([X, X, Y])
    C5b = Complex([X, X, Z])
    CC5 = CombinatorialComplex(final_states = [C5a, C5b], mechanisms = [mech_b])
    s5 = CC5.update_species()
    assert X in s5 and Y in s5 and Z in s5
    assert Complex([X, X]) in s5 and Complex([X, Z]) in s5 and Complex([X, Y]) in s5
    assert Complex([Y, Z]) not in s5


def test_CombinatorialComplex_update_reactions():
    mech_b = One_Step_Binding()
    ku, kb = 1, 1
    params = {"ku":ku, "kb":kb}

    #shortcut to write this faster
    def R(inputs, outputs):
        return Reaction.from_massaction(inputs, outputs, k_forward = kb, k_reverse = ku) 

    #Basic Case
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    CXYZ = Complex([X, Y, Z])
    CXY = Complex([X, Y])
    CXZ = Complex([X, Z])
    CYZ = Complex([Y, Z])

    CC1 = CombinatorialComplex(final_states = [CXYZ], mechanisms = [mech_b], parameters = params)

    r1 = CC1.update_reactions()
    r1_true = [R([X, Y], [CXY]), R([X, Z], [CXZ]),  R([Y, Z], [CYZ]), R([CXY, Z], [CXYZ]), R([CXZ, Y], [CXYZ]), R([CYZ, X], [CXYZ])]
    assert all([r in r1_true for r in r1]) and all([r in r1 for r in r1_true])

    #Set initial states
    CYZ = Complex([Y, Z])
    CXX = Complex([X, X])
    CXXY = Complex([X, X, Y])
    CXXZ = Complex([X, X, Z])
    CXXYZ = Complex([X, X, Y, Z])
    CC2 = CombinatorialComplex(final_states = [CXXYZ], initial_states =[CXX, CYZ], mechanisms = [mech_b], parameters = params)
    r2 = CC2.update_reactions()
    r2_true = [R([CXX, Y], [CXXY]), R([CXX, Z], [CXXZ]), R([CXXY, Z], [CXXYZ]), R([CXXZ, Y], [CXXYZ]), R([CYZ, X], [CXYZ]), R([CXYZ, X], [CXXYZ])]
    assert all([r in r2_true for r in r2]) and all([r in r2 for r in r2_true])
    
    #set intermediate states
    CC3 = CombinatorialComplex(final_states = [CXXYZ], intermediate_states =[CXX, CYZ], mechanisms = [mech_b], parameters = params)
    r3 = CC3.update_reactions()
    r3_true = [R([X, X], [CXX]), R([Y, Z], [CYZ])] + r2_true
    assert all([r in r3_true for r in r3]) and all([r in r3 for r in r3_true])

    #set initial and intermediate states
    CXXYYZ = Complex([X, X, Y, Y, Z])
    CXXYZZ = Complex([X, X, Y, Z, Z])
    CXXYYZZ = Complex([X, X, Y, Y, Z, Z])
    CC4 = CombinatorialComplex(final_states = [CXXYYZZ], intermediate_states =[CXXYZ], initial_states = [CXX, CYZ], mechanisms = [mech_b], parameters = params)
    r4 = CC4.update_reactions()
    r4_true =  [R([CXX, Y], [CXXY]), R([CXX, Z], [CXXZ]), R([CYZ, X], [CXYZ]), R([CXYZ, X], [CXXYZ]), R([CXXZ, Y], [CXXYZ]), R([CXXY, Z], [CXXYZ]), R([CXXYZ, Y], [CXXYYZ]), R([CXXYZ, Z], [CXXYZZ]), R([CXXYYZ, Z], [CXXYYZZ]), R([CXXYZZ, Y], [CXXYYZZ])]
    assert all([r in r4_true for r in r4]) and all([r in r4 for r in r4_true])

    #multiple final states
    CC5 = CombinatorialComplex(final_states = [CXXY, CXXZ], mechanisms = [mech_b], parameters = params)
    r5 = CC5.update_reactions()
    r5_true = [R([X, X], [CXX]), R([CXX, Y], [CXXY]), R([CXX, Z], [CXXZ]), R([X, Y], [CXY]), R([CXY, X], [CXXY]), R([X, Z], [CXZ]), R([CXZ, X], [CXXZ])]
    assert all([r in r5_true for r in r5]) and all([r in r5 for r in r5_true])







