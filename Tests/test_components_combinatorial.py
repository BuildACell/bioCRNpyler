#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import CombinatorialComplex, One_Step_Binding, Reaction, CombinatorialConformation
from biocrnpyler import Species, OrderedPolymerSpecies, Complex
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


def test_CombinatorialConformation_init():

    #Test getters and setters of properties via init

    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C = Complex([X, X])
    p1 = OrderedPolymerSpecies([X, Y, Z])
    pc1 = Complex([p1[0], p1[2]]).parent #this gets the PolymerConformation

    CC1 = CombinatorialConformation(initial_states = [p1], final_states = [pc1])
    
    assert len(CC1.initial_states) == 1 and len(CC1.final_states) == 1

    p1b = OrderedPolymerSpecies([X, Y, Complex([Z, Z])])
    CC1b = CombinatorialConformation(initial_states = [p1], final_states = [p1b])
    assert len(CC1.initial_states) == 1 and len(CC1.final_states) == 1

    p2 = OrderedPolymerSpecies([Z, Y, X])
    pc2 = Complex([pc1.polymers[0][1], p2[1]]).parent
    CC2 = CombinatorialConformation(initial_states = [pc1], final_states = [pc2])
    assert len(CC2.initial_states) == 1 and len(CC2.final_states) == 1

    CC3 = CombinatorialConformation(initial_states = [p1], final_states = [pc2], intermediate_states = [pc1])
    assert len(CC3.intermediate_states) == 1

    CC4 = CombinatorialConformation(initial_states = [p1], final_states = [pc2], excluded_states = [pc1])
    assert len(CC4.excluded_states) == 1

    CC5 = CombinatorialConformation(initial_states = [p1], final_states = [pc2], excluded_complexes = [C])
    assert len(CC5.excluded_complexes) == 1

    #The following should produce errors:
    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [C], final_states = [pc1])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [p1], final_states = [C])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [p1], final_states = [pc1], intermediate_states = [C])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [p1], final_states = [pc1], excluded_states = [C])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [p1], final_states = [pc1], excluded_complexes = [pc2])



def test_CombinatorialConformation_compute_complexes_to_add_to_polymer():
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    p0 = OrderedPolymerSpecies([X])
    p1 = OrderedPolymerSpecies([X, Y, Z])
    p2 = OrderedPolymerSpecies([Z, Y, X])
    p3 = OrderedPolymerSpecies([Complex([X, X]), Y, Complex([X, Z])])

    CC = CombinatorialConformation(initial_states = [], final_states = [])

    #An empty list is returned if the Polymers are the same ("Nothing to add")
    complexes_to_add = CC.compute_complexes_to_add_to_polymer(p0, p0)
    assert len(complexes_to_add) == 0

    #None is returned if the Polymers cannot be converted
    complexes_to_add = CC.compute_complexes_to_add_to_polymer(p0, p1)
    assert complexes_to_add == None
    complexes_to_add = CC.compute_complexes_to_add_to_polymer(p1, p2)
    assert complexes_to_add == None


    #Basic Case
    complexes_to_add = CC.compute_complexes_to_add_to_polymer(p1, p3)
    assert len(complexes_to_add) == 2
    assert (0, [X]) in complexes_to_add and (2, [X]) in complexes_to_add

    #Nested Complex Case
    p4 = OrderedPolymerSpecies([Complex([X, X]), Y, Complex([Z, Complex([X, X])])])
    complexes_to_add = CC.compute_complexes_to_add_to_polymer(p1, p4)
    assert len(complexes_to_add) == 2
    assert (0, [X]) in complexes_to_add and (2, [Complex([X, X])]) in complexes_to_add

    #These cases should give errors
    with pytest.raises(ValueError):
        CC.compute_complexes_to_add_to_polymer(X, p4)
    with pytest.raises(ValueError):
        CC.compute_complexes_to_add_to_polymer(p1, X)



def test_CombinatorialConformation_compute_polymer_mapping():
    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C = Complex([X, X])
    p1 = OrderedPolymerSpecies([X, Y, Z])
    p2 = OrderedPolymerSpecies([Z, Y, X])
    pc1 = Complex([p1[0], p1[2]]).parent #this gets the PolymerConformation

    #Create a PolymerConformation with two polymers
    pc2 = Complex([pc1.polymers[0][1], p2[1]]).parent


    CC = CombinatorialConformation(initial_states = [], final_states = [])

    #No Polymers to Add, but p1 is in pc1 so there is a mapping
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, pc1)
    assert len(polymers_to_add) == 0
    assert p1 in polymer_mapping and pc1.polymers[0] in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[p1] == pc1.polymers[0]
    assert polymer_mapping[pc1.polymers[0]] == p1
    assert polymer_mapping[p1, pc1.polymers[0]] == []

    #There is no path between these two polymers
    assert CC.compute_polymer_mapping(p1, p2) is None

    #There is no path from p2 to a conformation only containing p1
    assert CC.compute_polymer_mapping(p2, pc1) is None

    #One Polymer to add from a OrderedPolymerSpecies
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, pc2)
    assert len(polymers_to_add) == 1
    assert polymers_to_add[0] == p2
    assert p1 in polymer_mapping and pc2.polymers[0] in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[p1] == pc2.polymers[0]
    assert polymer_mapping[pc2.polymers[0]] == p1
    assert polymer_mapping[p1, pc2.polymers[0]] == []

    #One Polymer to add from a PolymerConformation
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(pc1, pc2)
    assert len(polymers_to_add) == 1
    assert str(polymers_to_add[0]) == str(p2)
    assert pc1.polymers[0] in polymer_mapping and pc2.polymers[0] in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[pc1.polymers[0]] == pc2.polymers[0]
    assert str(polymer_mapping[pc2.polymers[0]]) == str(p1)
    assert polymer_mapping[pc1.polymers[0], pc2.polymers[0]] == []

    #Create a polymers from p1 which has a new Complex in it
    p3 = Complex([p1[0], X, X]).parent
    p4 = Complex([p3[1], Y, Z]).parent

    #Polymer Mapping between Polymers which require additions
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, p3)
    assert len(polymers_to_add) == 0
    assert p1 in polymer_mapping and p3 in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[p1] == p3 and polymer_mapping[p3] == p1
    assert len(polymer_mapping[p1, p3]) == 1 and polymer_mapping[p1, p3][0] == (0, [X, X])

    #This case required additions in two places
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, p4)
    assert len(polymers_to_add) == 0
    assert p1 in polymer_mapping and p4 in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[p1] == p4 and polymer_mapping[p4] == p1
    assert len(polymer_mapping[p1, p4]) == 2 and polymer_mapping[p1, p4][0] == (0, [X, X]) and polymer_mapping[p1, p4][1] == (1, [Y, Z])
    

    #Try in a PolymerConformation with multiple Polymers
    pc3 = Complex([p4[1], p2[1]]).parent
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, pc3)
    assert len(polymers_to_add) == 1 and str(polymers_to_add[0]) == str(p2)
    assert p1 in polymer_mapping and pc3.polymers[1] in polymer_mapping and len(polymer_mapping) == 3
    assert polymer_mapping[p1] == pc3.polymers[1] and polymer_mapping[pc3.polymers[1]] == p1
    assert len(polymer_mapping[p1, pc3.polymers[1]]) == 2 and polymer_mapping[p1, pc3.polymers[1]][0] == (0, [X, X]) and polymer_mapping[p1, pc3.polymers[1]][1] == (1, [Y, Z])

    #Try between two polymer conformations each with multiple polymers
    polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(pc2, pc3)
    assert len(polymers_to_add) == 0
    assert pc2.polymers[0] in polymer_mapping and pc2.polymers[1] in polymer_mapping and pc3.polymers[0] in polymer_mapping and pc3.polymers[1] in polymer_mapping and len(polymer_mapping) == 6
    assert polymer_mapping[pc2.polymers[0]] == pc3.polymers[1] and polymer_mapping[pc3.polymers[1]] == pc2.polymers[0]
    assert polymer_mapping[pc2.polymers[1]] == pc3.polymers[0] and polymer_mapping[pc3.polymers[0]] == pc2.polymers[1]
    assert len(polymer_mapping[pc2.polymers[1], pc3.polymers[0]]) == 0
    assert len(polymer_mapping[pc2.polymers[0], pc3.polymers[1]]) == 2 and polymer_mapping[pc2.polymers[0], pc3.polymers[1]][0] == (0, [X, X]) and polymer_mapping[pc2.polymers[0], pc3.polymers[1]][1] == (1, [Y, Z])
    
    #ValueError Cases
    with pytest.raises(ValueError):
        CC.compute_polymer_mapping(X, pc3)

    with pytest.raises(ValueError):
        CC.compute_polymer_mapping(p1, X)

    #Not Implemented Errors (this happens when the Polymer Mapping is Ambigious - meaning there are multiple possible mappings)
    with pytest.raises(NotImplementedError):
        pc4 = Complex([p4[0], p3[0]]).parent
        polymer_mapping, polymers_to_add = CC.compute_polymer_mapping(p1, pc4)


def test_CombinatorialConformation_compute_complexes_to_add_to_conformation():
    CC = CombinatorialConformation(initial_states = [], final_states = [])

    X, Y, Z, S = Species("X"), Species("Y"), Species("Z"), Species("S")
    C = Complex([X, X])
    p1 = OrderedPolymerSpecies([X, Y, Z])
    p2 = OrderedPolymerSpecies([Z, Y, X])
    p3 = Complex([p1[0], X]).parent
    pc_c = Complex([p1[0], p1[2], S])
    pc1 = pc_c.parent #this gets the PolymerConformation
    pc2 = Complex([p1[0], p2[2], S]).parent #A polymer Conformation with two Polymers
    pc3 = Complex([pc1.polymers[0][1], p2[1], S, S]).parent #A polymer Conformation with two polymers and two Complexes
    pc_c2 = Complex([pc_c, pc_c.parent.polymers[0][1], S]) #add something to a Complex inside a PolymerConformation
    pc4 = pc_c2.parent
    pc5 = Complex([p3[1], p3[2]]).parent #A conformation which has a complex in a p1 and a complex between two monomers in p1

    #Starting with an OrderedPolymer Species and Creating an OrderedPolymerSpecies should result in None
    assert CC.compute_complexes_to_add_to_conformation(p1, p2) is None

    #Starting with a OrderedPolymerSpecies and creating a PolymerConformation
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(p1, pc1)
    assert complexes_to_add[0][0][0] == (p1, 0) and complexes_to_add[0][0][1] == (p1, 2) and len(complexes_to_add[0][1]) == 1 and complexes_to_add[0][1][0] == S

    #pc1 cannot be created from p2
    assert CC.compute_complexes_to_add_to_conformation(p2, pc1) is None


    #pc2 can be created from p1
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(p1, pc2)
    assert complexes_to_add[0][0][0] == (p1, 0) and complexes_to_add[0][0][1] == (pc2.polymers[1], 2) and len(complexes_to_add[0][1]) == 1 and complexes_to_add[0][1][0] == S

    #pc2 cannot be created from pc1
    assert CC.compute_complexes_to_add_to_conformation(pc1, pc2) is None

    #pc3 can be created from p1
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(p1, pc3)
    assert len(complexes_to_add) == 2
    assert complexes_to_add[0][0][0] == (p1, 0) and complexes_to_add[0][0][1] == (p1, 2) and len(complexes_to_add[0][1]) == 1 and complexes_to_add[0][1][0] == S
    assert complexes_to_add[1][0][0] == (p1, 1) and complexes_to_add[1][0][1] == (pc3.polymers[1], 1) and len(complexes_to_add[1][1]) == 2 and complexes_to_add[1][1][0] == S
    

    #pc3 can be created from pc1
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(pc1, pc3)
    assert complexes_to_add[0][0][0] == (pc1.polymers[0], 1) and complexes_to_add[0][0][1] == (pc3.polymers[1], 1) and len(complexes_to_add[0][1]) == 2 and complexes_to_add[0][1][0] == S
    
    
    

    #pc4 can be created from pc1 by adding to an existing complex
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(pc1, pc4)

    assert complexes_to_add[0][0][0][1].monomer_eq(pc_c)
    assert str(complexes_to_add[0][0][1]) == str((pc4.polymers[0], 1)) #str used to ignore parents/location
    assert len(complexes_to_add[0][1]) == 1 and complexes_to_add[0][1][0] == S

    #In this case, a Complex has to be added to p1 and a Complex has to be added between two monomers in p1
    cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(p1, pc5)
    assert len(complexes_to_add) == 1 and complexes_to_add[0] == ([(p1, 1), (p1, 2)], []) 

    #ValueError Cases
    with pytest.raises(ValueError):
        CC.compute_complexes_to_add_to_conformation(X, pc3)

    with pytest.raises(ValueError):
        CC.compute_complexes_to_add_to_conformation(p1, X)

    #Not Implemented Errors (this happens when the Polymer Mapping is Ambigious - meaning there are multiple possible mappings)
    with pytest.raises(NotImplementedError):
        pc5 = Complex([pc1.polymers[0][1], pc4.polymers[0][1]]).parent
        cm, complexes_to_add = CC.compute_complexes_to_add_to_conformation(p1, pc5)


def test_CombinatorialConformation_get_combinations_between():
    CC = CombinatorialConformation(initial_states = [], final_states = [])

    X, Y, Z, S = Species("X"), Species("Y"), Species("Z"), Species("S")
    p1 = OrderedPolymerSpecies([X, Y, Z])
    p2 = OrderedPolymerSpecies([Z, Y, X])
    p3 = Complex([p1[0], X]).parent #p1 with a Complex
    p4 = Complex([p3[1], Y]).parent #p1 with two Complexes
    pc_c = Complex([p1[0], p1[2], S])
    pc1 = pc_c.parent #this gets the PolymerConformation
    pc2 = Complex([p1[0], p2[2], S]).parent #A polymer Conformation with two Polymers
    pc3 = Complex([pc1.polymers[0][1], p2[1], S, S]).parent #A polymer Conformation with two polymers and two Complexes
    pc4 = Complex([p3[1], p3[2]]).parent #p1 with a Complex in the polymer and a Complex between two monomers
    pc5 = Complex([pc_c, p2[0], X]).parent #Form a Complex around a Complex already in an PolymerConformation
    pc6 = Complex([pc3.complexes[0], pc3.complexes[1]]).parent #Bind two complexes together from the same conformation


    #In this case, a species are added in one place to p1 to create a new polymer p3
    combos = CC.get_combinations_between(p1, p3)
    assert len(combos) == 1 and len(combos[0]) == 1 and combos[0][0] == ([(p1, 0)], [X])

    #In this case, species are added in two places to p1 to get a new polymer p4
    combos = CC.get_combinations_between(p1, p4)
    assert len(combos) == 2 and len(combos[0]) == 2 and len(combos[1]) == 2
    assert combos[0][0] == ([(p1, 0)], [X]) and combos[0][1] == ([(p1, 1)], [Y])
    assert combos[1][1] == ([(p1, 0)], [X]) and combos[1][0] == ([(p1, 1)], [Y])

    #In this case, p1 binds to itself to create pc1
    combos = CC.get_combinations_between(p1, pc1)
    assert len(combos) == 1 and len(combos[0]) == 1
    assert combos[0][0] == ([(p1, 0), (p1, 2)], [S])

    #In this case, an additional polymer p2 must be added to p1 as well
    combos = CC.get_combinations_between(p1, pc2)
    assert len(combos) == 1 and len(combos[0]) == 1
    assert combos[0][0] == ([(p1, 0), (pc2.polymers[1], 2)], [S])

    #In this case, there are two complexes between two polymers
    combos = CC.get_combinations_between(p1, pc3)
    assert len(combos) == 2 and len(combos[0]) == 2 and len(combos[1]) == 2
    assert combos[0][0] == ([(p1, 0), (p1, 2)], [S]) and combos[0][1] == ([(p1, 1), (pc3.polymers[1], 1)], [S, S])
    assert combos[1][1] == ([(p1, 0), (p1, 2)], [S]) and combos[1][0] == ([(p1, 1), (pc3.polymers[1], 1)], [S, S])

    #In this case, a PolymerConformation becomes a different PolymerConformation in a single step
    combos = CC.get_combinations_between(pc1, pc3)
    assert len(combos) == 1 and len(combos[0]) == 1 and combos[0][0] == ([(pc1.polymers[0], 1), (pc3.polymers[1], 1)], [S, S])

    #In this case, p1 gains a complex at index 0 and index 1 and 2 bind together forming a PolymerConformation in two steps
    combos = CC.get_combinations_between(p1, pc4)
    assert len(combos) == 2 and len(combos[0]) == 2 and len(combos[1]) == 2
    assert combos[0][0] == ([(p1, 1), (p1, 2)], []) and combos[0][1] == ([(p1, 0)], [X])
    assert combos[1][1] == ([(p1, 1), (p1, 2)], []) and combos[1][0] == ([(p1, 0)], [X])
    

    #In this case, a Polymer is bound to a complex already in a polymer conformation
    combos = CC.get_combinations_between(pc1, pc5)
    assert len(combos) == 1 and len(combos[0]) == 1 and combos[0][0] == ([(pc1, pc_c), (pc5.polymers[1], 0)], [X])

    #Bind two complexes together in conformations
    combos = CC.get_combinations_between(pc3, pc6)
    assert len(combos) == 1 and len(combos[0]) == 1 and combos[0][0] == ([(pc3, pc3.complexes[1]), (pc3, pc3.complexes[0])], [])
    
    assert CC.get_combinations_between(p1, p1) is None
    assert CC.get_combinations_between(p1, p2) is None
    assert CC.get_combinations_between(p2, pc1) is None
    assert CC.get_combinations_between(pc1, Complex([p2[0], p2[1]]).parent) is None





