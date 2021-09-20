#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import CombinatorialComplex, One_Step_Binding, Reaction, CombinatorialConformation
from biocrnpyler import Species, OrderedPolymerSpecies, Complex, PolymerConformation
import pytest


def test_CombinatorialConformation_init():

    #Test getters and setters of properties via init

    X, Y, Z = Species("X"), Species("Y"), Species("Z")
    C = Complex([X, X])
    p1 = OrderedPolymerSpecies([X, Y, Z])
    pc0 = PolymerConformation(polymer = p1)
    pc1 = Complex([pc0.polymers[0][0], pc0.polymers[0][2]]).parent #this gets the PolymerConformation

    CC1 = CombinatorialConformation(initial_states = [pc0], final_states = [pc1])
    
    assert len(CC1.initial_states) == 1 and len(CC1.final_states) == 1

    pc1b = Complex([pc0.polymers[0][2], Z]).parent
    CC1b = CombinatorialConformation(initial_states = [pc0], final_states = [pc1b])
    assert len(CC1.initial_states) == 1 and len(CC1.final_states) == 1

    pc1c = Complex([pc1b.polymers[0][0], pc1b.polymers[0][1]]).parent
    CC1c0 = CombinatorialConformation(initial_states = None, final_states = [pc1b, pc1c])
    assert len(CC1c0.final_states) == 2
    assert len(CC1c0.initial_states) == 1
    CC1c1 = CombinatorialConformation(initial_states = [pc0, pc1b], final_states = [pc1c])
    assert len(CC1c1.initial_states) == 2


    #The following should produce errors:
    p2 = OrderedPolymerSpecies([Z, Y, X])
    pc2 = Complex([pc1.polymers[0][1], PolymerConformation(polymer = p2).polymers[0][1]]).parent
    with pytest.raises(ValueError): #Multiple internal polymers
        CC = CombinatorialConformation(initial_states = [pc1], final_states = [pc2])

    with pytest.raises(ValueError):#Multiple internal polymers
        CC = CombinatorialConformation(initial_states = [pc0, pc1], final_states = [pc2])

    with pytest.raises(ValueError):#Multiple internal polymers
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [pc1], intermediate_states = [pc2])
    with pytest.raises(ValueError):#Multiple internal polymers
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [pc2], excluded_states = [pc1])
    with pytest.raises(ValueError):#Multiple internal polymers
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [pc1], excluded_states = [pc2])

    #Non PolymerConformations passed as initial/final/intermediate states
    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [C], final_states = [pc1])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [C])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [pc1], intermediate_states = [C])

    with pytest.raises(ValueError):
        CC = CombinatorialConformation(initial_states = [pc0], final_states = [pc1], excluded_states = [C])

    




def pass_test_CombinatorialConformation_compute_complexes_to_add_to_polymer():
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

def test_compute_species_changes():
    X, Y, Z, S = Species("X"), Species("Y"), Species("Z"), Species("S")
    C = Complex([X, X])
    p0 = OrderedPolymerSpecies([X, Y, Z, C])
    pc0 = PolymerConformation(polymer = p0)
    c1 = Complex([pc0.polymers[0][0], pc0.polymers[0][1]])
    pc1 = c1.parent
    ind_c1 = pc1.get_polymer_positions(c1, 0)
    c2 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], S, S])
    pc2 = c2.parent
    ind_c2 = pc2.get_polymer_positions(c2, 0)
    c3 = Complex([pc1.polymers[0][2], pc1.polymers[0][3], Z])
    pc3 = c3.parent
    ind_c3 = pc3.get_polymer_positions(c3, 0)
    c4 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], pc0.polymers[0][3], Z])
    pc4 = c4.parent
    ind_c4 = pc4.get_polymer_positions(c4, 0)


    print("pc0", pc0, "\npc1", pc1, "\npc2", pc2, "\npc3", pc3, "\npc4", pc4)


    CC = CombinatorialConformation(initial_states = [pc0], intermediate_states = [pc1], final_states = [pc2, pc3])

    #No additional species added
    SC, MC = CC.compute_species_changes(s0 = pc0, sf = pc1)
    assert (c1, ind_c1) in SC
    assert len(SC[(c1, ind_c1)]) == 0
    assert (c1, ind_c1) in MC
    assert len(MC[(c1, ind_c1)]) == 0

    #External Species added
    SC, MC = CC.compute_species_changes(s0 = pc0, sf = pc2)
    assert (c2, ind_c2) in SC
    assert len(SC[(c2, ind_c2)]) == 2 and S in SC[(c2, ind_c2)]
    assert all([len(MC[cf]) == 0 for cf in MC])

    #Multiple Complexes created
    SC, MC = CC.compute_species_changes(s0 = pc0, sf = pc3)
    assert (pc3.complexes[0], ind_c1) in SC and (pc3.complexes[1], ind_c3) in SC
    assert len(SC[pc3.complexes[0], ind_c1]) == 0
    assert len(SC[pc3.complexes[1], ind_c3]) == 1 and Z in SC[pc3.complexes[1], ind_c3]
    assert all([len(MC[cf]) == 0 for cf in MC])

    #A single complex created, from a conformation with complexes
    SC, MC = CC.compute_species_changes(s0 = pc1, sf = pc3)
    assert (pc3.complexes[0], ind_c1) not in SC and (pc3.complexes[1], ind_c3) in SC
    assert len(SC[pc3.complexes[1], ind_c3]) == 1 and Z in SC[pc3.complexes[1], ind_c3]
    assert len(MC[pc3.complexes[1], ind_c3]) == 0

    #Adding species to an existing Complex
    SC, MC = CC.compute_species_changes(s0 = pc1, sf = pc2)
    assert (pc1.complexes[0], ind_c1) not in SC and (pc2.complexes[0], ind_c2) in SC
    assert len(SC[pc2.complexes[0], ind_c2]) == 2 and S in SC[pc2.complexes[0], ind_c2]
    assert pc1.complexes[0] in MC[pc2.complexes[0], ind_c2]

    #adding species (including inside the polymer) to an existing complex
    SC, MC = CC.compute_species_changes(s0 = pc1, sf = pc4)
    assert (pc1.complexes[0], ind_c1) not in SC and (pc4.complexes[0], ind_c4) in SC
    assert Z in SC[pc4.complexes[0], ind_c4] and len(SC[pc4.complexes[0], ind_c4]) == 1
    assert pc1.complexes[0] in MC[pc4.complexes[0], ind_c4]

    #merging two complexes
    SC, MC = CC.compute_species_changes(s0 = pc3, sf = pc4)
    assert (pc3.complexes[0], ind_c1) not in SC and (pc3.complexes[1], ind_c3) not in SC and (pc4.complexes[0], ind_c4) not in SC
    assert pc3.complexes[0] in MC[pc4.complexes[0], ind_c4] and pc3.complexes[1] in MC[pc4.complexes[0], ind_c4]

    #Cases where compute_species_changes should return False

    #No changes between the conformations
    assert not CC.compute_species_changes(s0 = pc0, sf = pc0)

    #because s0 contains more species/complexes than sf
    assert not CC.compute_species_changes(s0 = pc3, sf = pc0)
    assert not CC.compute_species_changes(s0 = pc3, sf = pc1)
    assert not CC.compute_species_changes(s0 = pc2, sf = pc1)
    assert not CC.compute_species_changes(s0 = pc2, sf = pc3)
    assert not CC.compute_species_changes(s0 = pc3, sf = pc2)


def test_get_combinations_between():
    X, Y, Z, S = Species("X"), Species("Y"), Species("Z"), Species("S")
    C = Complex([X, X])
    p0 = OrderedPolymerSpecies([X, Y, Z, C])
    pc0 = PolymerConformation(polymer = p0)
    c1 = Complex([pc0.polymers[0][0], pc0.polymers[0][1]])
    pc1 = c1.parent
    c2 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], S, S])
    pc2 = c2.parent
    c3 = Complex([pc1.polymers[0][2], pc1.polymers[0][3], Z])
    pc3 = c3.parent
    c4 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], pc0.polymers[0][3], Z])
    pc4 = c4.parent
    c5a = Complex([pc0.polymers[0][0], pc0.polymers[0][2]])
    pc5a = c5a.parent
    c5b = Complex([pc5a.polymers[0][1], pc5a.polymers[0][3]])
    pc5b = c5b.parent


    print("pc0", pc0, "\npc1", pc1, "\npc2", pc2, "\npc3", pc3, "\npc4", pc4)


    #Arguments here don't matter in this test
    CC = CombinatorialConformation(initial_states = [], intermediate_states = [], final_states = [pc3])

    #No additional species added
    perms = CC.get_combinations_between(s0 = pc0, sf = pc1)
    assert len(perms) == 1

    #External Species added
    perms = CC.get_combinations_between(s0 = pc0, sf = pc2)
    assert len(perms) == 1

    #Multiple Complexes created
    perms = CC.get_combinations_between(s0 = pc0, sf = pc3)
    assert len(perms) == 4

    #A single complex created, from a conformation with complexes
    perms = CC.get_combinations_between(s0 = pc1, sf = pc3)
    assert len(perms) == 1

    #Adding species to an existing Complex
    perms = CC.get_combinations_between(s0 = pc1, sf = pc2)
    assert len(perms) == 1

    #adding species (including inside the polymer) to an existing complex
    perms = CC.get_combinations_between(s0 = pc1, sf = pc4)
    assert len(perms) == 1

    #merging two complexes
    perms = CC.get_combinations_between(s0 = pc3, sf = pc4)
    assert len(perms) == 1
   
    #Adding an excluded_state
    #excluded_state matters here
    CC = CombinatorialConformation(initial_states = [], excluded_states = [pc1], final_states = [pc3])
    perms = CC.get_combinations_between(s0 = pc0, sf = pc3)
    assert len(perms) == 2

    #Cases where compute_species_changes should return False
    #s0 contains more species/complexes than sf
    assert len(CC.get_combinations_between(s0 = pc3, sf = pc0)) == 0
    assert len(CC.get_combinations_between(s0 = pc3, sf = pc1)) == 0
    assert len(CC.get_combinations_between(s0 = pc2, sf = pc1)) == 0
    assert len(CC.get_combinations_between(s0 = pc2, sf = pc3)) == 0
    assert len(CC.get_combinations_between(s0 = pc3, sf = pc2)) == 0
    #Partially overlapping sets of monomers are bound
    assert len(CC.get_combinations_between(s0 = pc3, sf = pc5a))== 0
    assert len(CC.get_combinations_between(s0 = pc5a, sf = pc3))== 0
    #These ones are especially tricky because the same Monomers are bound, but in different sets of complexes
    assert len(CC.get_combinations_between(s0 = pc3, sf = pc5b))== 0
    assert len(CC.get_combinations_between(s0 = pc5b, sf = pc3))== 0


def test_update_species_and_reactions():
    params = {"kf":1, "kr":1}

    X, Y, Z, S = Species("X"), Species("Y"), Species("Z"), Species("S")
    C = Complex([X, X])
    p0 = OrderedPolymerSpecies([X, Y, Z, C])
    pc0 = PolymerConformation(polymer = p0)
    c1 = Complex([pc0.polymers[0][0], pc0.polymers[0][1]])
    pc1 = c1.parent
    c2 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], S, S])
    pc2 = c2.parent
    c3 = Complex([pc1.polymers[0][2], pc1.polymers[0][3], Z])
    pc3 = c3.parent
    c4 = Complex([pc0.polymers[0][0], pc0.polymers[0][1], pc0.polymers[0][2], pc0.polymers[0][3], Z])
    pc4 = c4.parent

    #No conformation changes are possible
    CC0 = CombinatorialConformation(initial_states = [], final_states = [pc0], parameters = params)
    species = CC0.update_species()
    reactions = CC0.update_reactions()
    assert len(species) == 0
    assert len(reactions) == 0

    #Only a single binding reaction can occur
    CC1 = CombinatorialConformation(initial_states = [], final_states = [pc1], parameters = params)
    species = CC1.update_species()
    reactions = CC1.update_reactions()
    assert len(species) == 2
    assert pc1 in species and pc0 in species
    assert len(reactions) == 1

    #This should be the same as above, by default
    CC1b = CombinatorialConformation(initial_states = [pc0], final_states = [pc1], parameters = params)
    species = CC1b.update_species()
    reactions = CC1b.update_reactions()
    assert len(species) == 2
    assert pc1 in species and pc0 in species
    assert len(reactions) == 1

    #Test multiple final states
    CC2 = CombinatorialConformation(initial_states = [pc0], final_states = [pc1, pc2], parameters = params)
    species = CC2.update_species()
    reactions = CC2.update_reactions()
    assert len(species) == 4
    assert pc1 in species and pc0 in species and pc2 in species and S in species
    assert len(reactions) == 2

    #Test multiple pathways
    CC3 = CombinatorialConformation(initial_states = [pc0], final_states = [pc3], parameters = params)
    species = CC3.update_species()
    reactions = CC3.update_reactions()
    assert len(species) == 5
    assert pc1 in species and pc0 in species and pc3 in species and Z in species
    assert len(reactions) == 4

    #Test intermediate states
    CC3b = CombinatorialConformation(initial_states = [pc0], intermediate_states = [pc1], final_states = [pc3], parameters = params)
    species = CC3b.update_species()
    reactions = CC3b.update_reactions()
    assert len(species) == 4
    assert pc1 in species and pc0 in species and pc1 in species and Z in species and pc3 in species
    assert len(reactions) == 2

    #Test excluded intermediate states
    CC3c = CombinatorialConformation(initial_states = [pc0], excluded_states = [pc1], final_states = [pc3], parameters = params)
    species = CC3c.update_species()
    reactions = CC3c.update_reactions()
    assert len(species) == 4
    assert pc0 in species and pc1 not in species and Z in species and pc3 in species
    assert len(reactions) == 2

    #Test adding a dead end intermediate species
    CC3c = CombinatorialConformation(initial_states = [pc0], intermediate_states = [pc1, pc2], final_states = [pc3], parameters = params)
    species = CC3c.update_species()
    reactions = CC3c.update_reactions()
    assert len(species) == 6
    assert pc0 in species and pc1 in species and Z in species and pc3 in species and pc2 in species and S in species
    assert len(reactions) == 3

    #Test expanding a pathway through intermediates
    CC4a = CombinatorialConformation(initial_states = [pc0], final_states = [pc4], parameters = params)
    species = CC4a.update_species()
    reactions = CC4a.update_reactions()
    assert len(species) == 3
    assert pc0 in species and pc4 in species
    assert len(reactions) == 1

    CC4b = CombinatorialConformation(initial_states = [pc0], intermediate_states = [pc3], final_states = [pc4], parameters = params)
    species = CC4b.update_species()
    reactions = CC4b.update_reactions()
    assert len(species) == 6
    assert pc0 in species and pc4 in species and pc1 in species and pc3 in species
    assert len(reactions) == 5


