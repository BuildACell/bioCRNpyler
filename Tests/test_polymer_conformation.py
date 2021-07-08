#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import *

from biocrnpyler import OrderedPolymerSpecies, PolymerConformation, Complex, Species

def test_polymer_conformation_no_complex_instantiation():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    pc = PolymerConformation(polymer = p)
    c1 = ComplexSpecies([p[0], p[1], Species("S")], called_from_complex = True)

    assert str(pc) == str(p) #these should have the same name!
    assert str(p) == str(pc.polymers[0])
    assert len(pc.complexes) == 0

    with pytest.raises(AssertionError):
        pc = PolymerConformation(polymer = [p, p])

    with pytest.raises(ValueError):
        pc = PolymerConformation(polymer = p, complexes = [c1])

    with pytest.raises(ValueError):
        pc = PolymerConformation()

def test_single_polymer_single_complex_instantiation():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p[0], p[1], Species("S")], called_from_complex = True)
    pc = PolymerConformation(complexes = [c1])
    #Test naming convention
    assert str(pc) == f"conformation__{p}_np0p0_{c1}_"

    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    c2 = ComplexSpecies([p[0], p2[0]], called_from_complex = True)

    #In these cases, TypeErrors should be raised
    with pytest.raises(ValueError):
        #Emtpy complexes is not allowed
        pc2 = PolymerConformation(complexes = [])

    with pytest.raises(ValueError):
        #set a PolymerConformation as a parent of the Species.
        S = Species("S")
        S.parent = pc
        c3 = ComplexSpecies([p[0], S], called_from_complex = True)
        pc2 = PolymerConformation(complexes = [c3])

    with pytest.raises(ValueError):
        #Cannot place an entire polymer into a Complex
        S = Species("S")
        S.parent = pc
        c3 = ComplexSpecies([p, S], called_from_complex = True)
        pc2 = PolymerConformation(complexes = [c3])


def test_multiple_polymer_multiple_complex_instantiation():
    p1 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    
    c1 = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    c2 = ComplexSpecies([p1[0], p2[2]], called_from_complex = True)
    c3 = ComplexSpecies([Species("S1"), Species("S2")], called_from_complex = True)
    
    
    pc1 = PolymerConformation(complexes = [c1, c2])
    
    #Check naming convention
    assert str(pc1) == f"conformation__{p1}_{p2}_p0p0_{c1}_p0p1_{c2}_"

    #check that order doesn't matter
    pc1_r = PolymerConformation(complexes = [c2, c1])
    assert pc1 == pc1_r

    #check that ComplexSpecies with the same string representation but different polymers are treated differently
    c1b = ComplexSpecies([p2[1], p1[0]], called_from_complex = True) #c1 and c1b have the same string representation - but connect different polymers
    pc1b = PolymerConformation(complexes = [c1b, c2])
    assert str(pc1b) == f"conformation__{p1}_{p2}_p0p1_{c1}_p0p1_{c2}_"
    assert pc1 != pc1b

    #Check that order doesn't matter for Complexes with the same string representation
    pc1br = PolymerConformation(complexes = [c2, c1b]) #reverse case
    assert pc1b == pc1br

    pc2 = PolymerConformation(complexes = [c1, c2])
    pc2b = PolymerConformation(complexes = [c1b, c2])
    assert pc2 != pc2b

    #In these cases, value errors should be raised
    with pytest.raises(ValueError):
        #c3 is not part of the conformation
        pc = PolymerConformation(complexes = [c1, c3])

    #try 2 polymers with complexes that look the same
    p3 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c4 = ComplexSpecies([p1[0], p3[0]])

    pc3 = PolymerConformation([c1, c4])
    assert str(pc3) == f"conformation__{p1}_{p1}_p0p0_{c1}_p0p1_{c4}_"

    #check that order doesn't matter for polymers which look the same
    pc3r = PolymerConformation([c4, c1])
    assert pc3 == pc3r

    with pytest.raises(ValueError):
        #duplicate Complexes are not allowed (same object case)
        pc4 = PolymerConformation([c1, c1])

    c1b = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    with pytest.raises(ValueError):
        #duplicate Complexes are not allowed (identical object case)
        pc4 = PolymerConformation([c1, c1b])


def test_from_polymer_conformation():
    #tests the classmethod .from_polymer_conformation

    #Produce a PolymerConformation
    p1 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    pc1 = PolymerConformation(complexes = [c1])

    #This Complex uses the polymer inside the pc1 so there is only a single polymer in the final PolymerConformation
    c1b = ComplexSpecies([pc1.polymers[0][0], pc1.polymers[0][2]], called_from_complex = True)
    c1bb = ComplexSpecies([p1[0], p1[2]])
    pc1b = PolymerConformation.from_polymer_conformation([pc1], [c1b])

    #This is correct (and checking ordering effects don't matter)
    assert pc1b == PolymerConformation([c1, c1bb]) == PolymerConformation([c1bb, c1]) 

    #This is incorrect because of polymers being copied into conformations
    assert pc1b != PolymerConformation([c1, c1b])

    #This Complex uses the old Polymer Instance which will represent the binding to a new, identical, PolymerSpecies
    #all three of the following therefore should be the same
    c1c = ComplexSpecies([p1[0], p1[2]], called_from_complex = True)
    pc1c = PolymerConformation.from_polymer_conformation([pc1], [c1c])

    p1d = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1d = ComplexSpecies([p1d[0], p1d[2]], called_from_complex = True)
    pc1d = PolymerConformation.from_polymer_conformation([pc1], [c1c])

    #Create directly
    pc1e = PolymerConformation([c1d, c1])
    pc1er = PolymerConformation([c1, c1d]) #Order shouldn't matter
    
    assert pc1c == pc1d == pc1e == pc1er

    #Test from multiple PolymerConformations
    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    c2 = ComplexSpecies([p2[0], p2[1]], called_from_complex = True)
    pc2 = PolymerConformation(complexes = [c2])
    #Create a PolymerConformation that links the two Polymers together
    c3 = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    pc3 = PolymerConformation([c3])

    #Add an additional binding site connecting pc2 and pc3
    c4 = ComplexSpecies([pc3.polymers[0][0], pc2.polymers[0][1]]) #get the conformations polymers due to copying
    pc4 = PolymerConformation.from_polymer_conformation([pc2, pc3], [c4])
    c4b = ComplexSpecies([p1[0], p2[1]])

    assert pc4 == PolymerConformation([c2, c3, c4b]) == PolymerConformation([c4b, c2, c3]) #And order doesn't matter
    assert pc4 != PolymerConformation([c2, c3, c4]) #This should be different due to copying

    #The following should produce errors
    with pytest.raises(TypeError):
        pc1b = PolymerConformation.from_polymer_conformation(pc1, [c1b])
    with pytest.raises(TypeError):
        pc1b = PolymerConformation.from_polymer_conformation([c1b], [c1b])

def test_from_polymer_replacement():
    #tests the classmethod .from_polymer_replacement

    #Produce a PolymerConformation
    p1 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p1[0], p1[1], Species("S"), Species("S")], called_from_complex = True)
    pc1 = PolymerConformation(complexes = [c1])

    #Produce a second Polymer
    p2 = OrderedPolymerSpecies([Species("m3"), Species("m2"), Species("m1")])
    #replace p1 with p2
    pc1_replaced = PolymerConformation.from_polymer_replacement(pc1, [pc1.polymers[0]], [p2])

    #This should be equivalent to creating a new PolymerConformation
    c2 = ComplexSpecies([p2[0], p2[1], Species("S"), Species("S")], called_from_complex = True)
    pc2 = PolymerConformation(complexes = [c2])
    assert pc1_replaced == pc2

    #Replacing a Polymer with the same polymer is valid, but doesn't do anything
    pc1_replaced_b = PolymerConformation.from_polymer_replacement(pc1, [pc1.polymers[0]], [p1])
    assert pc1_replaced_b == pc1

    
    #These conditions should raise ValueErrors
    with pytest.raises(TypeError):
        PolymerConformation.from_polymer_replacement(pc1, [pc1.polymers[0]], None)

    with pytest.raises(ValueError):
        PolymerConformation.from_polymer_replacement(pc1, [None], [p2])

    with pytest.raises(TypeError):
        PolymerConformation.from_polymer_replacement(pc1, [pc1.polymers[0]], [None])

    with pytest.raises(ValueError):
        PolymerConformation.from_polymer_replacement(pc1, [p1], [p2])

    with pytest.raises(ValueError):
        PolymerConformation.from_polymer_replacement(pc1, [pc1.polymers[0]], [p1, p2])

    with pytest.raises(TypeError):
        PolymerConformation.from_polymer_replacement(None, [pc1.polymers[0]], [p1, p2])


def test_get_complex():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p[0], p[1], Species("S")], called_from_complex = True)
    pc = PolymerConformation(complexes = [c1])

    assert str(pc.get_complex(c1)) == str(c1) #these are not technically equal because they have different parents

    
def test_polymer_conformation_ordered_polymer_species_name_equality():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    pc = PolymerConformation(polymer = p)
    assert str(p) == str(pc)
    


