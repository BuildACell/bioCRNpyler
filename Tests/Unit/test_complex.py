#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import OrderedPolymer, OrderedMonomer, OrderedPolymerSpecies, Species, ComplexSpecies, Complex, OrderedComplexSpecies, PolymerConformation

def test_complex_no_polymer():
    a = Species('A')
    b = Species('B')
    #Create a ComplexSpecies
    c = Complex([a,b])
    assert(c == ComplexSpecies([a, b], called_from_complex = True))

    #Ordered Complex
    truth = OrderedComplexSpecies([a,Complex([b,a]),a])
    testcomplx = Complex([a,Complex([b,a]),a],ordered=True)
    assert(truth == testcomplx)

def test_complex_with_polymer():
    a = Species('A')
    b = Species('B')
    c = Complex([a,b])
    d = OrderedPolymerSpecies([a,b,a])

    #Complex in an OrderedPolymerSpecies
    d_c = Complex([d[1],a])

    assert(isinstance(d_c, ComplexSpecies))
    truth = OrderedPolymerSpecies([a,Complex([b,a]),a])
    assert(d_c.parent == truth)

    assert(d[1] not in d_c) #different parents
    assert(b in d_c) #b is unbound when put in the Complex

    #Ordered case
    d_co = Complex([d[1],a], ordered = True)
    assert isinstance(d_co, OrderedComplexSpecies)
    truth_o = OrderedPolymerSpecies([a, Complex([b,a], ordered = True),a])
    assert d_co.parent==truth_o

    #Cannot complex two Species inside a PolymerSpecies without puting the PolymerSpecies into a Conformation first
    with pytest.raises(TypeError):
        c = Complex([Species("S"), d[0], d[2]])


def test_complex_with_single_polymer():
    a = Species('A')
    b = Species('B')
    p = OrderedPolymerSpecies([a,b,a])

    #This should just produce a ComplexSpecies around the PolymerSpecies
    c2 = Complex([p, Species("S")])
    assert c2 == ComplexSpecies([p, Species("S")], called_from_complex = True)

    #If the Polymer is in a Conformation, it cannot be Complexed.
    pc = PolymerConformation(polymer = p)
    assert str(pc.polymers[0]) == str(p)
    assert pc.polymers[0].parent == pc
    assert p.parent is None

    with pytest.raises(ValueError):
        c2 = Complex([pc.polymers[0], Species("S")])

    #A monomer form the polymer can still be complexed, however
    c2 = Complex([pc.polymers[0][0], Species("S")])
    assert isinstance(c2.parent, PolymerConformation)
    assert c2.parent != pc
    assert len(c2.parent.complexes) == 1
    assert len(c2.parent.polymers) == 1
    assert str(c2) == str(ComplexSpecies([pc.polymers[0][0], Species("S")], called_from_complex = True))


def test_complex_with_multiple_polymers():
    a = Species('A')
    b = Species('B')
    p = OrderedPolymerSpecies([a,b,a])
    pc0 = PolymerConformation(polymer = p)
    p2 = OrderedPolymerSpecies([b, a, b])
    pc2 = PolymerConformation(polymer = p2)

    #Bind one two monomers from one polymer
    #Polymers must be placed into a Conformation before being bound together
    with pytest.raises(TypeError):
        c = Complex([Species("S"), p[0], p[2]])

    c = Complex([Species("S"), pc0.polymers[0][0], pc0.polymers[0][2]])

    #Correct parent
    assert c.parent == PolymerConformation([ComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)])
    #correct complex returned
    assert str(c) == str(ComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)) 

    #Ordered Case
    oc = Complex([Species("S"), pc0.polymers[0][0], pc0.polymers[0][2]], ordered = True)
    #Correct parent
    assert oc.parent == PolymerConformation([OrderedComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)])
    #correct complex returned
    assert str(oc) == str(OrderedComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True))

    #Two polymers which bind together
    #Polymers must be placed into a Conformation before being bound together
    with pytest.raises(TypeError):
        c2 = Complex([Species("S"), p[0], p2[1]])

    c2 = Complex([Species("S"), pc0.polymers[0][0], pc2.polymers[0][1]]) 
    assert c2.parent == PolymerConformation([ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True)])
    assert str(c2) == str(ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True))

    #Two polymers already bound together in a Conformation interacting at a new location
    l1 = c2.parent.polymers[0][1]
    l2 = c2.parent.polymers[1][0]
    c3 = Complex([l1, l2], ordered = True)
    assert c3.parent == PolymerConformation([
        ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True),
        OrderedComplexSpecies([p[1], p2[0]], called_from_complex = True)
        ])
    assert str(c3) == str(OrderedComplexSpecies([p[1], p2[0]], called_from_complex = True))

def test_complex_with_polymer_replacement():
    #These tests show how the order of binding can matter. 
    #Best practices is to put everything into a PolymerConformation before doing Complex if PolymerConformations are being used.

    a = Species('A')
    b = Species('B')
    s = Species("S")
    p = OrderedPolymerSpecies([a,b,a])
    pc0 = PolymerConformation(polymer = p) #Put the polymer inside a conformation
    #Bind two monomers together
    pc = Complex([pc0.polymers[0][0], pc0.polymers[0][1]]).parent #get a PolymerConformation

    assert isinstance(pc, PolymerConformation)

    #create a Complex around an unbound element of p (from within pc)
    c = Complex([s, pc.polymers[0][2], s], ordered = True)
    assert str(c) == str(OrderedComplexSpecies([s, p[2], s], called_from_complex = True))
    assert str(pc.polymers[0]) == str(c.parent.polymers[0])
    assert len(c.parent.complexes) > len(pc.complexes)

    #Make the same thing with the replacement first
    #this gives a different final complex than doing things in the previous order
    p2 = OrderedPolymerSpecies([a,b,Complex([s, a, s], ordered = True)])
    pc2 = PolymerConformation(polymer = p2)
    assert str(c.parent) != str(p2)
    assert c.parent.parent != Complex([pc2.polymers[0][0], pc2.polymers[0][1]]).parent

def test_complex_with_a_complex_in_a_conformation():
    #This occurs when Complexes are formed around Complexes in PolymerConformations.
    #In these cases, the Complexes are merged to prevent nested Complexes inside of PolymerConformations.
    #Using ordered = True to test that order is preserved
    a = Species('A')
    b = Species('B')
    c = Species('C')
    p = OrderedPolymerSpecies([a,b,c])
    pc0 = PolymerConformation(polymer = p)
    pc = Complex([pc0.polymers[0][0], pc0.polymers[0][1]], ordered = True).parent #get a PolymerConformation
    c = pc.complexes[0] #get a complex from the PolymerConformation

    c2 = Complex([c, Species("S")], ordered = True) #Create a Complex with a Complex
    pc2 = c2.parent
    assert str(c2) == str(OrderedComplexSpecies([p[0], p[1], Species("S")])) #merging done correctly
    assert pc2 == Complex([pc0.polymers[0][0], pc0.polymers[0][1], Species("S")], ordered = True).parent #check the parent PolymerConformation
    assert len(pc2.complexes) == 1

    #Create a PolymerConformation with two complexes
    c3 = Complex([pc2.polymers[0][0], pc2.polymers[0][2], Species("S2")], ordered = True)
    pc3 = c3.parent
    assert len(pc3.complexes) == 2
    assert str(c3) == str(OrderedComplexSpecies([p[0], p[2], Species("S2")], called_from_complex = True))

    #merge the two complexes in pc3
    c4 = Complex([pc3.complexes[0], pc3.complexes[1], Species("S3")], ordered = True)
    assert len(c4.parent.complexes) == 1
    assert str(c4) == str(Complex([pc0.polymers[0][0], pc0.polymers[0][1], Species("S"), pc0.polymers[0][0], pc0.polymers[0][2], Species("S2"), Species("S3")], ordered = True))
    assert c4.parent == Complex([pc0.polymers[0][0], pc0.polymers[0][1], Species("S"), pc0.polymers[0][0], pc0.polymers[0][2], Species("S2"), Species("S3")], ordered = True).parent

    
def test_invalid_complex():
    with pytest.raises(TypeError):
        Complex("A")
    with pytest.raises(TypeError):
        Complex(species = "A")
