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
    assert(d_c.parent,truth)

    #Ordered case
    print("d", d)
    d_co = Complex([d[1],a], ordered = True)
    print("d", d)
    print("p", d_co.parent)
    assert isinstance(d_co, OrderedComplexSpecies)
    truth_o = OrderedPolymerSpecies([a, Complex([b,a], ordered = True),a])
    assert d_co.parent==truth_o

def test_complex_a_whole_polymer():
    a = Species('A')
    b = Species('B')
    p = OrderedPolymerSpecies([a,b,a])
    c = Complex([Species("S"), p[0], p[2]])
    pc = c.parent

    #This works as expected
    c2 = Complex([p, Species("S")])
    assert c2 == ComplexSpecies([p, Species("S")], called_from_complex = True)

    #If the Polymer is in a Conformation, it cannot be Complexed.
    with pytest.raises(TypeError):
        c2 = Complex([pc.polymers[0], Species("S")])
    

def test_complex_with_multiple_polymers():
    a = Species('A')
    b = Species('B')
    p = OrderedPolymerSpecies([a,b,a])
    p2 = OrderedPolymerSpecies([b, a, b])

    #Bind a single OrderedPolymerSpecies into a PolymerConformation
    c = Complex([Species("S"), p[0], p[2]])
    #Correct parent
    assert c.parent == PolymerConformation([ComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)])
    #correct complex returned
    assert str(c) == str(ComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)) 

    #Ordered Case
    oc = Complex([Species("S"), p[0], p[2]], ordered = True)
    #Correct parent
    assert oc.parent == PolymerConformation([OrderedComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True)])
    #correct complex returned
    assert str(oc) == str(OrderedComplexSpecies([Species("S"), p[0], p[2]], called_from_complex = True))

    #Two polymers which bind together
    c2 = Complex([Species("S"), p[0], p2[1]])
    assert c2.parent == PolymerConformation([ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True)])
    assert str(c2) == str(ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True))

    #Two polymers already bound together in a Conformation interacting at a new location
    l1 = c2.parent.polymers[0][1]
    print(l1, l1.parent, l1.parent.parent)
    l2 = c2.parent.polymers[1][0]
    print(l2, l2.parent, l2.parent.parent)
    c3 = Complex([l1, l2], ordered = True)
    assert c3.parent == PolymerConformation([
        ComplexSpecies([Species("S"), p[0], p2[1]], called_from_complex = True),
        OrderedComplexSpecies([p[1], p2[0]], called_from_complex = True)
        ])
    assert str(c3) == str(OrderedComplexSpecies([p[1], p2[0]], called_from_complex = True))

def test_complex_with_polymer_replacement():
    #This occurs when a Complex is formed around a Monomer in a Polymer in a PolymerConformation which does not include any other Monomers with parents.
    a = Species('A')
    b = Species('B')
    p = OrderedPolymerSpecies([a,b,a])
    pc = Complex([p[0], p[1]]).parent #get a PolymerConformation
    #create a Complex around an unbound element of p (from within pc)
    c = Complex([Species("S"), pc.polymers[0][2], Species("S")], ordered = True)
    assert str(c) == str(OrderedComplexSpecies([Species("S"), p[2], Species("S")], called_from_complex = True))

    #Make the same thing with the replacement first
    p2 = OrderedPolymerSpecies([a,b,Complex([Species("S"), a, Species("S")], ordered = True)])
    assert str(c.parent) == str(p2)
    assert c.parent.parent == Complex([p2[0], p2[1]]).parent
    

def test_complex_with_a_complex_in_a_conformation():
    #This occurs when Complexes are formed around Complexes in PolymerConformations.
    #In these cases, the Complexes are merged to prevent nested Complexes inside of PolymerConformations.
    #Using ordered = True to test that order is preserved
    a = Species('A')
    b = Species('B')
    c = Species('C')
    p = OrderedPolymerSpecies([a,b,c])
    pc = Complex([p[0], p[1]], ordered = True).parent #get a PolymerConformation
    c = pc.complexes[0] #get a complex from the PolymerConformation

    c2 = Complex([c, Species("S")], ordered = True) #Create a Complex with a Complex
    pc2 = c2.parent
    assert str(c2) == str(OrderedComplexSpecies([p[0], p[1], Species("S")])) #merging done correctly
    assert pc2 == Complex([p[0], p[1], Species("S")], ordered = True).parent #check the parent PolymerConformation
    assert len(pc2.complexes) == 1

    #Create a PolymerConformation with two complexes
    c3 = Complex([pc2.polymers[0][0], pc2.polymers[0][2], Species("S2")], ordered = True)
    pc3 = c3.parent
    assert len(pc3.complexes) == 2
    assert str(c3) == str(OrderedComplexSpecies([p[0], p[2], Species("S2")], called_from_complex = True))

    #merge the two complexes in pc3
    c4 = Complex([pc3.complexes[0], pc3.complexes[1], Species("S3")], ordered = True)
    assert len(c4.parent.complexes) == 1
    assert str(c4) == str(Complex([p[0], p[1], Species("S"), p[0], p[2], Species("S2"), Species("S3")], ordered = True))
    assert c4.parent == Complex([p[0], p[1], Species("S"), p[0], p[2], Species("S2"), Species("S3")], ordered = True).parent

    
def test_invalid_complex():
    with pytest.raises(TypeError):
        Complex("A")
    with pytest.raises(TypeError):
        Complex(species = "A")
