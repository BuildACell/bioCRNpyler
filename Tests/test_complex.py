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

def test_complex_with_polymer_replacement():
    pass

def test_complex_with_a_complex_in_a_conformation():
    pass
    
def test_invalid_complex():
    with pytest.raises(TypeError):
        Complex("A")
    with pytest.raises(TypeError):
        Complex(species = "A")
