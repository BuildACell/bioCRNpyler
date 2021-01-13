#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import *

from biocrnpyler import OrderedPolymerSpecies, PolymerConformation, Complex, Species

def test_single_polymer_single_complex_instantiation():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p[0], p[1]], called_from_complex = True)
    pc = PolymerConformation(complexes = [c1], called_from_complex = True)
    #Test naming convention
    assert str(pc) == f"conformation__{p}_0_0_{c1}_"

    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    c2 = ComplexSpecies([p[0], p2[0]], called_from_complex = True)

    #In these cases, TypeErrors should be raised
    with pytest.raises(TypeError):
        #Emtpy complexes is not allowed
        pc2 = PolymerConformation(complexes = [], called_from_complex = True)

    with pytest.raises(ValueError):
        #set a PolymerConformation as a parent of the Species.
        S = Species("S")
        S.parent = pc
        c3 = ComplexSpecies([p[0], S], called_from_complex = True)
        pc2 = PolymerConformation(complexes = [c3])


def test_multiple_polymer_multiple_complex_instantiation():
    p1 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    

    c1 = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    c2 = ComplexSpecies([p1[0], p2[2]], called_from_complex = True)
    c3 = ComplexSpecies([Species("S1"), Species("S2")], called_from_complex = True)
    
    
    pc1 = PolymerConformation(complexes = [c1, c2], called_from_complex = True)
    
    #Check naming convention
    assert str(pc1) == f"conformation__{p1}_{p2}_0_0_{c1}_0_1_{c2}_"

    #check that order doesn't matter
    pc1_r = PolymerConformation(complexes = [c2, c1], called_from_complex = True)
    assert pc1 == pc1_r

    #check that ComplexSpecies with the same string representation but different polymers are treated differently
    c1b = ComplexSpecies([p2[1], p1[0]], called_from_complex = True) #c1 and c1b have the same string representation - but connect different polymers
    pc1b = PolymerConformation(complexes = [c1b, c2], called_from_complex = True)
    assert str(pc1b) == f"conformation__{p1}_{p2}_0_1_{c1}_0_1_{c2}_"
    assert pc1 != pc1b

    pc2 = PolymerConformation(complexes = [c1, c2], called_from_complex = True)
    pc2b = PolymerConformation(complexes = [c1b, c2], called_from_complex = True)
    assert pc2 != pc2b

    #In these cases, value errors should be raised
    with pytest.raises(ValueError):
        #c3 is not part of the conformation
        pc = PolymerConformation(complexes = [c1, c3], called_from_complex = True)

    #try 2 polymers with complexes that look the same
    p3 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c4 = ComplexSpecies([p1[0], p3[0]])

    pc3 = PolymerConformation([c1, c4], called_from_complex = True)
    assert str(pc3) == f"conformation__{p1}_{p1}_0_1_{c4}_0_0_{c1}_"

    with pytest.raises(ValueError):
        #duplicate Complexes are not allowed (same object case)
        pc4 = PolymerConformation([c1, c1], called_from_complex = True)

    c1b = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    with pytest.raises(ValueError):
        #duplicate Complexes are not allowed (identical object case)
        pc4 = PolymerConformation([c1, c1b], called_from_complex = True)



