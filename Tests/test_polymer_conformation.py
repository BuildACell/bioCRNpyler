#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import *

from biocrnpyler import OrderedPolymerSpecies, PolymerConformation, Complex, Species

def test_single_polymer_single_complex_instantiation():
    p = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    c1 = ComplexSpecies([p[0], p[1]], called_from_complex = True)
    pc = PolymerConformation(polymers = [p], complexes = [c1], called_from_complex = True)
    #Test naming convention
    assert str(pc) == f"conformation__{p}_0_0_{c1}_"

    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    c2 = ComplexSpecies([p[0], p2[0]], called_from_complex = True)

    #In these cases, TypeErrors should be raised
    with pytest.raises(TypeError):
        pc = PolymerConformation([], [], called_from_complex = True)
    with pytest.raises(TypeError):
        pc = PolymerConformation([p], [], called_from_complex = True)
    with pytest.raises(TypeError):
        pc = PolymerConformation([], [c1], called_from_complex = True)

    #In these cases, ValueErrors should be raised.
    with pytest.raises(ValueError):
        #p2 should be in the Polymers list
        pc = PolymerConformation(polymers = [p], complexes = [c2], called_from_complex = True)



def test_multiple_polymer_multiple_complex_instantiation():
    p1 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m3")])
    p2 = OrderedPolymerSpecies([Species("m1"), Species("m2"), Species("m4")])
    c1 = ComplexSpecies([p1[0], p1[1]], called_from_complex = True)
    c2 = ComplexSpecies([p1[0], p2[2]], called_from_complex = True)
    c3 = ComplexSpecies([Species("S1"), Species("S2")], called_from_complex = True)
    
    pc1 = PolymerConformation(polymers = [p1, p2], complexes = [c1, c2], called_from_complex = True)
    
    #Check naming convention
    assert str(pc1) == f"conformation__{p1}_{p2}_0_0_{c1}_0_1_{c2}_"

    #check that order doesn't matter
    pc1_r = PolymerConformation(polymers = [p2, p1], complexes = [c2, c1], called_from_complex = True)
    assert pc1 == pc1_r

    #check that ComplexSpecies with the same string representation but different polymers are treated differently
    c1b = ComplexSpecies([p2[1], p1[0]], called_from_complex = True) #c1 and c1b have the same string representation - but connect different polymers
    pc1b = PolymerConformation(polymers = [p1, p2], complexes = [c1b, c2], called_from_complex = True)
    assert str(pc1b) == f"conformation__{p1}_{p2}_0_1_{c1}_0_1_{c2}_"
    assert pc1 != pc1b

    pc2 = PolymerConformation(polymers = [p1, p2], complexes = [c1, c2], called_from_complex = True)
    pc2b = PolymerConformation(polymers = [p1, p2], complexes = [c1b, c2], called_from_complex = True)
    assert pc2 != pc2b

    #In these cases, value errors should be raised
    with pytest.raises(ValueError): 
        #p2 should not be in the Polymer list because no monomers have it as a parent
        pc = PolymerConformation(polymers = [p1, p2], complexes = [c1], called_from_complex = True)

    with pytest.raises(ValueError):
        #c3 is not part of the conformation
        pc = PolymerConformation(polymers = [p1], complexes = [c1, c3], called_from_complex = True)


