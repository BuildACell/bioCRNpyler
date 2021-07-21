#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Enzyme
import pytest


def test_enzyme():
    substrates = 'S1'
    products = 'P1'
    enzyme = 'E1'

    e = Enzyme(enzyme=enzyme, substrates=substrates, products=products)
    assert any([s.name == 'S1' for s in e.substrates])
    assert any([p.name == 'P1' for p in e.products])
    assert enzyme == e.enzyme.name

    assert e.get_species().name == 'E1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        e.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        e.update_reactions()
