#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import IntegralMembraneProtein, MembraneChannel, MembranePump
import pytest


def test_IntegralMembraneProtein():
    membrane_protein = 'MP1'
    products = 'P1'

    mp = IntegralMembraneProtein(membrane_protein=membrane_protein, product=products)
    assert membrane_protein == mp.membrane_protein.name
    assert products== mp.product.name

    assert mp.get_species().name == 'MP1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        mp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        mp.update_reactions()

def test_MembraneChannel():
    integral_membrane_protein = 'IMP1'
    substrates = 'S1'

    imp = MembraneChannel(integral_membrane_protein, substrate=substrates)
    assert integral_membrane_protein == imp.integral_membrane_protein.name
    assert substrates == imp.substrate.name
    assert substrates == imp.product.name

    assert imp.get_species().name == 'IMP1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        imp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        imp.update_reactions()

def test_MembranePump(): 
    membrane_pump = 'MPump1'
    substrates = 'S1'

    imp = MembranePump(membrane_pump, substrate=substrates)
    assert membrane_pump == imp.membrane_pump.name
    assert substrates == imp.substrate.name
    assert substrates == imp.product.name

    assert imp.get_species().name == 'MPump1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        imp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        imp.update_reactions()