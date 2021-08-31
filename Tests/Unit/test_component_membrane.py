#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import IntegralMembraneProtein
import pytest


def test_IntegralMembraneProtein():
    membrane_protein = 'MP1'
    products = 'P1'

    mp = IntegralMembraneProtein2(membrane_protein=membrane_protein, product=products)
    assert membrane_protein == mp.membrane_protein.name
    assert products== mp.product.name

    assert mp.get_species().name == 'MP1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        mp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type catalysis in Component'):
        mp.update_reactions()
