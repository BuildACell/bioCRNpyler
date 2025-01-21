#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import DiffusibleMolecule, IntegralMembraneProtein, MembraneChannel, MembranePump, MembraneSensor
import pytest

def test_DiffusibleMolecule():
    diffusion_molecule = 'DP'

    dm = DiffusibleMolecule(substrate=diffusion_molecule)
    assert diffusion_molecule == dm.substrate.name
    assert diffusion_molecule== dm.product.name

    assert dm.get_species().name == 'DP'

    with pytest.raises(KeyError, match='Unable to find mechanism of type diffusion in Component'):
        dm.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type diffusion in Component'):
        dm.update_reactions()

def test_IntegralMembraneProtein():
    membrane_protein = 'MP1'
    products = 'P1'

    imp = IntegralMembraneProtein(membrane_protein=membrane_protein, product=products)
    assert membrane_protein == imp.membrane_protein.name
    assert products== imp.product.name

    assert imp.get_species().name == 'MP1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type membrane_insertion in Component'):
        imp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type membrane_insertion in Component'):
        imp.update_reactions()

def test_MembraneChannel():
    integral_membrane_protein = 'IMP1'
    substrates = 'S1'

    mc = MembraneChannel(integral_membrane_protein, substrate=substrates)
    assert integral_membrane_protein == mc.integral_membrane_protein.name
    assert substrates == mc.substrate.name
    assert substrates == mc.product.name

    assert mc.get_species().name == 'IMP1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type transport in Component'):
        mc.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type transport in Component'):
        mc.update_reactions()

def test_MembranePump(): 
    membrane_pump = 'MPump1'
    substrates = 'S1'

    mp = MembranePump(membrane_pump, substrate=substrates)
    assert membrane_pump == mp.membrane_pump.name
    assert substrates == mp.substrate.name
    assert substrates == mp.product.name

    assert mp.get_species().name == 'MPump1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type transport in Component'):
        mp.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type transport in Component'):
        mp.update_reactions()

def test_MembraneSensor():
    membrane_sensor = 'MSensor1'
    response_protein = 'RP1'
    assigned_substrate = 'Sub_A1'
    signal_substrate = 'Sub_S1'

    ms = MembraneSensor(membrane_sensor, response_protein=response_protein, 
                        assigned_substrate=assigned_substrate, signal_substrate=signal_substrate)
    assert response_protein == ms.response_protein.name
    assert assigned_substrate == ms.assigned_substrate.name
    assert signal_substrate == ms.signal_substrate.name

    assert ms.get_species().name == 'MSensor1'

    with pytest.raises(KeyError, match='Unable to find mechanism of type membrane_sensor in Component'):
            ms.update_species()

    with pytest.raises(KeyError, match='Unable to find mechanism of type membrane_sensor in Component'):
            ms.update_reactions()