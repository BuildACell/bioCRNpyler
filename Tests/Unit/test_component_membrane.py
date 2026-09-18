#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest

from biocrnpyler import (
    DiffusibleMolecule,
    IntegralMembraneProtein,
    MembraneChannel,
    MembraneCarrier,
    MembranePump,
    MembraneSensor,
)

def test_DiffusibleMolecule():
    diffusion_molecule = 'DP'

    dm = DiffusibleMolecule(substrate=diffusion_molecule)

    # Iterate over list
    subs_name= []
    for sub in dm.substrate:
        subs_name.append(sub.name)
    prod_name= []
    for prod in dm.product:
        prod_name.append(prod.name)
    species_name= []
    for species in dm.get_species():
        species_name.append(species.name)

    assert ['DP'] == subs_name
    assert ['DP'] == prod_name

    assert species_name == ['DP']

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion in Component',
    ):
        dm.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion in Component',
    ):
        dm.update_reactions()

def test_IntegralMembraneProtein():
    membrane_protein = 'MP1'
    products = 'P1'

    imp = IntegralMembraneProtein(
        membrane_protein=membrane_protein, product=products
    )
    assert membrane_protein == imp.membrane_protein.name
    assert products == imp.product.name

    assert imp.get_species().name == 'MP1'

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type membrane_integration in Component',
    ):
        imp.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type membrane_integration in Component',
    ):
        imp.update_reactions()

def test_MembraneChannel():
    membrane_channel = 'IMP1'
    substrates = 'S1'

    mc = MembraneChannel(membrane_channel, substrate=substrates)
    assert membrane_channel == mc.membrane_channel.name

    # Iterate over list
    subs_name= []
    for sub in mc.substrate_in:
        subs_name.append(sub.name)
    prod_name= []
    for prod in mc.substrate_out:
        prod_name.append(prod.name)

    assert [substrates] == subs_name
    assert [substrates] == prod_name

    assert mc.get_species().name == 'IMP1'

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion in Component',
    ):
        mc.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion in Component',
    ):
        mc.update_reactions()

def test_MembraneCarrier():
    membrane_carrier = 'IMP1'
    substrates = 'S1'

    mc = MembraneCarrier(membrane_carrier, substrate=substrates)
    assert membrane_carrier == mc.membrane_carrier.name

    # Iterate over list
    subs_name= []
    for sub in mc.substrate_in:
        subs_name.append(sub.name)
    prod_name= []
    for prod in mc.substrate_out:
        prod_name.append(prod.name)

    assert [substrates] == subs_name
    assert [substrates] == prod_name

    assert mc.get_species().name == 'IMP1'

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion or transport in Component',
    ):
        mc.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type diffusion or transport in Component',
    ):
        mc.update_reactions()

def test_MembranePump():
    membrane_pump = 'MPump1'
    substrates = 'S1'

    mp = MembranePump(membrane_pump, substrate=substrates, direction='exporter')
    assert membrane_pump == mp.membrane_pump.name

    # Iterate over list
    subs_name= []
    for sub in mp.substrate:
        subs_name.append(sub.name)
    prod_name= []
    for prod in mp.product:
        prod_name.append(prod.name)

    assert [substrates] == subs_name
    assert [substrates] == prod_name

    assert mp.get_species().name == 'MPump1'

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type transport in Component',
    ):
        mp.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type transport in Component',
    ):
        mp.update_reactions()

def test_MembraneSensor():
    membrane_sensor = 'MSensor1'
    response_protein = 'RP1'
    assigned_substrate = 'Sub_A1'
    signal_substrate = 'Sub_S1'

    ms = MembraneSensor(
        membrane_sensor,
        response_protein=response_protein,
        assigned_substrate=assigned_substrate,
        signal_substrate=signal_substrate,
    )
    assert response_protein == ms.response_protein.name
    assert assigned_substrate == ms.assigned_substrate.name
    assert signal_substrate == ms.signal_substrate.name

    assert ms.get_species().name == 'MSensor1'

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type membrane_sensor in Component',
    ):
        ms.update_species()

    with pytest.raises(
        KeyError,
        match='Unable to find mechanism of type membrane_sensor in Component',
    ):
        ms.update_reactions()
