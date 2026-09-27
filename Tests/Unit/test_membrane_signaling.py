#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

# Test Membrane Transport Mechanisms
from biocrnpyler import (
    Complex,
    Component,
    Sensor_TwoComponentSystem,
    ParameterKey,
    Species,
)

def contains(element, nested_array):
    """Recursively checks if an element is in a nested list."""
    return any(
        contains(element, sublist)
        if isinstance(sublist, list)
        else element == sublist
        for sublist in nested_array
    )

def total_length(nested_array):
    """Recursively counts the total number of elements in a nested list."""
    count = 0
    for item in nested_array:
        if isinstance(item, list):
            count += total_length(item)  # Recursively count sublist elements
        else:
            count += 1  # Count individual elements
    return count

class test_sensor_twocomponentsystem:
    tcs = Sensor_TwoComponentSystem()
    MS = Species('MS1')
    MS.ATP = 2
    RP = Species('RP1')
    sub_assign = Species('S1')
    sub_signal = Species('S2')
    product = Species('RP_active')
    energy = Species('E1')
    waste = Species('W1')

    # Create empty dictionary for complexes
    complex_dict = {}
    # Complex1
    complex_dict['Activated_MS'] = Complex([sub_signal, MS])
    # Complex2
    complex_dict['ATP:Activated_MS'] = Complex(
        [MS.ATP * [energy], complex_dict['Activated_MS']]
    )
    # Complex3
    complex_dict['ADP:Activated_MS:sub'] = Complex(
        [complex_dict['Activated_MS'], MS.ATP * [waste], sub_assign]
    )
    # Complex4
    complex_dict['Activated_MS:sub'] = Complex(
        [complex_dict['Activated_MS'], sub_assign]
    )
    # Complex5
    complex_dict['Activated_MS:sub:RP'] = Complex(
        [complex_dict['Activated_MS:sub'], RP]
    )
    # Complex7
    complex_dict['Activated_RP'] = Complex(
        [RP, sub_assign]
    )
    # Complex6
    complex_dict['Activated_MS:Activated_RP'] = Complex(
        [complex_dict['Activated_MS'], complex_dict['Activated_RP']]
    )

    # Test Update Species
    assert (
        total_length(
            tcs.update_species(
                MS, RP, sub_assign, sub_signal, product, energy, waste
            )
        )
        == 14
    )
    assert contains(
        complex_dict['Activated_MS'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['ATP:Activated_MS'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['ADP:Activated_MS:sub'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['Activated_MS:sub'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['Activated_MS:sub:RP'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['Activated_RP'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )
    assert contains(
        complex_dict['Activated_MS:Activated_RP'],
        tcs.update_species(
            MS, RP, sub_assign, sub_signal, product, energy, waste
        ),
    )

    # Test Update Reactions
    # Define sensor parameter dictionary and component
    sensor_param_dict = {
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='kb_sigMS',
        ): 2e-3,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_sigMS',
        ): 2e-10,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='kb_autoPhos',
        ): 2e-3,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_autoPhos',
        ): 2e-10,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='k_hydro',
        ): 1e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_waste',
        ): 1e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='kb_phosRP',
        ): 2e-3,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_phosRP',
        ): 2e-10,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='k_phosph',
        ): 1e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_activeRP',
        ): 2e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='kb_activeRP',
        ): 2e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_activeRP',
        ): 2e-1,
        ParameterKey(
            mechanism='sensor_two_component_signaling',
            part_id=None,
            name='ku_dephos',
        ): 2e-10,
    }
    sensor_params = Component('sensor_params', parameters=sensor_param_dict)

    assert (
        len(
            tcs.update_reactions(
                MS,
                RP,
                sub_assign,
                sub_signal,
                product,
                energy,
                waste,
                component=sensor_params,
            )
        )
        == 9
    )

    assert (
        len(
            tcs.update_reactions(
                MS,
                RP,
                sub_assign,
                sub_signal,
                product,
                energy,
                waste,
                component=sensor_params,
            )
        )
        == 9
    )
