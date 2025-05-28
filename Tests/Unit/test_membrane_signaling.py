
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex, ParameterKey, Component

#Test Membrane Transport Mechanisms
from biocrnpyler import Membrane_Signaling_Pathway_MM

def contains(element, nested_array):
    """Recursively checks if an element is in a nested list."""
    return any(
        contains(element, sublist) if isinstance(sublist, list) else element == sublist
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

class test_membrane_signaling_MM():
    tcs = Membrane_Signaling_Pathway_MM()
    MSP = Species("MSP1")
    MSP.ATP=2
    RP = Species("RP1")
    sub_assign = Species("S1")
    sub_signal = Species("S2")
    product= Species('RP_active')
    energy = Species("E1")
    waste = Species("W1")
    
    #Create empty dictionary for complexes
    complex_dict={}
    #Complex1
    complex_dict['Activated_MP']=Complex([sub_signal, MSP])
    #Complex2
    complex_dict['ATP:Activated_MP']=Complex([MSP.ATP*[energy], complex_dict['Activated_MP']])
    #Complex3
    complex_dict['ADP:Activated_MP:Sub']=Complex([complex_dict['Activated_MP'], MSP.ATP*[waste], sub_assign])
    #Complex4
    complex_dict['Activated_MP:Sub']=Complex([complex_dict['Activated_MP'], sub_assign])
    #Complex5
    complex_dict['Activated_MP:Sub:RP']=Complex([complex_dict['Activated_MP:Sub'], RP])
    #Complex6
    complex_dict['Activated_MP:RP:Sub']=Complex([complex_dict['Activated_MP'], RP, sub_assign])
    
    #Test Update Species
    assert total_length(tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste)) == 12
    assert contains(complex_dict['Activated_MP'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    assert contains(complex_dict['ATP:Activated_MP'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    assert contains( complex_dict['ADP:Activated_MP:Sub'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    assert contains(complex_dict['Activated_MP:Sub'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    assert contains(complex_dict['Activated_MP:Sub:RP'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    assert contains(complex_dict['Activated_MP:RP:Sub'], tcs.update_species(MSP, RP, sub_assign, sub_signal, product, energy, waste))
    
    #Test Update Reactions
    #Define sensor parameter dictionary and component
    sensor_param_dict = {
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "kb_sigMS"):2e-3, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_sigMS"):2e-10, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "kb_autoPhos"): 2e-3, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_autoPhos"):2e-10,  
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "k_hydro"):1e-1,
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_waste"):1e-1, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "kb_phosRP"):2e-3, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_phosRP"): 2e-10, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "k_phosph"):1e-1, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_activeRP"):2e-1, 
        ParameterKey(mechanism = "two_component_membrane_signaling", part_id = None, name = "ku_dephos"):2e-10,}
    sensor_params=Component("sensor_params",parameters = sensor_param_dict)

    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, product, energy, waste, component=sensor_params)) == 8

    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, product, energy, waste, component=sensor_params)) == 8