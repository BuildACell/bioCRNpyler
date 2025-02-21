
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex, ParameterKey, Component

#Test Membrane Transport Mechanisms
from biocrnpyler import Simple_Diffusion, Membrane_Protein_Integration, Simple_Transport, Facilitated_Transport_MM, Primary_Active_Transport_MM

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

class test_simple_diffusion():
    sd = Simple_Diffusion()
    substrate = Species("DMi")
    product = Species("DMo")
    c_fake = Species("C")

    #Test Update Species
    assert len(sd.update_species(substrate, product))==2
    assert substrate in sd.update_species(substrate, product)
    assert product in sd.update_species(substrate, product)

    #Test Update Reactions
    assert len(sd.update_reactions(substrate, product, k_diff = 1.0,)) == 1
    assert len(sd.update_reactions(substrate, product, k_diff = 1.0, complex_species = c_fake,)) == 1

class test_membrane_integration():
    mpi = Membrane_Protein_Integration()
    MP = Species("MP1")
    MP.size=2
    IMP = Species("IMP1")
    c1 = Complex([MP]*MP.size)
    c_fake = Species("C")
    
    #Test Update Species
    assert len(mpi.update_species(MP, IMP)) == 3
    assert c1 in mpi.update_species(MP, IMP)
    assert c_fake in mpi.update_species(MP, IMP, complex = c_fake)
    
    #Test Update Reactions
    #Define sensor parameter dictionary and component
    insertion_param_dict = {
        ParameterKey(mechanism = "membrane_protein_integration", part_id = None, name = "kb_oligmor"):2e-3, 
        ParameterKey(mechanism = "membrane_protein_integration", part_id = None, name = "ku_oligmor"):2e-10, 
        ParameterKey(mechanism = "membrane_protein_integration", part_id = None, name = "kex"): 2e-3, 
        ParameterKey(mechanism = "membrane_protein_integration", part_id = None, name = "kcat"):2e-10,}
    insertion_params=Component("insertion_params",parameters = insertion_param_dict)

    assert len(mpi.update_reactions(MP, IMP, component=insertion_params)) == 2
    assert len(mpi.update_reactions(MP, IMP, component=insertion_params, complex_species = c_fake,)) == 2

class test_simple_transport():
    st = Simple_Transport()
    MC = Species("MC1")
    MC.attributes=['Passive']
    substrate = Species("S1")
    product = Species("P1")
    c_fake = Species("C")
    
    #Test Update Species
    assert len(st.update_species(MC, substrate, product)) == 3
    
    #Test Update Reactions
    assert len(st.update_reactions(MC, substrate, product, k_trnsp = 1.0)) == 1
    assert len(st.update_reactions(MC, substrate, product, k_trnsp = 1.0,complex_species = c_fake,)) == 1

class test_facilitated_transport_MM():
    ft = Facilitated_Transport_MM()
    MC = Species("MC1")
    MC.attributes=['Importer']
    substrate = Species("S1")
    product = Species("P1")
    c1 =  Complex([substrate, MC])
    c2 =  Complex([product, MC])
    c_fake = Species("C")
    
    #Test Update Species
    assert total_length(ft.update_species(MC, substrate, product)) == 5
    assert contains(c1,ft.update_species(MC, substrate, product))
    assert contains(c2 , ft.update_species(MC, substrate, product))
    
    #Test Update Reactions
    #Define sensor parameter dictionary and component
    transport_param_dict = {
        ParameterKey(mechanism = "facilitated_membrane_protein_transport", part_id = None, name = "kb_subMC"):2e-3, 
        ParameterKey(mechanism = "facilitated_membrane_protein_transport", part_id = None, name = "ku_subMC"):2e-10, 
        ParameterKey(mechanism = "facilitated_membrane_protein_transport", part_id = None, name = "k_trnspMC"): 2e-3, 
        ParameterKey(mechanism = "facilitated_membrane_protein_transport", part_id = None, name = "ku_prodMC"):2e-10,}
    transport_params=Component("transport_params",parameters = transport_param_dict)

    #Test Update Reactions
    assert len(ft.update_reactions(MC, substrate, product, component=transport_params)) == 4
    assert len(ft.update_reactions(MC, substrate, product, component=transport_params,complex_species = c_fake,))  == 4

class test_active_transport_MM():
    pat = Primary_Active_Transport_MM()
    MP = Species("MC1")
    MP.ATP=2
    MP.attributes=['Exporter']
    substrate = Species("S1")
    product = Species("P1")
    energy = Species("E1")
    waste = Species("W1")
    c1 =  Complex([substrate, MP])
    c2 =  Complex([MP.ATP*[energy], c1])
    c3 =  Complex([MP.ATP*[energy], product, MP])
    c4 =  Complex([MP.ATP*[waste], MP])
    c_fake = Species("C")
    
    #Test Update Species
    assert total_length(pat.update_species(MP, substrate, product, energy, waste)) == 9
    assert contains(c1, pat.update_species(MP, substrate, product, energy, waste))
    assert contains(c2,  pat.update_species(MP, substrate, product, energy, waste))
    assert contains(c3 , pat.update_species(MP, substrate, product, energy, waste))
    assert contains(c4 ,pat.update_species(MP, substrate, product, energy, waste))
    
    #Test Update Reactions
    #Define sensor parameter dictionary and component
    transport_param_dict = {
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "kb_subMP"):2e-3, 
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "ku_subMP"):2e-10, 
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "kb_subMPnATP"): 2e-3, 
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "ku_subMPnATP"): 2e-1,
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "k_trnspMP"):2e-10,
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "ku_prodMP"):2e-10,
        ParameterKey(mechanism = "active_membrane_protein_transport", part_id = None, name = "ku_MP"):2e-10,}
    transport_params=Component("transport_params",parameters = transport_param_dict)

    assert len(pat.update_reactions(MP, substrate, product, energy, waste, component=transport_params)) == 7
    assert len(pat.update_reactions(MP, substrate, product, energy, waste, component=transport_params,complex_species = c_fake,)) == 7