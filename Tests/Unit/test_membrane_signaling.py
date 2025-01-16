
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex

#Test Membrane Transport Mechanisms
from biocrnpyler import Membrane_Signaling_Pathway_MM

class test_membrane_signaling_MM():
    tcs = Membrane_Signaling_Pathway_MM()
    MSP = Species("MSP1")
    MSP.ATP=2
    RP = Species("RP1")
    sub_assign = Species("S1")
    sub_signal = Species("S2")
    energy = Species("E1")
    waste = Species("W1")
    c1 =  Complex([sub_signal, MSP])
    c2 =  Complex([MP.ATP*[energy], c1])
    c3 =  Complex([c1, MP.ATP*[waste], sub_assign])
    c4 =  Complex([c1, sub_assign])
    c5 = Complex([c4, RP])
    c6 = Complex([c1, RP, sub_assign])
    c7 =  Complex([RP, sub_assign])

    c_fake = Species("C")
    
    #Test Update Species
    assert len(tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)) == 13
    assert c1 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c2 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c3 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c4 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c5 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c6 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)
    assert c7 in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste)

    assert c_fake in tcs.update_species(MSP, RP, sub_assign, sub_signal, energy, waste, complex = c_fake)
    
    #Test Update Reactions
    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, energy, waste, kb1 = 2e-3, ku1 = 2e-10, 
        kb2 = 2e-3, ku2 = 2e-10, k_hydro = 1e-1, ku3 = 2e-1, kb4 = 2e-3, ku4 = 2e-10, 
        k_phosph = 1e-1, ku5 = 2e-1,ku6 = 2e-10)) == 8

    assert len(tcs.update_reactions(MSP, RP, sub_assign, sub_signal, energy, waste, kb1 = 2e-3, ku1 = 2e-10,
        kb2 = 2e-3, ku2 = 2e-10, k_hydro = 1e-1, ku3 = 2e-1, kb4 = 2e-3, ku4 = 2e-10, 
        k_phosph = 1e-1, ku5 = 2e-1,ku6 = 2e-10, complex_species = c_fake,)) == 8