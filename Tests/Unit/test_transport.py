
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex

#Test Membrane Transport Mechanisms
from biocrnpyler import Membrane_Protein_Integration, Passive_Transport, Facilitated_Passive_Transport, Primary_Active_Transport

class test_membrane_integration():
    mpi = Membrane_Protein_Integration()
    MP = Species("MP1")
    MP.size=2
    IMP = Species("IMP1")
    c1 = Complex([MP]*MP.size)
    c_fake = Species("C")
    
    #Test Update Species
    assert len(mpi.update_species(MP, IMP)) == 4
    assert c1 in mpi.update_species(MP, IMP)
    assert c_fake in mpi.update_species(MP, IMP, complex = c_fake)
    
    #Test Update Reactions
    assert len(mpi.update_reactions(MP, IMP, kb = 1.0, kd = 1.0, ku = 1.0, kcat = 1.0, kex = 1.0)) == 4
    assert len(mpi.update_reactions(MP, IMP, kb = 1.0, kd = 1.0, ku = 1.0, kcat = 1.0, kex = 1.0, complex_species = c_fake,)) == 4

class test_passive_transport():
    pt = Passive_Transport2()
    MC = Species("MC1")
    substrate = Species("S1")
    product = Species("P1")

    c_fake = Species("C")
    
    #Test Update Species
    assert len(pt.update_species(MC, substrate, product)) == 5
    assert c_fake in pt.update_species(MC, substrate, product, complex = c_fake)
    
    #Test Update Reactions
    assert len(pt.update_reactions(MC, substrate, product, k_channel = 1.0)) == 1
    assert len(pt.update_reactions(MC, substrate, product, k_channel = 1.0, complex_species = c_fake,)) == 1

class test_facilitated_passive_transport():
    fpt = Facilitated_Passive_Transport2()
    MC = Species("MC1")
    substrate = Species("S1")
    product = Species("P1")
    c1 =  Complex([substrate, MC])
    c2 =  Complex([product, MC])
    c_fake = Species("C")
    
    #Test Update Species
    assert len(fpt.update_species(MC, substrate, product)) == 5
    assert c1 in fpt.update_species(MC, substrate, product)
    assert c2 in fpt.update_species(MC, substrate, product)
    assert c_fake in fpt.update_species(MC, substrate, product, complex = c_fake)
    
    #Test Update Reactions
    assert len(fpt.update_reactions(MC, substrate, product, kb = 1.0, ku = 1.0, k1 = 1.0)) == 4
    assert len(fpt.update_reactions(MC, substrate, product, kb = 1.0, ku = 1.0, k1 = 1.0, complex_species = c_fake,)) == 4

class test_active_transport():
    pat = Primary_Active_Transport2()
    MP = Species("MC1")
    MP.ATP=2
    substrate = Species("S1")
    product = Species("P1")
    energy = Species("E1")
    waste = Species("W1")
    c1 =  Complex([substrate, MP])
    c2 =  Complex([product, MP])
    c3 =  Complex([MP.ATP*[energy], c1])
    c4 =  Complex([MP.ATP*[waste], c2])
    
    c_fake = Species("C")
    
    #Test Update Species
    assert len(pat.update_species(MP, substrate, product, energy, waste)) == 7
    assert c1 in pat.update_species(MP, substrate, product, energy, waste)
    assert c2 in pat.update_species(MP, substrate, product, energy, waste)
    assert c3 in pat.update_species(MP, substrate, product, energy, waste)
    assert c4 in pat.update_species(MP, substrate, product, energy, waste)

    assert c_fake in pat.update_species(MP, substrate, product, energy, waste, complex = c_fake)
    
    #Test Update Reactions
    assert len(pat.update_reactions(MP, substrate, product, energy, waste, kb = 1.0, ku = 1.0, kcat = 1.0, kex = 1.0)) == 7
    assert len(pat.update_reactions(MP, substrate, product, energy, waste, kb = 1.0, ku = 1.0, kcat = 1.0, kex = 1.0, complex_species = c_fake,)) == 7

