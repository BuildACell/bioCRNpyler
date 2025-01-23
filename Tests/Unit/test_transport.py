
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import Species, Complex

#Test Membrane Transport Mechanisms
from biocrnpyler import Simple_Diffusion, Membrane_Protein_Integration, Simple_Transport, Facilitated_Transport_MM, Primary_Active_Transport_MM

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
    assert len(mpi.update_species(MP, IMP)) == 4
    assert c1 in mpi.update_species(MP, IMP)
    assert c_fake in mpi.update_species(MP, IMP, complex = c_fake)
    
    #Test Update Reactions
    assert len(mpi.update_reactions(MP, IMP, kb1 = 1.0, ku1 = 1.0, kb2 = 1.0, ku2 = 1.0, kcat = 1.0, kex = 1.0)) == 2
    assert len(mpi.update_reactions(MP, IMP, kb1 = 1.0, ku1 = 1.0, kb2 = 1.0, ku2 = 1.0, kcat = 1.0, kex = 1.0, complex_species = c_fake,)) == 2

class test_simple_transport():
    st = Simple_Transport()
    MC = Species("MC1")
    MC.attributes=['Passive']
    substrate = Species("S1")
    product = Species("P1")
    
    #Test Update Species
    assert len(st.update_species(MC, substrate, product)) == 3
    
    #Test Update Reactions
    assert len(st.update_reactions(MC, substrate, product, k_channel = 1.0)) == 1

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
    assert len(ft.update_species(MC, substrate, product)) == 5
    assert c1 in ft.update_species(MC, substrate, product)
    assert c2 in ft.update_species(MC, substrate, product)
    assert c_fake in ft.update_species(MC, substrate, product, complex = c_fake)
    
    #Test Update Reactions
    assert len(ft.update_reactions(MC, substrate, product, k1 = 1.0, ku1 = 1.0, ku2 = 1.0, k_trnsp = 1.0)) == 4
    assert len(ft.update_reactions(MC, substrate, product, k1 = 1.0, ku1 = 1.0, ku2 = 1.0, k_trnsp = 1.0, complex_species = c_fake,)) == 4

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
    c2 =  Complex([product, MP])
    c3 =  Complex([MP.ATP*[energy], c1])
    c4 =  Complex([MP.ATP*[waste], c2])
    
    c_fake = Species("C")
    
    #Test Update Species
    assert len(pat.update_species(MP, substrate, product, energy, waste)) == 9
    assert c1 in pat.update_species(MP, substrate, product, energy, waste)
    assert c2 in pat.update_species(MP, substrate, product, energy, waste)
    assert c3 in pat.update_species(MP, substrate, product, energy, waste)
    assert c4 in pat.update_species(MP, substrate, product, energy, waste)

    assert c_fake in pat.update_species(MP, substrate, product, energy, waste, complex = c_fake)
    
    #Test Update Reactions
    assert len(pat.update_reactions(MP, substrate, product, energy, waste, k1 = 1.0, ku1 = 1.0, k2 = 1.0, ku2 = 1.0, k_trnsp = 1.0, ku3= 1.0, ku4= 1.0)) == 7
    assert len(pat.update_reactions(MP, substrate, product, energy, waste, k1 = 1.0, ku1 = 1.0, k2 = 1.0, ku2 = 1.0, k_trnsp = 1.0, ku3= 1.0, ku4= 1.0, complex_species = c_fake,)) == 7