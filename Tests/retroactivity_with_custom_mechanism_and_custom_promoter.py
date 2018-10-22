import txtl
import component
import dna_assembly
import mechanism

#Parameters
kb, ku, ktx, ktl, kdeg = 100, 10, 2.0, 1.0, 5.
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg":kdeg}

#A constituitively expressed reporter
reference_assembly = dna_assembly.DNAassembly(name = "ref", promoter = "P", rbs = "BCD")
#A constiuitively expressed load (RNA and Protein)
full_load_assembly = dna_assembly.DNAassembly(name = "Load", promoter = "P", rbs = "BCD")
#A constiutively transcribed (but not translated) load
RNA_load_assembly = dna_assembly.DNAassembly(name = "Txload", promoter = "P", rbs = None)

#Load genes on orthogonal polymerases
T7 = component.Protein("T7") #Create a new protein T7
#Create a custom promoter with a custom mechanism that uses T7 instead of RNAP
T7P = dna_assembly.Promoter("T7P", mechanisms={
    "transcription":mechanism.Transcription_MM(name = "T7_transcription_mm", rnap=T7)})
#A load assembly with the custom T7 promoter
T7_load_assembly = dna_assembly.DNAassembly(name = "T7Load", promoter = T7P, rbs = "BCD")

#Each new assembly requires its own promoter instance - so here I create another one
T7P = dna_assembly.Promoter("T7P", mechanisms={
    "transcription":mechanism.Transcription_MM(name = "T7_transcription_mm", rnap=T7)})
#A load assembly with the custom T7 promoter and no RBS
T7RNA_load_assembly = dna_assembly.DNAassembly(name = "T7TxLoad", promoter = T7P, rbs = None)

#Add all the assemblies to a mixture
components = [reference_assembly, full_load_assembly, T7_load_assembly, T7, RNA_load_assembly, T7RNA_load_assembly]
myMixture = txtl.BasicExtract(name = "txtl", parameters = parameters, components = components)


print("\nMixture with Assembly\n", repr(myMixture))
print("\n")
for comp in components:
    print(repr(comp))
    print("Species:", repr(comp.update_species()))
    print("Reactions:")
    for rxn in comp.update_reactions():
        print("\t",repr(rxn))
    print("\n")


myCRN = myMixture.compile_crn()
print("species in my crn", len(myCRN.species))
print("reactions in my crn", len(myCRN.reactions))
print("\n"+repr(myCRN)+"\n")

#TODO ADD Code to generate below plots with BioSCRAPE
"""
THE BELOW CODE USES A SSA SIMULATOR DEVELOPED BY WILLIAM POOLE
species, rxns = myCRN.pyrepr()
print(species, len(species))
for rxn in rxns:
    for s in rxn[0]:
        if str(s) not in species:
            print(s, "not in species!", rxn[0], str(s) in species)
simCRN = CRN_Simulator.CRN(species, rxns)
steps = 300000

import pylab as plt
plt.figure(figsize = (16, 8))
plt.subplot(221)
plt.title("Load on a RNAP Promoter")
for dna_Load in [0, 1.0, 2.0, 5., 10., 25., 50, 100]:
    print("Simulating for dna_Load=", dna_Load)
    x0_dict = {"protein_T7": 10., "protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_ref':5., 'dna_Load':dna_Load}
    x0 = myCRN.initial_condition_vector(x0_dict)

    CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)
    plt.plot(t, CD["protein_ref"], label = "Load = "+str(dna_Load))

plt.xlim(0, 100)
plt.xlabel("time")
plt.ylabel("Reference Protein")
plt.legend()

plt.subplot(222)
plt.title("Load on a T7 Promotoer")
for dna_Load in [0, 1.0, 5., 10., 25., 50, 100]:
    print("Simulating for dna_T7Load=", dna_Load)
    x0_dict = {"protein_T7": 10., "protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_ref':5., 'dna_T7Load':dna_Load}
    x0 = myCRN.initial_condition_vector(x0_dict)

    CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)
    plt.plot(t, CD["protein_ref"], label = "Load = "+str(dna_Load))
plt.xlim(0, 100)
plt.xlabel("time")
plt.ylabel("Reference Protein")
plt.legend()

plt.subplot(223)
plt.title("Load on a RNAP Promotoer, No RBS")
for dna_Load in [0, 1.0, 2.0, 5., 10., 25., 50, 100]:
    print("Simulating for dna_TxLoad=", dna_Load)
    x0_dict = {"protein_T7": 10., "protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_ref':5., 'dna_TxLoad':dna_Load}
    x0 = myCRN.initial_condition_vector(x0_dict)

    CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)
    plt.plot(t, CD["protein_ref"], label = "Load = "+str(dna_Load))
plt.xlim(0, 100)
plt.xlabel("time")
plt.ylabel("Reference Protein")
plt.legend()

plt.subplot(224)
plt.title("Load on a T7 Promotoer, No RBS")
for dna_Load in [0, 1.0, 2.0, 5., 10., 25., 50, 100]:
    print("Simulating for dna_T7TxLoad=", dna_Load)
    x0_dict = {"protein_T7": 10., "protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_ref':5., 'dna_T7TxLoad':dna_Load}
    x0 = myCRN.initial_condition_vector(x0_dict)

    CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)
    plt.plot(t, CD["protein_ref"], label = "Load = "+str(dna_Load))
plt.xlim(0, 100)
plt.xlabel("time")
plt.ylabel("Reference Protein")
plt.legend()

plt.show()
"""