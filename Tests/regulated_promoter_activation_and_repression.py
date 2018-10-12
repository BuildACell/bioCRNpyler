import txtl
import component
import dna_assembly
import mechanism
import CRN_Simulator
#Parameters
kb, ku, ktx, ktl, kdeg = 200, 10, 2.0, .25, 1.5
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg":kdeg,
              "cooperativity":2,
              ('cooperative_binding', 'repressor', 'kb'):1000, ('cooperative_binding',"repressor", 'ku'):5.0,
              ('transcription', 'repressor', 'kb'):1, ('transcription',"repressor", 'ku'):1000.0,
              ('cooperative_binding', 'activator', 'kb'):1000, ('cooperative_binding', "activator", 'ku'):5.0,
              ('transcription', 'activator', 'kb'): 1000, ('transcription', "activator", 'ku'): 1.0,
              "ktx_leak":ktx, "kb_leak":kb/50, "ku_leak":ku*50}

P_reg = dna_assembly.RegulatedPromoter("P_regulated", regulators=["activator", "repressor"], leak = True)

reg_rep_assembly = dna_assembly.DNAassembly(name = "reporter", promoter = P_reg, rbs = "BCD")

activator = component.Protein("activator")
repressor = component.Protein("repressor")

components = [reg_rep_assembly, activator, repressor]
myMixture = txtl.BasicExtract(name = "txtl", parameters = parameters, components = components)

myCRN = myMixture.compile_crn()




print("\n"+repr(myCRN))

print("\nTotal Species", len(myCRN.species))
print("\nTotal Reactions", len(myCRN.reactions))

import pylab as plt
plt.figure(figsize = (12, 8))

species, rxns = myCRN.pyrepr()
print(species, len(species))
simCRN = CRN_Simulator.CRN(species, rxns)
steps = 500000

print("Simulating with repressor")
x0_dict = {"protein_RNAP":10., "protein_RNAase":20.0, "Ribo":100.,
               'dna_reporter':5., 'protein_activator':0, 'protein_repressor':100}
x0 = myCRN.initial_condition_vector(x0_dict)
CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)
plt.subplot(211)
rna = CD["rna_reporter"]+CD["rna_reporter:complex_Ribo"]+CD["rna_reporter:protein_RNAase"]
plt.plot(t, rna, label = "repressor present")
plt.subplot(212)
plt.plot(t, CD["protein_reporter"], label = "repressor present")

print("Simulating with activator")
x0_dict = {"protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_reporter':5., 'protein_activator':100, 'protein_repressor':0}
x0 = myCRN.initial_condition_vector(x0_dict)
CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)

plt.subplot(211)
rna = CD["rna_reporter"]+CD["rna_reporter:complex_Ribo"]+CD["rna_reporter:protein_RNAase"]
plt.plot(t, rna, label = "activator present")
plt.subplot(212)
plt.plot(t, CD["protein_reporter"], label = "activator present")

print("Simulating with no activator or repressor")
x0_dict = {"protein_RNAP":10., "protein_RNAase":5.0, "Ribo":100.,
               'dna_reporter':5., 'protein_activator':0, 'protein_repressor':0}
x0 = myCRN.initial_condition_vector(x0_dict)
CD, t, rxn_list = simCRN.simulate_cython(x0, steps, return_count_dict=True)

plt.subplot(211)
rna = CD["rna_reporter"]+CD["rna_reporter:complex_Ribo"]+CD["rna_reporter:protein_RNAase"]
plt.plot(t, rna, label = "No activator and no repressor")

plt.subplot(212)
plt.plot(t, CD["protein_reporter"], label = "No activator and no repressor")

plt.subplot(211)
plt.title("Stochastic Simulation of an Activatable/Repressible DNA Assembly")
plt.xlim(0, 200)
plt.ylabel("reporter mRNA")
plt.xlabel("time")
plt.legend()

plt.subplot(212)
plt.xlim(0, 200)
plt.ylabel("reporter protein")
plt.xlabel("time")
plt.legend()
plt.show()