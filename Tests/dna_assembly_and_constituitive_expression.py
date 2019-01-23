from biocrnpyler import * 

kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg}
myMixture = BasicExtract(name = "txtl", parameters = parameters)

A1 = DNAassembly(name = "G1", promoter = "pBest",
                              rbs = "BCD2", transcript = "T1", protein = "GFP")

#Note: Protein and Transcript strings (or chemical_reaction_network.specie objects) are optional parameters
#DNAassemblies default to using their name for their transcript and protein products.

myMixture.add_components(A1)
myCRN = myMixture.compile_crn()

print("\n"+repr(A1))
print("\n"+repr(myMixture))
print("\n"+repr(myCRN))

import pylab as plt
import numpy as np

x0_dict = {"protein_RNAP":10., "protein_RNAase":10.0, "protein_Ribo":100.,
               'dna_G1':5.}

timepoints = np.arange(0, 4, .001)
print("calling sim with simulate_with_bioscrape")
sim_results_det = myCRN.simulate_with_bioscrape(timepoints, initial_condition_dict = x0_dict)
sim_results_sto = myCRN.simulate_with_bioscrape(timepoints, initial_condition_dict = x0_dict, stochastic = True)


plt.figure()
plt.plot(timepoints, sim_results_det["rna_T1"], color = "blue", label = "T1 (deterministic)")
plt.plot(timepoints, sim_results_det["protein_GFP"], color = "green", label = "GFP (deterministic)")

plt.plot(timepoints, sim_results_sto["rna_T1"], ":", color = "blue", label = "T1 (stochastic)")
plt.plot(timepoints, sim_results_sto["protein_GFP"], ":", color = "green", label = "GFP (stochastic)")
plt.legend()
plt.show()

