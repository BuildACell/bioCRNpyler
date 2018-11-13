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

x0_dict = {"protein_RNAP":10., "protein_RNAase":10.0, "complex_Ribo":100.,
               'dna_G1':5.}

timepoints = np.arange(0, 10, .01)
print("calling sim with simulate_with_bioscrape")
sim_results, model = myCRN.simulate_with_bioscrape(timepoints, initial_condition_dict = x0_dict)
results = sim_results.py_get_result()

T1_ind = model.get_species_index("rna_T1")
GFP_ind = model.get_species_index("protein_GFP")
plt.figure()
plt.plot(timepoints, results[:, T1_ind], label = "T1")
plt.plot(timepoints, results[:, GFP_ind], label = "GFP")
plt.legend()
plt.show()

