from biocrnpyler import *

import bioscrape
import numpy as np
import pylab as plt
import time as pytime

kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg}
myMixture = BasicExtract(name = "txtl", parameters = parameters)

A1 = DNAassembly(name = "G1", promoter = "P",
                              rbs = "BCD", transcript = "T1", protein = "Reporter")

myMixture.add_components(A1)
CRN = myMixture.compile_crn()
print("\n"+repr(CRN))

file_name = "bioscrape_test.xml"
f = CRN.write_sbml_file(file_name)

#Initial Condition Dict: repr(specie) --> concentration. Default is 0
#Note in SBML all the species names have ":" replaced with "_".
initial_condition_dict = {repr(A1.dna):2, "protein_Ribo":10, "protein_RNAP":5, "protein_RNAase":2.5}

#Bioscrape Simulation
timepoints = np.arange(0, 50, .01)
print("Simulating via SBML")
start = pytime.clock()
result = CRN.simulate_with_bioscrape_deterministic_via_sbml(timepoints, f, initial_condition_dict)
end = pytime.clock()
print("\telapsed time=", end-start)
#result = sim_result.py_get_result()
#rep_ind = model.get_species_index("protein_Reporter")
#rna_rep_ind = model.get_species_index("rna_T1")
#rna_rep_ribo_ind = model.get_species_index("complex_rna_T1_protein_Ribo")
#rna_rep_ase_ind = model.get_species_index("complex_rna_T1_protein_RNAase")
#protein_Ribo_ind = model.get_species_index("protein_Ribo")

print('Simulating Directly')
start = pytime.clock()
result_d = CRN.simulate_with_bioscrape(timepoints, initial_condition_dict)
end = pytime.clock()
print("\telapsed time=", end-start)
#result_d = sim_result_direct.py_get_result()
#rep_ind_d = model_direct.get_species_index("protein_Reporter")
#rna_rep_ind_d = model_direct.get_species_index("rna_T1")
#rna_rep_ribo_ind_d = model_direct.get_species_index("complex_rna_T1_protein_Ribo")
#rna_rep_ase_ind_d = model_direct.get_species_index("complex_rna_T1_protein_RNAase")
#protein_Ribo_ind_d = model_direct.get_species_index("protein_Ribo")



print("Plotting")
plt.figure()
plt.subplot(121)
plt.plot(timepoints, result["rna_T1"], color = 'red', label = "rna_Reporter")
plt.plot(timepoints, result["complex_rna_T1_protein_Ribo"], color = 'blue', label = "complex_rna_T1:protein_Ribo")
plt.plot(timepoints, result["complex_rna_T1_protein_RNAase"], color ='cyan', label = "complex_rna_T1:protein_RNAase")
plt.plot(timepoints, result["protein_Ribo"], color ='pink', label = "protein_Ribo")

plt.plot(timepoints, result_d["rna_T1"], ":", color = 'darkred', label = "rna_Reporter_d")
plt.plot(timepoints, result_d["complex_rna_T1_protein_Ribo"], ":", color = 'darkblue',label = "complex_rna_T1:protein_Ribo_d")
plt.plot(timepoints, result_d["complex_rna_T1_protein_RNAase"], ":", color ='green', label = "complex_rna_T1:protein_RNAase_d")
plt.plot(timepoints, result_d["protein_Ribo"], ":", color ='purple', label = "protein_Ribo_d")

plt.legend()

plt.subplot(122)
plt.plot(timepoints, result["protein_Reporter"], color = 'green', label = "Reporter")

plt.plot(timepoints, result_d["protein_Reporter"],":", color = 'darkgreen', label = "Reporter_d")

plt.legend()
plt.show()
