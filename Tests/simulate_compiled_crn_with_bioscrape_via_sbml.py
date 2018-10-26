import bioscrape
import numpy as np
import pylab as plt
import txtl
import dna_assembly

kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg}
myMixture = txtl.BasicExtract(name = "txtl", parameters = parameters)

A1 = dna_assembly.DNAassembly(name = "G1", promoter = "P",
                              rbs = "BCD", transcript = "T1", protein = "Reporter")

myMixture.add_components(A1)
CRN = myMixture.compile_crn()
print("\n"+repr(CRN))

file_name = "bioscrape_test.xml"
f = CRN.write_sbml_file(file_name)

#Initial Condition Dict: repr(specie) --> concentration. Default is 0
#Note in SBML all the species names have ":" replaced with "_".
initial_condition_dict = {repr(A1.dna):2, "complex_Ribo":10, "protein_RNAP":5, "protein_RNAase":2.5}

#Bioscrape Simulation
timepoints = np.arange(0, 50, .01)
print("Simulating")
sim_result, model= CRN.simulate_with_bioscrape_deterministic(timepoints, f, initial_condition_dict)

result = sim_result.py_get_result()
rep_ind = model.get_species_index("protein_Reporter")
rna_rep_ind = model.get_species_index("rna_T1")
rna_rep_ribo_ind = model.get_species_index("complex_rna_T1_complex_Ribo")
rna_rep_ase_ind = model.get_species_index("complex_rna_T1_protein_RNAase")
complex_ribo_ind = model.get_species_index("complex_Ribo")

print("Plotting")
plt.figure()
plt.subplot(121)
plt.plot(timepoints, result[:, rna_rep_ind], label = "rna_Reporter")
plt.plot(timepoints, result[:, rna_rep_ribo_ind], label = "complex_rna_T1:complex_Ribo")
plt.plot(timepoints, result[:, rna_rep_ase_ind], label = "complex_rna_T1:protein_RNAase")
plt.plot(timepoints, result[:, complex_ribo_ind], label = "complex_Ribo")


print("result[-1, rna_rep_ind]", result[-1, rna_rep_ind])
print("result[-1, rna_rep_ribo_ind]", result[-1, rna_rep_ribo_ind])
print("result[-1, rna_rep_ase_ind]", result[-1, rna_rep_ase_ind])
print("result[-1, complex_ribo_ind]", result[-1, complex_ribo_ind])

plt.legend()

plt.subplot(122)
plt.plot(timepoints, result[:, rep_ind], label = "Reporter")
plt.legend()
plt.show()
