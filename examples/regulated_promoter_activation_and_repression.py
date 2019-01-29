from biocrnpyler import *
import numpy as np
import pylab as plt
#Parameters
kb, ku, ktx, ktl, kdeg = 200, 10, 2.0, 50.0, 1.5
#(mechanism.name, part_id, param_name)
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg":kdeg,
              "cooperativity":4,
              ('translation_mm', 'BCD', 'ku'):ku, ('translation_mm', 'BCD', 'kb'): kb, ('translation_mm', 'BCD', 'ktl'):ktl,
              ('transcription_mm', 'activator', "ktx"): ktx, ('transcription_mm', 'repressor', "ktx"): ktx,
              ('one_step_cooperative_binding', 'repressor', 'kb'):1000, ('one_step_cooperative_binding',"repressor", 'ku'):5.0,
              ('transcription_mm', 'repressor', 'kb'):1, ('transcription_mm',"repressor", 'ku'):1000.0,
              ('one_step_cooperative_binding', 'activator', 'kb'):1000, ('one_step_cooperative_binding', "activator", 'ku'):5.0,
              ('transcription_mm', 'activator', 'kb'): 1000, ('transcription_mm', "activator", 'ku'): 1.0,
              ('transcription_mm', 'P_regulated', "kb_leak"): kb/10,('transcription_mm', 'P_regulated', "ku_leak"): ku*10,
              ('transcription_mm', 'P_regulated', "ktx_leak"):ktx}


P_reg = RegulatedPromoter("P_regulated", regulators=["activator", "repressor"], leak=True)

reg_rep_assembly = DNAassembly(name="reporter", promoter=P_reg, rbs="BCD")

activator = Protein("activator")
repressor = Protein("repressor")

components = [reg_rep_assembly, activator, repressor]
myMixture = BasicExtract(name="txtl", parameters=parameters, components=components)

myCRN = myMixture.compile_crn()
print("\n"+repr(myCRN))

timepoints = np.arange(0, 20, .01)


print("\n" + repr(myCRN))

# TODO Convert simulation code to bioscrape
"""
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
>>>>>>> 77dd0983a072372ff43f89d0e3ec4f9f9d252e46:examples/regulated_promoter_activation_and_repression.py

#Simulator without repressor or activator
x0_dict = {"protein_RNAP":10., "protein_RNAase":50.0, "protein_Ribo":1000.,
               'dna_reporter':5., 'protein_activator':0, 'protein_repressor':0}
results = myCRN.simulate_with_bioscrape(timepoints, x0_dict, stochastic = False)

#Simulator with repressor
x0_dict = {"protein_RNAP":10., "protein_RNAase":50.0, "protein_Ribo":1000.,
               'dna_reporter':5., 'protein_activator':0, 'protein_repressor':20}
results_rep = myCRN.simulate_with_bioscrape(timepoints, x0_dict, stochastic = False)

x0_dict = {"protein_RNAP":10., "protein_RNAase":50.0, "protein_Ribo":1000.,
               'dna_reporter':5., 'protein_activator':20, 'protein_repressor':0}
results_act = myCRN.simulate_with_bioscrape(timepoints, x0_dict, stochastic = False)

rna_tot = np.sum(results[myCRN.get_all_species_containing(reg_rep_assembly.transcript, return_as_strings=True)], 1)
rna_rep_tot = np.sum(results_rep[myCRN.get_all_species_containing(reg_rep_assembly.transcript, return_as_strings=True)], 1)
rna_act_tot = np.sum(results_act[myCRN.get_all_species_containing(reg_rep_assembly.transcript, return_as_strings=True)], 1)
protein_tot = np.sum(results[myCRN.get_all_species_containing(reg_rep_assembly.protein, return_as_strings=True)], 1)
protein_rep_tot = np.sum(results_rep[myCRN.get_all_species_containing(reg_rep_assembly.protein, return_as_strings=True)], 1)
protein_act_tot = np.sum(results_act[myCRN.get_all_species_containing(reg_rep_assembly.protein, return_as_strings=True)], 1)



plt.subplot(211)
plt.title("Deterministic Simulation of an Activatable/Repressible DNA Assembly")
plt.ylabel("reporter mRNA")
plt.xlabel("time")
plt.plot(timepoints, rna_tot, label = "no repressor or activator")
plt.plot(timepoints, rna_rep_tot, label = "repressor present")
plt.plot(timepoints, rna_act_tot, label = "activator present")
plt.legend()

plt.subplot(212)
plt.ylabel("reporter protein")
plt.xlabel("time")
plt.plot(timepoints, protein_tot, label = "no repressor or activator")
plt.plot(timepoints, protein_rep_tot, label = "repressor present")
plt.plot(timepoints, protein_act_tot, label = "activator present")
plt.legend()
<<<<<<< HEAD:Tests/regulated_promoter_activation_and_repression.py

plt.show()
=======
plt.show()"""