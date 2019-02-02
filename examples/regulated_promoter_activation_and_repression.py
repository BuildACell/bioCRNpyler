from biocrnpyler import *
import numpy as np
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
myMixture = BasicExtract(name="txtl", parameters=parameters, components=components, parameter_warnings=False)

myCRN = myMixture.compile_crn()
print("\n"+repr(myCRN))

time = np.arange(0, 20, .01)

import pylab as plt

x0 = {"protein_activator":0, "protein_repressor":0, "dna_reporter":10, "protein_Ribo":100, "protein_RNAP":20, "protein_RNAase":10}
R_const = myCRN.simulate_with_bioscrape(time, stochastic = False, initial_condition_dict = x0)

x0 = {"protein_activator":0, "protein_repressor":50, "dna_reporter":10, "protein_Ribo":100, "protein_RNAP":20, "protein_RNAase":10}
R_repressed = myCRN.simulate_with_bioscrape(time, stochastic = False, initial_condition_dict = x0)

x0 = {"protein_activator":50, "protein_repressor":0, "dna_reporter":10, "protein_Ribo":100, "protein_RNAP":20, "protein_RNAase":10}
R_active = myCRN.simulate_with_bioscrape(time, stochastic = False, initial_condition_dict = x0)

plt.figure()
plt.plot(time, R_const["protein_reporter"], label = "Constituitive Expression")
plt.plot(time, R_repressed["protein_reporter"], label = "Repressed Expression")
plt.plot(time, R_active["protein_reporter"], label = "Activated Expression")
plt.legend()
plt.show()