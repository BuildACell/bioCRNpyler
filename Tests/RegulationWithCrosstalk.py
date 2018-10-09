import txtl
import component
import dna_assembly
import mechanism
import CRN_Simulator
import numpy as np
#Parameters
kb, ku, ktx, ktl, kdeg = 500, 10, 2.0, .25, 1.5
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg":kdeg,
              "cooperativity":2, "ktx_leak":ktx, "kb_leak":kb/50, "ku_leak":ku*50}

N = 5
promoters_names = ["P"+str(i) for i in range(N)]

exp_range = 3
cross_talk_exp_matrix = np.random.randn(size = (N, N))
for i in range(N):
    pass




P_reg = dna_assembly.RegulatedPromoter("P_regulated", regulators=["activator", "repressor"], leak = True)

reg_rep_assembly = dna_assembly.DNAassembly(name = "reporter", promoter = P_reg, rbs = "BCD")

activator = component.Protein("activator")
repressor = component.Protein("repressor")

components = [reg_rep_assembly, activator, repressor]
myMixture = txtl.BasicExtract(name = "txtl", parameters = parameters, components = components)

myCRN = myMixture.compile_crn()