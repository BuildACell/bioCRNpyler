from biocrnpyler import *

#Create Param Dict
kb, ku, ktx, ktl, kdeg, cooperativity = 100, 10, 3, 2, 1, 1
kb_dcas_grna, ku_dcas_rgna = 100, .01 #binding constants for dCas9 and Guide RNAs
kb_dcas_dna, ku_dcas_dna = 100, .1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg, "cooperativity":cooperativity, #default params
              ("dCas9_binding", "kb"):kb_dcas_grna, ("dCas9_binding", "ku"):ku_dcas_rgna, #binding constants for gRNA and dCas9
              ("dCas9_dna_binding", kb): kb_dcas_dna, ("dCas9_dna_binding", ku): ku_dcas_dna,
              "kb_leak": kb, "ku_leak": ku, "ktx_leak": ktx}#binding constants for DNA and dcas-guide complex

#Create an assembly to express dCas9
const_dCas_assmebly = DNAassembly("dCas9", promoter = "P", rbs = "BCD")
#Get the dCas Species. Note: This could also be defined above and passed into the assmebly as protein=...
dCas = const_dCas_assmebly.protein
#Create a Guide RNA Guide1
gRNA = guideRNA("guide1", dCas=dCas)
#Create an assembly to express the guide RNA
const_gRNA_assembly = DNAassembly("gRNA", transcript = gRNA, promoter = "P", rbs = None)

#Get the guideRNA:dCas9 complex
repressor = gRNA.get_dCasComplex()
#set the unbinding  rate of the repressor dna complex with rnap to be very high (preventing transcription)
parameters[("transcription", repressor.name, "ku")] = 10000
#set the RNAP transcription rate lower when gRNA:dCAS9 is bound
parameters[("transcription", repressor.name, "ktx")] = .1

#Create a Promoter regulated by the repressor
P_reg = RegulatedPromoter("P_regulated", regulators=[repressor], leak = True)
P_reg.default_mechanisms['binding'] = Reversible_Bimolecular_Binding("dCas9_dna_binding")
#Create an assembly with the regulated promoter
reg_assembly = DNAassembly(name = "reporter", promoter = P_reg, rbs = "BCD")

#Create a list of components to add to the mixture (these could also be added one-by-one with Mixture.add_component(...)
components = [const_dCas_assmebly, const_gRNA_assembly, gRNA, reg_assembly]

#Create a BasicExtract Mixture
reaction_mix = BasicExtract("txtl", components = components, parameters=parameters)

#Compile a CRN
CRN = reaction_mix.compile_crn()
print(repr(CRN))


print("Simulating with BioSCRAPE")
#Initial Condition Dict: repr(specie) --> concentration. Default is 0
#Note in SBML all the species names have ":" replaced with "_".
x0_no_dcas = {repr(const_dCas_assmebly.dna):0, repr(reg_assembly.dna):1, repr(const_gRNA_assembly.dna):10, "protein_Ribo":20, "protein_RNAP":10, "protein_RNAase":10}
x0_with_dcas = {repr(const_dCas_assmebly.dna):2., repr(reg_assembly.dna):1, repr(const_gRNA_assembly.dna):10, "protein_Ribo":20, "protein_RNAP":10, "protein_RNAase":10}

#Bioscrape Simulation
import numpy as np
import pylab as plt
timepoints = np.arange(0, 10, .01)

sim_no_cas = CRN.simulate_with_bioscrape(timepoints, x0_no_dcas, stochastic = False)
sim_with_cas = CRN.simulate_with_bioscrape(timepoints, x0_with_dcas, stochastic = False)


plt.figure(figsize = (10, 10))
plt.subplot(311)
plt.title("Protein Products")
plt.plot(timepoints, sim_no_cas["protein_reporter"], color = "red", label = "Reporter: No dCas9")
plt.plot(timepoints, sim_with_cas["protein_reporter"],":", color = "red", label = "Reporter: With dCas9")
plt.plot(timepoints, sim_no_cas[repr(dCas)], color = "blue", label = "dCas9: No dCas9")
plt.plot(timepoints, sim_with_cas[repr(dCas)],":", color = "blue", label = "dCas9: With dCas9")
plt.legend()

plt.subplot(312)
plt.title("RNAs")
plt.plot(timepoints, sim_no_cas[repr(gRNA.get_specie())], color = "cyan", label = "gRNA: No dCas9")
plt.plot(timepoints, sim_with_cas[repr(gRNA.get_specie())],":", color = "cyan", label = "gRNA: With dCas9")
plt.legend()

plt.subplot(313)
plt.plot(timepoints, sim_no_cas[repr(reg_assembly.dna)], color = "red", label = "dna_reporter: No dCas9")
plt.plot(timepoints, sim_with_cas[repr(reg_assembly.dna)],":", color = "red", label = "dna_reporter: With dCas9")
plt.plot(timepoints, sim_no_cas["complex_rna_guide1_protein_dCas9_dna_reporter"], color = "blue", label = "repressed dna_reporter: No dCas9")
plt.plot(timepoints, sim_with_cas["complex_rna_guide1_protein_dCas9_dna_reporter"],":", color = "blue", label = "repressed dna_reporter: With dCas9")

plt.plot()
plt.legend()
plt.show()

