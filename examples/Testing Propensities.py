from biocrnpyler.chemical_reaction_network import Species, Reaction, ComplexSpecies, ChemicalReactionNetwork

print("Start")

#Names of different supported propensities
propensity_types = ['hillpositive', 'proportionalhillpositive', 'hillnegative', 'proportionalhillnegative', 'massaction', 'general']

kb = 100
ku = 10
kex = 1.
kd = .1

G = Species(name = "G", material_type = "dna") #DNA
A = Species(name = "A", material_type = "protein") #Activator
GA = ComplexSpecies([G, A, A]) #Activated Gene
X = Species(name = "X", material_type = "protein")

rxnd = Reaction([X], [], kd)

#Massaction Unregulated
species1 = [G, A, GA, X]


rxn0_1 = Reaction([G, A, A], [GA], k=kb, k_rev = ku)
rxn0_2 = Reaction([GA], [GA, X], k=kex)
CRN0 = ChemicalReactionNetwork(species1, [rxn0_1, rxn0_2, rxnd])

rxn1_1 = Reaction([G, A, A], [GA], k=kb, k_rev = ku)
rxn1_2 = Reaction([G], [G, X], k=kex)
CRN1 = ChemicalReactionNetwork(species1, [rxn1_1, rxn1_2, rxnd])


#hill positive
species2 = [G, A, X]
rxn2_1 = Reaction([G], [G, X], propensity_type = "hillpositive", propensity_params = {"k":kex, "n":2, "K":float(kb/ku), "s1":A})
CRN2 = ChemicalReactionNetwork(species2, [rxn2_1, rxnd])

#proportional hill positive
rxn3_1 = Reaction([G], [G, X], propensity_type = "proportionalhillpositive", propensity_params = {"k":kex, "n":2, "K":float(kb/ku), "s1":A, "d":G})
CRN3 = ChemicalReactionNetwork(species2, [rxn3_1, rxnd])

#hill Negative
rxn4_1 = Reaction([G], [G, X], propensity_type = "hillnegative", propensity_params = {"k":kex, "n":2, "K":float(kb/ku), "s1":A})
CRN4 = ChemicalReactionNetwork(species2, [rxn4_1, rxnd])

#proportional hill negative
rxn5_1 = Reaction([G], [G, X], propensity_type = "proportionalhillnegative", propensity_params = {"k":kex, "n":2, "K":float(kb/ku), "s1":A, "d":G})
CRN5 = ChemicalReactionNetwork(species2, [rxn5_1, rxnd])





import numpy as np
import pylab as plt
x0 = {repr(G):2, repr(A):10}
timepoints = np.linspace(0, 100, 200)
fname = "CRN.xml"
CRNs = [CRN0, CRN1, CRN2, CRN4, CRN3, CRN5]
LSs = ["-", "--", ":"]
plt.figure()
for i in range(6):

	plt.subplot(3, 2, i+1)

	CRN = CRNs[i]
	CRN.write_sbml_file(file_name = fname)
	print("Saved")

	from bioscrape.types import Model
	from bioscrape.simulator import py_simulate_model
	from bioscrape.sbmlutil import *
	M = Model(sbml_filename = fname)

	A_list = [0, 1, 2, 5, 10]
	for ind in range(len(A_list)):
		x0[repr(A)] = A_list[ind]

		M.set_species(x0)

		R = py_simulate_model(timepoints, Model = M)

		plt.plot(timepoints, R["protein_X"], label ="A="+str(A_list[ind]),color = 'blue', alpha = (ind+1)/len(A_list))

	txt = ""
	for rxn in CRN.reactions:
		txt += repr(rxn)+"\n"
	plt.title(txt[:-1], fontsize = 8)
	plt.legend()

plt.show()
	#CRN.simulate_with_bioscrape(timepoints, initial_condition_dict = x0)
	#CRN.simulate_with_bioscrape_via_sbml(timepoints, file = fname, initial_condition_dict = x0)