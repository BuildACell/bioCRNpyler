from biocrnpyler import *

#Create a simple CRN
s1 = Species(name="G",  material_type="DNA") #A species "G" of type "DNA"
s1double = Species("G",  material_type ="DNA") #G double species of the same name should be automatically removed

s2 = Species(name="R",  material_type="protein") #A species "R" of type "DNA"
s3 = Species(name="R",  material_type ="dimer") #Two species with the same name but different types are formally different
s4 = Species(name="R",  material_type="protein", attributes = ["deg"]) #Two species with the same name and type but different attributes are formally different

#Create Reactions
#DNA_G --> DNA_G + protein_R @ k=1
r1 = Reaction([s1], [s1, s2], k=1.0)

#2 protein_R <--> dimer_R @ kf = 100, krev = 10
r2 = Reaction([s2, s2], [s3], k=100., k_rev=10.) #Note duplicates in inputs

# dimer_R <--> 2 protein_R_deg
r3 = Reaction([s3], [s4], k=100., k_rev = 10., output_coefs=[2]) #Note use of output_coefs instead of duplicates in outputs

#protein_R_deg --> NULL
r4 = Reaction([s4], [], k = .1) #Empty set denoted as empty brackets
species = [s1, s1double, s2, s3, s4]
rxns = [r1, r2, r3, r4]
CRN = ChemicalReactionNetwork(species, rxns)
print(repr(CRN))

print("\nWriting SBML")
file_name = "sbml_test_file.xml"
f = CRN.write_sbml_file(file_name)

