import txtl
import dna_assembly

kb, ku, ktx, ktl, kdeg = 100, 10, 3, 2, 1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg}
myMixture = txtl.BasicExtract(name = "txtl", parameters = parameters)

A1 = dna_assembly.DNAassembly(name = "G1", promoter = "pBest",
                              rbs = "BCD2", transcript = "T1", protein = "GFP")

#Note: Protein and Transcript strings (or chemical_reaction_network.specie objects) are optional parameters
#DNAassemblies default to using their name for their transcript and protein products.

myMixture.add_components(A1)
myCRN = myMixture.compile_crn()

print("\n"+repr(A1))
print("\n"+repr(myMixture))
print("\n"+repr(myCRN))









#print("\nmyMixture.parameters", myMixture.parameters)
#print("\ndna_assembly.parameters", A1.parameters)
