from biocrnpyler import *

kb, ku, ktx, ktl, kdeg, kdil = 100, 10, 3, 2, 1, .1
parameters = {"kb":kb, "ku":ku, "ktx":ktx, "ktl":ktl, "kdeg": kdeg, "kdil":kdil}

#Creates a global mechanism that acts on all species generated except for
# those with the type or attribute "genome"
dilution_mechanism = Dilution(filter_dict = {"genome":False})

#Add this mechanism to a dictionary which is passed into the Mixture txtl.BasicExtract
global_mechanisms = {"dilution":dilution_mechanism}
myMixture = BasicExtract(name = "txtl", parameters = parameters, global_mechanisms = global_mechanisms)

#Creates a dna assembly. This assembly is type "dna" so it will be degraded
A_dna = DNAassembly(name = "G1", promoter = "pBest",
                              rbs = "BCD2")

#Create another dna assembly but set its internal specie's type to "genome" so it will not be degraded
#Note: this only protects the dna_G2 species encoded by this assembly. Protein and mRNA products will
#still be degraded by dilution. This could be overcome by creating custom transcript and protein species
#with some attribute that is passed into the filter_dict.
A_genome = DNAassembly(name = "G2", promoter = "pBest",
                              rbs = "BCD2")
A_genome.dna.type = "genome" #Note: that the code A_genome.dna.attributes.append("genome") also works here.

myMixture.add_components(A_dna)
myMixture.add_components(A_genome)
myCRN = myMixture.compile_crn()
print(repr(myCRN))