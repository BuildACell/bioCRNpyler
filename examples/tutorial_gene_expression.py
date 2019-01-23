# Tube level interface of biocrnpyler for users
from biocrnpyler import *
import numpy as np

# Set up the standard TXTL tubes by creating a Mixture, give any name you'd like for the mixture
mix = txtl(Mixture('txtl_gfp'))

# Specify what extract and buffer to use
tube1 = mix.extract('test_extract')

# Buffer with energy models will be available in future releases
# tube2 = mix.buffer('stdbuffer')

# Create genes to add to the mix as follows
gene1 = DNAassembly(name = "G1", promoter = 'pBest', rbs = 'BCD2', protein = "GFP", initial_conc = 10)
mix.add_dna(gene1)

gene2 = DNAassembly(name = "G2", promoter = 'pBest', rbs = 'BCD2', protein = "tetR", initial_conc = 10)
mix.add_dna(gene2)

# Combine all of the tubes together to get the model
well1 = mix.combine_tubes()

# Create an SBML file containing the model
filename = 'geneexpr.xml'
well1.write_sbml_file(filename)

# Run a simulation (using bioscrape) and plot the result
# (Optional) Specify the type of simulation (deterministic or stochastic)
timepoints = np.linspace(0,14,100)
simdata = well1.runsim_bioscrape(timepoints, filename, simtype = "stochastic")

# OR, you could use RoadRunner for simulation of the SBML model, simply call
# simdata = well1.runsim_roadrunner(timepoints, filename)
# OR, just use any other simulator of your choice with the SBML model stored in the file above.