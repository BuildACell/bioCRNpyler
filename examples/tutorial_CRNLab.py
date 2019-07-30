# Tube level interface of biocrnpyler for users
from biocrnpyler import *
import numpy as np

# Set up the standard TXTL tubes by creating a Mixture, give any name you'd like for the mixture
# from .crnlab import CRNLab

txtl = CRNLab("txtl_gfp")

txtl.mixture("mixture1", extract = "test_extract", extract_volume = 1e-6)

# OR, Optionally, do the following, 
# initial_concentration_dict = {"protein_Ribo":10, "protein_RNAP":5, "protein_RNAase":2.5} # If not, default values will be used
# Specify what extract and buffer to use
# txtl.extract("test_extract", init = initial_condition_dict)
# Buffer with energy models will be available in future releases
# txtl.buffer("stdbuffer")

# Create genes to add to the mix as follows
# gene1 = DNAassembly(name = "G1", promoter = "pBest", rbs = "BCD2", protein = "GFP", final_conc = 10)
# txtl.add_dna(gene1)

# Or simply add a new dna in one line as follows.
txtl.add_dna(name = "G2", promoter = "pBest", rbs = "BCD2", protein = "tetR", initial_conc = 10, volume = 1e-7)

# Combine all of the tubes together to get the model
well1 = txtl.get_model()
print(well1)
# Create an SBML file containing the model
filename = "geneexpr.xml"
txtl.write_sbml_file(filename)

# Run a simulation (using bioscrape) and plot the result
# (Optional) Specify the type of simulation (deterministic or stochastic)
timepoints = np.linspace(0,14*60*60,100)
simdata = well1.runsim_bioscrape(timepoints, filename, simtype = "stochastic", 
                                # plot_show = True)
                                species_to_plot= ["protein_tetR"], plot_show = True)

# OR, you could use RoadRunner for simulation of the SBML model, simply call
# simdata = well1.runsim_roadrunner(timepoints, filename)
# OR, just use any other simulator of your choice with the SBML model stored in the file above.

 