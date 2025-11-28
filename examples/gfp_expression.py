# gfp_expression.py - simple GFP expression in cells, extract, PURE
# RMM, 15 Nov 2025
#
# This file uses the various mixtures (and mechanisms) to express GFP and
# compare results across cells, extract, and PURE.  It is intended to allow
# checking whether default parameter values do something reasonable and how
# the different mixtures and mechanisms compare.

import matplotlib.pyplot as plt
import numpy as np

import biocrnpyler as bcp

# Set up some "units" for use later
nM, uM, mM = 1e-3, 1, 1e3  # default units for concentrations
sec, min, hrs = 1, 60, 3600  # default units for time

# Create a DNA assembly for strong expression of GFP
gfp_dna = bcp.DNAassembly(
    name='gfp', promoter='pconst', rbs='rbs_strong', protein='GFP'
)

# Simulation parmaters
initial_conditions_dict = {'dna_gfp': 1 * nM}
timepts = np.linspace(0, 6 * hrs, 1000)

#
# Extract-based, one step expression
#
expr_mixture = bcp.ExpressionExtract(
    name='expression',
    components=[gfp_dna],
    parameter_file='mixtures/extract_parameters.tsv',
)
expr_crn = expr_mixture.compile_crn()
expr_res = expr_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions_dict
)

#
# Extract-based, simple expression
#
simple_mixture = bcp.SimpleTxTlExtract(
    name='simple',
    components=[gfp_dna],
    parameter_file='mixtures/extract_parameters.tsv',
)
simple_crn = simple_mixture.compile_crn()
simple_res = simple_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions_dict
)

#
# Extract-based, with machinery
#
regular_mixture = bcp.TxTlExtract(
    name='regular',
    components=[gfp_dna],
    parameter_file=[
        'mixtures/extract_parameters.tsv',
    ],
)
regular_crn = regular_mixture.compile_crn()
regular_res = regular_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions_dict
)

#
# Extract-based, with energy
#
energy_mixture = bcp.EnergyTxTlExtract(
    name='energy',
    components=[gfp_dna],
    parameter_file=[
        'mixtures/extract_parameters.tsv',
    ],
)
energy_crn = energy_mixture.compile_crn()
energy_res = energy_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions_dict
)

#
# Comparison of extract-based expression mixtures
#

plt.figure(1)
plt.clf()
plt.plot(timepts / min, expr_res['protein_GFP'] / uM, 'k', label='GFP, expr')
plt.plot(timepts / min, simple_res['protein_GFP'] / uM, label='GFP, simple')
plt.plot(timepts / min, regular_res['protein_GFP'] / uM, label='GFP, regular')
plt.plot(timepts / min, energy_res['protein_GFP'] / uM, label='GFP, energy')

plt.title("Extract Mixture Comparisions - Protein")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend()

#
# mRNA comparisons
#
plt.figure(2)
plt.clf()
plt.plot(timepts / min, simple_res['rna_gfp'] / uM, label='mRNA, simple')
plt.plot(
    timepts/min,
    regular_res['complex_protein_Ribo_rna_gfp_'],
    label='mRNA:Ribo, reg')
plt.plot(
    timepts / min,
    energy_res['complex_metabolite_ATP_4x_metabolite_amino_acids_protein_Ribo_rna_gfp_'] / uM,
    label='mRNA, energy')

plt.title("Extract Mixture Comparisions - RNA")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend()

#
# Analysis of energy-based mixture
#

plt.figure(3)
plt.clf()
plt.plot(timepts / min, energy_res['metabolite_ATP'] / mM, 'b-', label='ATP')
plt.plot(timepts / min, energy_res['metabolite_ADP'] / mM, 'b:', label='ADP')
plt.plot(timepts / min, energy_res['metabolite_NTPs'] / mM, 'r', label='NTPs')
plt.plot(
    timepts / min,
    energy_res['metabolite_amino_acids'] / mM / 10,
    'k',
    label='AAs/10',
)
plt.plot(
    timepts / min,
    energy_res['metabolite_Fuel_3PGA'] / 10 / mM,
    label='3PGA/10',
)
# plt.plot(timepts/min, energy_res['metabolite_Fuel_3PGA']/mM, label='3PGA')
plt.plot(timepts / min, energy_res['protein_GFP'] / mM, 'g', label='GFP')

plt.title("Resource Utilization: EnergyTxTlExtract")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [mM]")
plt.legend()
