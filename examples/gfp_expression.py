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
from biocrnpyler.utils.units import nM, uM, mM, sec, min, hrs

# Create a DNA assembly for strong expression of GFP
gfp_dna = bcp.DNAassembly(
    name='gfp', promoter='pconst', rbs='rbs_strong', protein='GFP'
)

# Simulation parameters
initial_conditions = {'dna_gfp': 1 * nM}
timepts = np.linspace(0, 6 * hrs, 1000)

#
# Extract-based, one step expression
#
expr_mixture = bcp.ExpressionExtract(
    name='expression',
    components=[gfp_dna],
)
expr_crn = expr_mixture.compile_crn()
expr_res = expr_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions
)

#
# Extract-based, simple expression
#
simple_mixture = bcp.SimpleTxTlExtract(
    name='simple',
    components=[gfp_dna],
)
simple_crn = simple_mixture.compile_crn()
simple_res = simple_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions
)

#
# Extract-based, with machinery
#
regular_mixture = bcp.TxTlExtract(
    name='regular',
    components=[gfp_dna],
)
regular_crn = regular_mixture.compile_crn()
regular_res = regular_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions
)

#
# Extract-based, with energy
#
energy_mixture = bcp.EnergyTxTlExtract(
    name='energy',
    components=[gfp_dna],
)
energy_crn = energy_mixture.compile_crn()
energy_res = energy_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions
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
    timepts / min,
    regular_res['complex_protein_Ribo_rna_gfp_'],
    label='mRNA:Ribo, reg',
)
plt.plot(
    timepts / min,
    energy_res[
        'complex_metabolite_ATP_4x_metabolite_amino_acids_protein_Ribo_rna_gfp_'
    ] / uM,
    label='mRNA, energy',
)

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
    energy_res['metabolite_Fuel_3PGA'] / mM / 10,
    label='3PGA/10',
)
# plt.plot(timepts/min, energy_res['metabolite_Fuel_3PGA']/mM, label='3PGA')
plt.plot(timepts / min, energy_res['protein_GFP'] / mM, 'g', label='GFP')

plt.title("Resource Utilization: EnergyTxTlExtract")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [mM]")
plt.legend()

#
# Effect of a competing gene on the extract mixtures
#
# Figures 1-3 use 1 nM of DNA, where none of the machinery pools in
# TxTlExtract is close to limiting.  mRNA settles at roughly ktx/kdil = 60
# times the DNA concentration, so 1 nM of DNA yields only tens of nM of
# transcript against 1.1 uM of ribosomes, and a second gene changes
# nothing.  Raising both genes to 10 nM brings the transcript pool within
# reach of the ribosomes and separates the three mixtures:
#
#   simple   models no machinery at all, so the load is invisible
#   regular  shares ribosomes between the genes, so GFP drops
#   energy   has a fixed metabolite pool, so the two genes split it evenly
#
# The absolute yields are worth as much attention as the ratios.  Without
# an energy model the extracts run to hundreds of uM, far past what a real
# TX-TL reaction reaches, while the energy mixture saturates at a
# realistic level.  The vertical axis is truncated so that the energy
# traces stay readable, so final values are quoted in the legend.
#
plt.figure(4)
plt.clf()

load_conc = 10 * nM
load_conditions = {'dna_gfp': load_conc}
load_cfp_conditions = {'dna_gfp': load_conc, 'dna_cfp': load_conc}

# Add some additional DNA that will utilize resources
cfp_dna = bcp.DNAassembly(
    name='cfp', promoter='pconst', rbs='rbs_strong', protein='CFP'
)

# Re-run the single gene mixtures at the higher DNA concentration, so that
# the only difference within each pair below is the presence of CFP
simple_load_res = simple_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_conditions
)
regular_load_res = regular_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_conditions
)
energy_load_res = energy_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_conditions
)

# Simple mixture should not be affected
cfp_simple_mixture = bcp.SimpleTxTlExtract(
    name='energy',
    components=[gfp_dna, cfp_dna],
    parameter_file=[
        'mixtures/extract_parameters.tsv',
    ],
)
cfp_simple_crn = cfp_simple_mixture.compile_crn()
cfp_simple_res = cfp_simple_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_cfp_conditions
)

# Regular mixture shares RNAP, ribosomes and RNase across both genes, and
# at 10 nM the ribosomes are loaded enough for that to cost GFP
cfp_regular_mixture = bcp.TxTlExtract(
    name='energy',
    components=[gfp_dna, cfp_dna],
    parameter_file=[
        'mixtures/extract_parameters.tsv',
    ],
)
cfp_regular_crn = cfp_regular_mixture.compile_crn()
cfp_regular_res = cfp_regular_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_cfp_conditions
)

# Energy mixture caps on its metabolite pool rather than on DNA, so two
# genes of equal strength split that pool and GFP halves
cfp_energy_mixture = bcp.EnergyTxTlExtract(
    name='energy',
    components=[gfp_dna, cfp_dna],
    parameter_file=[
        'mixtures/extract_parameters.tsv',
    ],
)
cfp_energy_crn = cfp_energy_mixture.compile_crn()
cfp_energy_res = cfp_energy_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=load_cfp_conditions
)


def plot_cfp_pair(alone_res, cfp_res, label):
    """Plot GFP with and without CFP, quoting the final value of each."""
    alone = alone_res['protein_GFP'] / uM
    loaded = cfp_res['protein_GFP'] / uM
    final_alone = alone.to_numpy()[-1]
    final_loaded = loaded.to_numpy()[-1]
    lines = plt.plot(
        timepts / min,
        alone,
        '--',
        label=f'GFP, {label} ({final_alone:.0f} uM)',
    )
    plt.plot(
        timepts / min,
        loaded,
        color=lines[0].get_color(),
        label=(
            f'GFP, {label} w/ CFP ({final_loaded:.0f} uM, '
            f'{final_loaded / final_alone:.2f}x)'
        ),
    )


plot_cfp_pair(simple_load_res, cfp_simple_res, 'simple')
plot_cfp_pair(regular_load_res, cfp_regular_res, 'regular')
plot_cfp_pair(energy_load_res, cfp_energy_res, 'energy')

plt.ylim(0, 200)
plt.title("Extract Mixture Comparisions w/ CFP (10 nM each)")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend(fontsize=8)

#
# Comparison with PURE
#

pure_mixture = bcp.PURE(
    name='regular',
    components=[gfp_dna],
    include_machinery=True,
    include_resources=True,
    include_energy=True,
)
pure_crn = pure_mixture.compile_crn()
pure_res = pure_crn.simulate_with_bioscrape_via_sbml(
    timepts, initial_condition_dict=initial_conditions
)

plt.figure(5)
plt.clf()
plt.plot(timepts / min, energy_res['protein_GFP'] / uM, label='GFP, TX-TL')
plt.plot(timepts / min, pure_res['protein_GFP'] / uM, label='GFP, PURE')

plt.title("Mixture Comparisions - TX-TL vs PURE")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend()

#
# PURE mixtures
#

pure_timepts = np.linspace(0, 180 * min)

# The PURE comparisons below use 1 nM of each gene
cfp_initial_conditions = {'dna_gfp': 1 * nM, 'dna_cfp': 1 * nM}

simple_mixture = bcp.PURE(
    name='simple', components=[gfp_dna, cfp_dna],
    include_machinery=False, include_resources=False, include_energy=False)
simple_crn = simple_mixture.compile_crn()
simple_res = simple_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=cfp_initial_conditions
)

machinery_mixture = bcp.PURE(
    name='machinery', components=[gfp_dna, cfp_dna],
    include_machinery=True, include_resources=False, include_energy=False)
machinery_crn = machinery_mixture.compile_crn()
machinery_res = machinery_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=cfp_initial_conditions
)

resource_mixture = bcp.PURE(
    name='resources', components=[gfp_dna, cfp_dna],
    include_machinery=True, include_resources=True, include_energy=False)
resource_crn = resource_mixture.compile_crn()
resource_res = resource_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=cfp_initial_conditions
)

energy_mixture = bcp.PURE(
    name='energy', components=[gfp_dna, cfp_dna],
    include_machinery=True, include_resources=True, include_energy=True)
energy_crn = energy_mixture.compile_crn()
energy_res = energy_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=cfp_initial_conditions
)

energy_mixture_gfp = bcp.PURE(
    name='energy', components=[gfp_dna],
    include_machinery=True, include_resources=True, include_energy=True)
energy_crn_gfp = energy_mixture_gfp.compile_crn()
energy_res_gfp = energy_crn_gfp.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=initial_conditions
)

plt.figure(6)
plt.clf()

plt.plot(
    pure_timepts / min, simple_res['protein_GFP'] / uM,
    label='GFP+CFP, simple')
plt.plot(
    pure_timepts / min, machinery_res['protein_GFP'] / uM,
    label='GFP+CFP, machinery')
plt.plot(
     pure_timepts / min, resource_res['protein_GFP'] / uM,
     label='GFP+CFP, resources')
plt.plot(
    pure_timepts / min, energy_res['protein_GFP'] / uM,
    label='GFP+CFP, energy')
plt.plot(
    pure_timepts / min, energy_res_gfp['protein_GFP'] / uM,
    label='GFP only, energy')
plt.plot(
    timepts[0:50] / min, pure_res['protein_GFP'][0:50] + 0.1 / uM, '--',
    label='GFP only + 0.1, basic')

plt.ylim([0, 12])

plt.title("PURE Mixture Comparisions")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend()

#
# Activators and repressors
#

TetR = bcp.Protein('TetR')
aTc = bcp.Species('aTc')
TetR_inactive = bcp.ChemicalComplex(
    [TetR.species, aTc], mechanisms={'binding': bcp.One_Step_Binding()}
)
ptet_repressed = bcp.RepressiblePromoter('ptet', TetR)
dna_GFP_repressed = bcp.DNAassembly(
    'GFP', promoter=ptet_repressed, rbs='RBS', protein='GFP', length=714
)

initial_conditions = {'dna_GFP': 1 * nM, 'protein_TetR': 30 * uM}
repressed_mixture = bcp.PURE(
    name='repressed',
    components=[dna_GFP_repressed, TetR_inactive],
    include_machinery=True,
    include_resources=True,
    include_energy=True,
)
repressed_crn = repressed_mixture.compile_crn()

TetR = bcp.Protein('TetR')
aTc = bcp.Species('aTc')
TetR_inactive = bcp.ChemicalComplex(
    [TetR.species, aTc], mechanisms={'binding': bcp.One_Step_Binding()}
)
ptet_regulated = bcp.RegulatedPromoter('ptet', TetR)
dna_GFP_regulated = bcp.DNAassembly(
    'GFP', promoter=ptet_regulated, rbs='RBS', protein='GFP'
)
regulated_parameters = {
    ('transcription', 'ptet_leak', 'ktx'): 50,          # unbound
    ('transcription', 'ptet_TetR', 'ktx'): 0.001,       # bound
}

regulated_mixture = bcp.PURE(
    name='regulated',
    components=[dna_GFP_regulated, TetR_inactive],
    include_machinery=True,
    include_resources=True,
    include_energy=True,
    parameters=regulated_parameters,
)
regulated_crn = regulated_mixture.compile_crn()

plt.figure(7)
plt.clf()

offset = -0.01
TetR_0 = initial_conditions['protein_TetR']
for aTc_0 in np.linspace(0, 50*uM, 5):
# aTc_0 = 0
# for TetR_0 in np.linspace(0, 20*uM, 5):
    repressed_res = repressed_crn.simulate_with_bioscrape_via_sbml(
        pure_timepts, initial_condition_dict=initial_conditions
        | {'aTc': aTc_0} | {'protein_TetR': TetR_0}
    )
    lines = plt.plot(
        pure_timepts / min,
        repressed_res['protein_GFP'] / uM + offset,
        label=f"aTc={aTc_0 / uM} uM, TetR={TetR_0 / uM} uM",
    )

    regulated_res = regulated_crn.simulate_with_bioscrape_via_sbml(
        pure_timepts, initial_condition_dict=initial_conditions
        | {'aTc': aTc_0} | {'protein_TetR': TetR_0}
    )
    plt.plot(
        pure_timepts / min,
        regulated_res['protein_GFP'] / uM + offset,
        color=lines[0].get_color(),
        linestyle='--',
    )

    offset += 0.05

plt.title("Represssed (-) versus Regulated (--)")
plt.xlabel("Time [min]")
plt.ylabel("Concentration [uM]")
plt.legend()

#
# PURE debugging
#

repressed_res = repressed_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=initial_conditions
    | {'aTc': 37.5 * uM} | {'protein_TetR': 30 * uM}
)

plt.figure(8)
plt.clf()
bcp.plot_gene_expression_data(
    repressed_res, repressed_crn, dna_GFP_repressed,
    trace_offset=[0.01, 0.002, 0.1, 0.1])
plt.suptitle("Gene expression, repressed", fontsize='large')
plt.tight_layout()

regulated_res = regulated_crn.simulate_with_bioscrape_via_sbml(
    pure_timepts, initial_condition_dict=initial_conditions
    | {'aTc': 37.5 * uM} | {'protein_TetR': 30 * uM}
)

plt.figure(9)
plt.clf()
bcp.plot_gene_expression_data(
    regulated_res, regulated_crn, dna_GFP_regulated,
    trace_offset=[0.01, 0.002, 0.1, 0.1])
plt.suptitle("Gene expression, regulated", fontsize='large')
plt.tight_layout()

#
# Species distributed across complexes
#
# Show how a species is distributed between its free form and the complexes
# that contain it.  The ribosome is split between free ribosome and the
# translation complexes; TetR is sequestered by aTc.
#

plt.figure(10)
plt.clf()

ribosome = bcp.Species('Ribo', material_type='protein')

plt.subplot(1, 2, 1)
bcp.plot_all_species_containing(
    repressed_res, repressed_crn, ribosome,
    show_material=False, show_complex_material=False,
    legend_fontsize='x-small'
)
plt.title("Ribosome", fontsize='small')

plt.subplot(1, 2, 2)
bcp.plot_all_species_containing(
    repressed_res, repressed_crn, [TetR, aTc],
    show_material=False, show_complex_material=False,
    legend_fontsize='x-small'
)
plt.title("TetR and aTc", fontsize='small')

plt.suptitle("Species containing a given species", fontsize='large')
plt.tight_layout()
