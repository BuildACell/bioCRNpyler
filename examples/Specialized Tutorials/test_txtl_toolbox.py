import matplotlib.pyplot as plt
import numpy as np

from biocrnpyler.components import DNAassembly
from biocrnpyler.mixtures import EnergyTxTlExtract

A = DNAassembly(
    'A', promoter='P', rbs='rbs', initial_concentration=1 * 10**-6
)
E = EnergyTxTlExtract(
    components=[A], parameter_file='txtl_toolbox_parameters.txt'
)

CRN = E.compile_crn()

print(CRN.pretty_print())
maxtime = 30000
timepoints = np.arange(0, maxtime, 100)
R = CRN.simulate_with_bioscrape_via_sbml(timepoints)
plt.subplot(121)
plt.plot(timepoints, R[str(E.ntps.get_species())], label=E.ntps.get_species())
plt.plot(
    timepoints,
    R[str(E.amino_acids.get_species())],
    label=E.amino_acids.get_species(),
)
plt.plot(timepoints, R[str(E.fuel.get_species())], label=E.fuel.get_species())
plt.xticks(
    np.arange(0, maxtime, 3600),
    [str(i) for i in range(0, int(np.ceil(maxtime / 3600)))],
)
plt.legend()

plt.subplot(122)
plt.plot(timepoints, R[str(A.transcript)], label=A.transcript)
plt.plot(timepoints, R[str(A.protein)], label=A.protein)
plt.xticks(
    np.arange(0, maxtime, 3600),
    [str(i) for i in range(0, int(np.ceil(maxtime / 3600)))],
)
plt.legend()
