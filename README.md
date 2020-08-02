# BioCRNPyler -- Biomolecular Chemical Reaction Network Compiler
## Python toolbox to create CRN models in SBML for biomolecular mechanisms

[![Build Status](https://travis-ci.com/BuildACell/BioCRNPyler.svg?branch=master)](https://travis-ci.com/BuildACell/BioCRNPyler)
[![codecov](https://codecov.io/gh/BuildACell/BioCRNPyler/branch/master/graph/badge.svg)](https://codecov.io/gh/BuildACell/BioCRNPyler)
[![PyPI version](https://badge.fury.io/py/biocrnpyler.svg)](https://badge.fury.io/py/biocrnpyler)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/BuildACell/BioCRNPyler/master?filepath=%2Fexamples%2F)

BioCRNPyler is a Python package for the creation, manipulation,
and study of the structure, dynamics, and functions
of complex networks.

- **Website:** https://github.com/BuildACell/BioCRNPyler
- **Mailing list:** TBA
- **Source:** https://github.com/BuildACell/BioCRNPyler
- **Bug reports:** https://github.com/BuildACell/BioCRNPyler/issues
- **Documentation** [biocrnpyler.readthedocs.io](https://readthedocs.org/projects/biocrnpyler/)

# Simple example

Building a simple reaction network

```python
from biocrnpyler import *
# let's build the following CRN
# A -->[k1] 2B
# B -->[k2] B+D
# Species
A = Species("A")
B = Species("B")
C = Species("C")
D = Species("D")

#Reaction Rates
k1 = 3.
k2 = 1.4

#Reaciton Objects
R1 = Reaction.from_massaction([A], [B, B], k_forward = k1)
R2 = Reaction.from_massaction([B], [C, D], k_forward = k2)

#Make a CRN
CRN = ChemicalReactionNetwork(species = [A, B, C, D], reactions = [R1, R2])
print(CRN)
```
More advanced examples can be found in the [example](https://github.com/BuildACell/BioCRNPyler/tree/master/examples) folder, 
here's the first file in the Tutorial series: [Building CRNs](https://github.com/BuildACell/BioCRNPyler/blob/master/examples/1.%20Building%20CRNs%20Directly.ipynb)

# Installation


Install the latest version of BioCRNPyler::

    $ pip install biocrnpyler

Install with all optional dependencies::

    $ pip install biocrnpyler[all]

Further details about the installation process can be found in the [BioCRNPyler wiki](https://github.com/BuildACell/BioCRNPyler/wiki#installation).
# Bugs
Please report any bugs that you find [here](https://github.com/BuildACell/BioCRNPyler/issues).
Or, even better, fork the repository on [GitHub](https://github.com/BuildACell/BioCRNPyler),
and create a pull request (PR). We welcome all changes, big or small, and we
will help you make the PR if you are new to `git` (just ask on the issue and/or
see `docs/CONTRIBUTING.md`).

# Versions

BioCRNpyler versions:

* 1.0.0 (latest stable release): To install run `pip install biocrnpyler` 
* 0.2.1 (alpha release): To install run `pip install biocrnpyler==0.2.1`

# License
Released under the BSD 3-Clause License (see `LICENSE`)

Copyright (c) 2020, Build-A-Cell. All rights reserved.

