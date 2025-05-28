#############################################################
BioCRNPyler - Biomolecular Chemical Reaction Network Compiler
#############################################################

BioCRNpyler (pronounced "bio compiler") is a Python package for the
creation, manipulation, and study of the structure, dynamics, and
functions of complex biochemical networks.  BioCRNpyler compiles
high-level design specifications for biochemical processes --
represented using a modular library of biochemical parts, mechanisms,
and contexts -- to CRN implementations.

.. rubric:: Features

- Automated constructions of CRN representations to Systems Biology
  Markup Language (SBML) models compatible with numerous simulators
- Library of modular biochemical components allows for different
  architectures and implementations of biochemical circuits to be
  represented succinctly with design choices propagated throughout the
  underlying CRN automatically
- High-level design specifications can be embedded into diverse
  biomolecular environments, such as cell-free extracts and *in vivo*
  systems
- Includes a parameter database that allows users to rapidly prototype
  large models using very few parameters that can be customized later
  in the development cycle

.. rubric:: Links:

- Mailing list: [SBTools Google
  Group](https://groups.google.com/g/sbtools/) Email:
  sbtools@googlegroups.com
- Source code: https://github.com/BuildACell/BioCRNPyler
- Overview paper: [BioCRNpyler: Compiling Chemical Reaction Networks
  from Biomolecular Parts in Diverse
  Contexts](https://doi.org/10.1371/journal.pcbi.1009987)
- Bug reports: https://github.com/BuildACell/BioCRNPyler/issues
- Slack: Join the #biocrnpyler channel on SBTools slack: Ask on the
  public SBTools Google group to be added or send a message to one of
  the maintainers.

.. rubric:: How to cite

An `article <https://doi.org/10.1371/journal.pcbi.1009987>`_ about the
package is available in PLOS Computational Biology. If the BioCRNpyler
package helped you in your research, please cite::

  @article{biocrnpyler2022,
    title={{BioCRNpyler}: {Compiling} chemical reaction networks from
           biomolecular parts in diverse contexts},
    author={W. Poole and A. Pandey and A. Shur and Z. A. Tuza and
            R. M. Murray},
    journal = {PLOS Computational Biology},
    year = {2022},
    month = {04},
    volume = {18},
    url = {https://doi.org/10.1371/journal.pcbi.1009987},
    pages = {1-19}
  }

or the GitHub site: https://github.com/BuildACell/BioCRNPyler.


.. toctree::
   :caption: User Guide
   :maxdepth: 1
   :numbered: 2

   intro
   tutorial
   membrane_features
   Reference Manual <reference>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
