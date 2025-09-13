***********************
The BioCRNpyler Library
***********************

This chapter contains documentation on the classes and functions that
make up the BioCRNplyer package.  All objects and functions available
in the package can be access through the subpackages in which they are
contained.  For convenience, a list of low-level (core) classes,
components, mechanisms, and mixtures is also included here, in
individual sections of this chapter.

.. automodule:: biocrnpyler
    :members:
    :undoc-members:
    :show-inheritance:

Subpackages
===========

The BioCRNpyler package is organized as a set of subpackages that
define the differrent functions and objects used to model a system.
Information on the individual subpackages can be accessed via the
table below.  Information on all components, mechanisms, and mixtures
contained in the package are described further down on the page.

.. currentmodule:: biocrnpyler

.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:

   core
   components
   mechanisms
   mixtures
   utils

Low-Level Classes
=================

Low-level chemical reaction networks can be implemented by defining
species and reactions directly.  The following classes are used to
implement this functionality.

.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:

   ChemicalReactionNetwork
   Complex
   Parameter
   Reaction
   Species

More detailed information about specialized subclasses for
representing species can be found in the
:mod:`biocrnpyler.core.species` module.

Components
==========

.. automodule:: biocrnpyler.components
    :members:
    :undoc-members:
    :show-inheritance:

The following subsections provide a list of all components currently
available in the BioCRNpyler package.

.. toctree::
    :maxdepth: 2
    :caption: Components

.. include:: _autogen_components.rst

Mechanisms
==========

.. automodule:: biocrnpyler.mechanisms
    :members:
    :undoc-members:
    :show-inheritance:

The following subsections provide a list of all mechanisms currently
available in the BioCRNpyler package.

.. toctree::
    :maxdepth: 2
    :caption: Mechanisms

.. include:: _autogen_mechanisms.rst

Mixtures
========

.. automodule:: biocrnpyler.mixtures
    :members:
    :undoc-members:
    :show-inheritance:

The following subsections provide a list of all mixtures currently
available in the BioCRNpyler package.

.. toctree::
    :maxdepth: 2
    :caption: Mixtures

.. include:: _autogen_mixtures.rst

