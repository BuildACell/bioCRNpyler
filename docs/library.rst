***********************
The BioCRNpyler Library
***********************

This chapter contains documentation on the classes and functions that
make up the BioCRNplyer package.  All objects and functions available
in the package can be access through the subpackages in which they are
contained.  For convenience, a list of low-level (core) classes,
components, mechanisms, and mixtures is also included here, in
individual sections of this chapter.

Core Classes
============

Prefix: `biocrnpyler.core`

.. automodule:: biocrnpyler.core

The core classes in BioCRNpyler define the low-level objects that are
used to specify chemical reaction networks as well as defining the
base classes for components, mechanisms as mixtures.

Chemical Reaction Networks (CRNs)
---------------------------------

Low-level chemical reaction networks can be implemented by defining
species and reactions directly.  The following classes are used to
implement this functionality.

.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:

   ChemicalReactionNetwork
   Compartment
   ModelParameter
   Reaction
   Species

Species types
-------------

A number of different species are used internally to keep track of
different types of molecular constructs.  These are normally not
accessed at the user level, but are useful when defining components
and mechanisms.

.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:

   Complex
   ComplexSpecies
   MonomerCollection
   NamedPolymer
   OrderedComplexSpecies
   OrderedMonomer
   OrderedPolymer
   OrderedPolymerSpecies
   PolymerConformation
   WeightedSpecies


Propensities
------------

    Propensities define the rate laws for chemical reactions in a CRN.
    Different propensity types implement different kinetic models such as
    mass action, Hill functions, and custom formulas. Propensities can be
    deterministic (ODE) or stochastic (Gillespie).

.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:

   Propensity
   GeneralPropensity
   Hill
   HillNegative
   HillPositive
   MassAction
   ProportionalHillPositive
   ProportionalHillNegative

Parameter databases
-------------------

Parameters are organized into databases that allow hierarchical searching.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   ParameterDatabase
   ParameterEntry
   ParameterKey


Base classes
------------

The base classes define the primary objects used to represent a model
and compile it into a chemical reaction network.  These classes are
not called directly, but are utilized for the components, mechanisms,
and mixtures that make up the BioCRNpyler library.

.. autosummary::
   :toctree: generated/
   :nosignatures:

   Component
   Mechanism
   Mixture
   Parameter

Components
==========

Prefix: `biocrnpyler.components`

.. automodule:: biocrnpyler.components

The following subsections provide a list of all components currently
available in the BioCRNpyler package.

.. include:: _autogen_components.rst


Mechanisms
==========

Prefix: `biocrnpyler.mechanisms`

.. automodule:: biocrnpyler.mechanisms

The following subsections provide a list of all mechanisms currently
available in the BioCRNpyler package.

.. include:: _autogen_mechanisms.rst


Mixtures
========

Prefix: `biocrnpyler.mixtures`

.. automodule:: biocrnpyler.mixtures

The following subsections provide a list of all mixtures currently
available in the BioCRNpyler package.

.. include:: _autogen_mixtures.rst
