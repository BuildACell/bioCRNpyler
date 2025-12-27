.. currentmodule:: biocrnpyler

.. _parameters_ref:

**********
Parameters
**********

BioCRNpyler has a flexible parameter system allowing for models to be
rapidly produced using very few parameters common across many
reactions and then refined so that mechanisms, components, and parts
have unique reaction parameters derived from the literature or
experiments. This parameter system is designed to work automatically,
so all users have to do is specify a parameter file of the proper
format or a parameter dictionary.  Parameters have default values that
can be overridden at the component, mixture, or part level.

Describing Parameters
=====================

Parameters are defined according to the mechanism, part name, and
parameter name.  This allows the same parameter symbol to be re-used
in different mechanisms and allows different values for the same
parameter depending on the particular part.  So, as an example, there
might be an overall degradation parameter that holds for all RNA in a
mixture, but that is modified for a given type of RNA (e.g., an
aptamer with strong secondary structure) and is further modified for
specific types of aptamers (differentiated based on the part name).

Internally, parameters are stored as a dictionary, where the key is a
tuple::

  parameters = {
      ('mechanism_1', 'part_id', 'parameter_A'): value,
      ('mechanism_2', 'part_id', 'parameter_B'): value,
      ...
  }

The mechanism string describes the underlying process (RNA
degradation, transcription, etc), the part ID string indicates the
specific part (or None if the parameter value should be used for all
parts), and the parameter string is the name of the parameter used by
the component or mixture.

Parameters can be specified directly as dictionaries, but they are
commonly stored in files.  Parameter files can be .tsv (or .txt) or
.csv. BioCRNpyler will strip out all lines starting with a '#' to
allow commenting of parameter files.  The first (non-comment) line of
the file should contain column headings. The following headings are
required (in any order): `mechanism_id`, `part_id`, `param_name`, and
`param_val` (spaces can be substituted for underscores and headings are
not case sensitive).

The `mechanism_id` is the mechanism type or the name of the specific
mechanism that will use this parameter. For example, parameters
associated with all transcription mechanisms would set `mechanism_id`
to 'transcription', while parameters specific to Michaelis-Menten
transcription (`~mechanisms.Transcription_MM`) are specified using
'transcription_mm' for `mechanism_id`.

The `part_id` refers to the name supplied by the component that will
use this mechanism, for example 'ptet' for a TetR repressed promoter
and 'ptet_leak' for leak reactions of that promoter.  Some components
also use a part type to define parameters for different classes of
reactions and the `part_id` will be checked for this part type as
well.  An example is the `~components.dna.RegulatedPromoter` component,
which uses the 'dna_protein' part type for the 'binding' mechanism,
and the `~components.ChemicalComplex` component, which uses the
'chemical_complex' part type for 'binding' reactions.  A named
`part_id` has priority over the part type.

The `param_name` field refers to the name of the model parameter, for
example 'ktx', 'kb', or 'ku'. The value of each entry is case
sensitive and underscores are different from spaces.

All parameter headings (if they are not empty) are required to
be alpha-numeric (plus single internal underscores) and start with a
letter.

The easiest way to examine which parameter values are used in a model
and where BioCRNpyler found them is with the pretty_print function of
a `~core.ChemicalReactionNetwork` or a `~core.Reaction` with
`show_rates` and `show_keys` both set to True::

  CRN.pretty_print(show_rates=True, show_keys=True)

The `show_rates` keyword toggles whether reaction rates are shown. The
parameter values used in the rates will be shown below the rate
function. The `show_keys` keyword toggles whether the `~core.ParameterKey`
used to search and find these rates is displayed. The printed output
will show the search_key for the mechanism and component and indicate
the found_key where the parameter was found after applying the
defaulting mechanism.

Parameter Value Defaulting
==========================

Parameters can be passed into components and mixtures as both files
and/or dictionaries. Parameters into components will be used by
reactions created by that component::

    Component(... parameters=parameters, parameter_file='file_name', ...)

If a component cannot find its own parameters, it will use parameters
passed into the mixture::

    Mixture(... parameters=parameters, parameter_file='file_name', ...)

At both the mixture and component levels, BioCRNpyler does not require
an exact keyword match. When mechanisms are called by components to
produce reactions, they will initially look for::

    ParameterKey(
        mechanism='mechanism_name', part_id='part_id',
	name='parameter_name') --> param_val

If that particular parameter key cannot be found, the software will
default to the following keys in this order::

    ParameterKey(
        mechanism='mechanism_type', part_id='part_id', name='parameter_name')
    ParameterKey(
        mechanism=None, part_id='part_id' , name='parameter_name')
    ParameterKey(
        mechanism='mechanism_type', part_id='part_type', name='parameter_name')
    ParameterKey(
        mechanism=None, part_id='part_type' , name='parameter_name')
    ParameterKey(
        mechanism='mechanism_name', part_id=None , name='parameter_name')
    ParameterKey(
        mechanism='mechanism_type', part_id=None , name='parameter_name')
    ParameterKey(
        mechanism= None, part_id=None , name='parameter_name')

The 'mechanism_name' field refers to the `name` attribute of a
`~core.Mechanism` while the 'mechanism_type' refers to the `type`
attribute of a `~core.Mechanism`. Either of these can be used as a
mechanism ID.  Similarly, the 'part_id' field refers to the name
assigned by the `~core.Component` that creates the reactions while the
`part_type` refers to the type of reaction (e.g. 'dna_protein' versus
'chemical_complex' for the 'binding' mechanisms). These naming
conventions allow for models to be constructed using default
parameter values and for parameters to be shared between different
mechanisms and/or components.

Parameters passed into a component supersede mixture level parameters
in reactions produced by that component.

Components and mixtures can both have multiple parameter
files by passing in a list of filenames instead of a single filename
to the parameter_file keyword. Components use parameters loaded from
their file(s) before defaulting to the file(s) supplied to a
mixture. The last file in any list will take precedent and overwrite
parameter files which were written earlier.

When parameters are specified via parameter files, a set of
directories is searched to find the first matching filename:

    1. The current directory
    2. The list of directories in the `BCP_PATH` environment variable
    3. The directory containing the BioCRNpyler source code files

Parameter values for commonly used components are included in the
source code distribution.  The following parameter files are available:

.. list-table::
   :header-rows: 1

   * - Path
     - Description
   * - components/tetr_parameters.tsv
     - Promoter and binding parameters for TetR and aTc
   * - mechanisms/binding_parameters.tsv
     - Binding and unbinding rates, including cooperative binding
   * - mechanisms/transport_parameters.tsv
     - Membrane transport parameters
   * - mixtures/cell_parameters.tsv
     - Gene expression parameters for dilution mixtures (cells)
   * - mixtures/extract_parameters.tsv
     - Gene expression parameters for cell-based extracts
   * - mixtures/pure_parameters.tsv
     - Gene expression parameters for PURE


Parameter Units
===============

The parameter database in BioCRNpyler can also be used to import units
for model parameters. Simply write your units in the parameter file
(.csv/.tsv) as a new column titled "unit". To manually create a
Parameter object with a unit, you can run::

    new_parameter = Parameter(
        parameter_name='None', parameter_value='1.0', unit='M')

BioCRNpyler supports the following unit definitions. The string in
parenthesis is the unit identifier that you should use for that
particular unit:

1. Time units: second ('second' or 'sec'), minute ('minute' or 'min'),
   hour ('hour' or 'hrs').
2. Concentration units: nanomolar ('nM'), micromolar ('uM'),
   millimolar ('mM'), molar ('M').
3. Volume units: nanolitre ('nL'), microlitre ('uL'), millilitre ('mL'),
   litre ('L').
4. Reciprocal units: precede the unit name with '/' or 'per\_'.
   E.g. 'nM/sec' or 'nM per_second'.

The default units are 'uM', 'uL', and 'sec'.  To allow easy
conversation from other units, the `utils.units` module contains a set of
conversion factors that can be loaded using the following import
command::

  from biocrnpyler.utils.units import (
      nM, uM, mM, M, nL, uL, mL, L, sec, min, hrs
  )

With these variables in place, conversions to different units can be
handled by multiplication and division.  For example::

  initial_conditions = {'dna_GFP': 1*nM}
  plt.plot(results['time']/min, results['protein_GFP'/uM])

To add support for any new unit of your choice, simply go to
`utils/units.py` in the `biocrnpyler/` source code and add your own
function definition for a unit (you can use the existing functions to
get an idea about how to create your own new definition).

.. note::

   Units are not yet full implemented in BioCRNpyler.  While they can
   be specified when defining parameters, unit specifications are not
   yet implemented in conversion to SBML, which can result in errors
   if inconsistent units are used.  Parameters specified in files
   distributed with BioCRNpyler use uM, uL, and sec for all parameter
   values.
