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
can be overriden at the component, mixture, or Part level.

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
degradation, transcripton, etc), the part ID string indicates the
specific part (or None if the parameter value should be used for all
parts), and the parameter string is the name of the parameter used by
the component or mixture.

Parameters can be specified directly as dictionaries, but they are
commonly stored in files.  Parameter files can be .tsv (or .txt) or
.csv. The first line of the file should contain column headings. The
following headings are required (in any order): mechanism_id, part_id,
param_name, and param_val (spaces can be substituted for underscores
and headings are not case sensitive). The mechanism_id is the name of
the mechanism or the kind of mechanism that will use this parameter,
for example 'transcription' or 'transcription_mm' for Mechalis-Menten
transcription would go in this column. The part_id refers to the name
supplied by the component that will use this mechanism, for example
'ptet' for a TetR repressed promoter and 'ptet_leak' for leak
reactions of that promoter. The param_name field refers to the name of
the model parameter, for example 'ktx', 'kb', or 'ku'. The value of
each entry is case sensitive and underscores are different from
spaces. All these parameter headings (if they are not empty) are
required to be alpha-numeric (plus single internal underscores) and
start with a letter.

The easiest way to examine which parameter values are used in a Model
and where BioCRNpyler found them is with the pretty_print function of
a `ChemicalReactionNetwork` or a `Reaction` with `show_rates` and
`show_keys` both set to True::

  CRN.pretty_print(show_rates=True, show_keys=True)

The `show_rates` keyword toggles whether reaction rates are shown. The
parameter values used in the rates will be shown below the rate
function. The `show_keys` keyword toggles whether the ParameterKey
used to search and find these rates is displayed. The printed output
will show the search_key for the mechanism and component and indicate
the found_key where the parameter was found after applying the
defaulting mechanism.

Parameter Value Defaulting
==========================

Parameters can be passed into components and mixtures as both files
and/or dictionaries. Parameters into components will be used by
reactions created by that component:

    Component(... parameters = parameters, parameter_file = "file_name" ...)

If a component cannot find its own parameters, it will use parameters
passed into the mixture:

    Mixture(... parameters=parameters, parameter_file='file_name', ...)

At both the mixture and component levels, BioCRNpyler does nto require
an exact keyword match. When mechanisms are called by components to
produce reactions, they will initially look for

    ParameterKey(
        mechanism='mechanism_name', part_id='part_id',
	name='parameter_name') --> param_val

If that particular parameter key cannot be found, the software will
default to the following keys in this order:

    ParameterKey(
        mechanism='mechanism_type', part_id='part_id', name='parameter_name')
    ParameterKey(
        mechanism=None, part_id='part_id' , name='parameter_name')
    ParameterKey(
        mechanism='mechanism_name', part_id=None , name='parameter_name')
    ParameterKey(
        mechanism='mechanism_type', part_id=None , name='parameter_name')
    ParameterKey(
        mechanism= None, part_id=None , name='parameter_name')

As a note, mechanism_name refers to the .name variable of a
Mechanism. mechanism_type refers to the .type variable of a
Mechanism. Either of these can be used as a mechanism_id. This allows
for models to be constructed easily using default parameter values and
for parameters to be shared between different mechanisms and/or
Components.

Parameters passed into a component supercede mixture level parameters
in reactions produced by that component.

Components and mixtures can both have one more multiple parameter
files by passing in a list of filenames instead of a single filename
to the parameter_file keyword. Components use parameters loaded from
their file(s) before defaulting to the file(s) supplied to a
mixture. The last file in any list will take precedent and overwrite
parameter files which were written earlier.

When parameters are specified via parameter files, a set of
directories is search to find the first matching filename:

    1. The current directory
    2. The list of directores in the `BCP_PATH` environment variable
    3. The directory containing the BioCRNpyler source code files

Parameter values for commonly used components are included in the
source code distribution.  The following parameter files are available:

.. list-table::
   :header-rows: 1

   * - Path
     - Description
   * - /components/tetr_parameters.tsv
     - Promoter and binding parameters for TetR and aTc
   * - mechanisms/transport_parameters.tsv
     - Membrane transport parameters
   * - mixtures/pure_parameters.tsv
     - Gene expression parameters for PURE


Parameter Units
===============

The Parameter database in BioCRNpyler can also be used to import units
for model parameters. Simply write your units in the parameter file
(.csv/.tsv) as a new column titled "unit". To manually create a
Parameter object with a unit, you can run::

    new_parameter = Parameter(
        parameter_name='None', parameter_value='1.0', unit='M')

BioCRNpyler supports the following unit definitions. The string in
paranthesis is the unit identifier that you should use for that
particular unit:

1. Time units: second ("second"), minute ("minute"), hour ("hour").
2. Concentration units: nanomolar ("nM"), micromolar ("uM"),
   millimolar ("mM"), molar ("M").
3. Volume units: nanolitre ("nL"), microlitre ("uL"), millilitre ("mL"),
   litre ("L").
4. Common parameter units: per second ("per_second"), per hour ("per_hour"),
   1/(molar x second) ("litre_per_mole_per_second"),
   1/(molar x hour) ("litre_per_mole_per_hour")

To add support for any new unit of your choice, simply go to
`utils/units.py` in the `biocrnpyler/` source code and add your own
function definition for a unit (you can use the existing functions to
get an idea about how to create your own new definition).
