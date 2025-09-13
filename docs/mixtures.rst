.. currentmodule:: biocrnpyler

.. _mixtures_ref:

********
Mixtures
********

A Mixture in BioCRNpyler defines the *context* in which Components are 
compiled into a chemical reaction network (CRN). While 
:ref:`components<components_ref>` describe *what* biomolecular parts are 
present, and :ref:`mechanisms<mechanisms_ref>` define *how* biological 
processes are modeled, the Mixture ties these together by specifying *which* 
Mechanisms are available, *which* Components are present, and *what* 
parameters to use.

This separation of design from context is central to BioCRNpyler's 
flexibility. The same set of Components can be compiled in different 
Mixtures to produce models with varying levels of detail, different 
kinetic assumptions, or different biological environments (e.g. cell-free 
extract vs. in vivo).

Defining and Using Mixtures
----------------------------

A Mixture is defined by specifying:

- A `name` identifying the modeling context.
- A list of Components to include.
- A dictionary of Mechanisms to make available during compilation.
- A parameter database shared by all Components and Mechanisms.

During compilation, each Component in the Mixture calls the relevant 
Mechanisms to generate its species and reactions. The Mixture ensures 
consistent parameter handling, mechanism availability, and naming 
conventions across the entire model.

For example, a simple cell-free system might include transcription and 
translation Mechanisms, along with Components representing DNA 
assemblies:

::

    tx = SimpleTranscription()
    tl = SimpleTranslation()

    dna_part = DNAassembly(
        name='GFP_expression',
        promoter='P_lac',
        rbs='RBS_standard',
        protein='GFP'
    )

    simple_mixture = Mixture(
        name='cell_free',
        components=[dna_part],
        mechanisms={
            'transcription': tx,
            'translation': tl
        }
    )

When compiled, this Mixture will automatically generate all species and 
reactions needed to model gene expression under the specified 
assumptions.

Controlling Model Detail
-------------------------

One of the most powerful features of Mixtures is that they determine 
*model resolution*. By changing which Mechanisms are included, you 
can quickly shift between simplified and detailed representations 
without changing your Components.

For example, you can define gene expression in two ways depending on 
the Mechanisms you supply in the Mixture. A very simple model might 
use a single combined mechanism:

::

    simple_mixture = Mixture(
        name='simple_expression',
        mechanisms={
            'transcription': OneStepGeneExpression(),
	    'translation`: EmptyMechanism()
        }
    )

This would generate a one-step reaction:

.. math::

   DNA \rightarrow DNA + Protein

For more detail, you can use separate Mechanisms to model transcription,  
translation, mRNA degradation, protein degradation, and dilution:

::

    detailed_mixture = Mixture(
        name='detailed_expression',
        mechanisms={
            'transcription': SimpleTranscription(),
            'translation': SimpleTranslation(),
            'rna_degradation': SimpleDegradation(),
            'protein_degradation': SimpleDegradation(),
            'dilution': Dilution()
        }
    )

This approach models RNA explicitly, includes degradation pathways,  
and accounts for loss due to cell growth and division.

Parameters and Defaults
------------------------

Mixtures also manage model parameters, providing a central database 
that Components and Mechanisms can query during compilation. This 
parameter database stores values such as rate constants, Hill 
coefficients, or binding affinities in a structured way.

Parameters in BioCRNpyler are identified by *ParameterKeys*, which 
specify:

- The Mechanism type (e.g. 'transcription', 'translation')
- The parameter name (e.g. 'k', 'K', 'n')
- The Component name or part name (optional for specificity)

When a Mechanism needs a parameter value during compilation, BioCRNpyler 
uses a defaulting hierarchy to search the parameter database. The search 
tries to find the most specific match first, falling back to more general 
entries if needed. This allows you to define highly specific parameters 
for a single Component, as well as broad defaults that apply across the 
entire model.

For example, suppose you define the following parameters in a Mixture:

::

    params = {
        ('transcription', 'k', 'GFP_expression'): 0.2,
        ('transcription', 'k', None): 0.1,
        ('translation', 'k'): 0.5
    }

Here:

- The rate 0.2 will be used for transcription of the 'GFP_expression' 
  Component specifically.
- The rate 0.1 will be used for any other transcription Component 
  without its own specific value.
- The translation rate is 0.5 for all translation Mechanisms.

You can supply this database when defining the Mixture:

::

    mix = Mixture(
        name='cell_free',
        components=[dna_part],
        mechanisms={
            'transcription': SimpleTranscription(),
            'translation': SimpleTranslation()
        },
        parameters=params
    )

Additionally, Components can have their own local parameter databases, 
which override the Mixture's parameters for that specific Component. 
This design lets you easily manage complex parameter sets while 
maintaining clear, reusable defaults across the entire model.

Predefined Mixtures
--------------------

BioCRNpyler includes several predefined Mixture classes designed to 
represent common experimental contexts. These Mixtures come with 
appropriate Mechanisms, Components, and default parameters already 
configured, making it easy to set up standard modeling scenarios 
quickly. Users can either use these as-is or subclass them to create 
custom variations.

Extract-Based Mixtures
~~~~~~~~~~~~~~~~~~~~~~~

Extract-based Mixtures are designed to model cell-free transcription-
translation (TX-TL) systems, such as in vitro expression reactions. 
These environments lack the complexity of living cells but are widely 
used for prototyping circuits and parts. Extract Mixtures typically 
include gene expression Mechanisms, global dilution or degradation, 
and suitable parameter sets.

BioCRNpyler includes several extract-based Mixture classes:

- `Extract`: A flexible base class for defining extract systems.
- `SimpleExtract`: Uses OneStepGeneExpression for very simple models 
  where transcription and translation are collapsed into a single step.
- `ExpressionExtract`: Includes separate transcription and translation 
  Mechanisms to explicitly model mRNA as an intermediate.
- `CombinatorialExtract`: Supports combinatorial assembly of DNA 
  parts and dynamic generation of components through enumeration.

For example, you can create a simple extract-based Mixture in code as:

::

    extract_mixture = SimpleExtract(
        name='cell_free_extract',
        components=[dna_part]
    )

This Mixture typically includes the following Mechanisms:

+----------------+----------------------------------------+
| Mechanism Type | Description                            |
+================+========================================+
| expression     | One-step gene expression combining     |
|                | transcription and translation          |
+----------------+----------------------------------------+
| dilution       | Global mechanism modeling loss due to  |
|                | cell-free extract degradation or decay |
+----------------+----------------------------------------+

During compilation, it will find all DNA species in the model and
generate protein production reactions for them automatically.    

This Mixture automatically includes the OneStepGeneExpression Mechanism 
for all DNA Components, a global dilution Mechanism to model loss over 
time, and a parameter database appropriate for cell-free systems. During 
compilation, it will find all DNA species in the model and generate 
protein production reactions for them automatically.

You can also create an extract-based Mixture with
separate transcription and translation Mechanisms using
ExpressionExtract::

    extract_mixture = ExpressionExtract(
        name='cell_free_expression',
        components=[dna_part]
    )
    
This Mixture includes the following Mechanisms:

+----------------+------------------------------------------+
| Mechanism Type | Description                              |
+================+==========================================+
| transcription  | SimpleTranscription Mechanism generating |
|                | mRNA from DNA                            |
+----------------+------------------------------------------+
| translation    | SimpleTranslation Mechanism producing    |
|                | protein from mRNA                        |
+----------------+------------------------------------------+
| dilution       | Global mechanism modeling loss due to    |
|                | extract degradation or decay             |
+----------------+------------------------------------------+

During compilation, this Mixture will create mRNA intermediates and
model gene expression as two sequential steps, allowing more detailed
exploration of transcriptional and translational dynamics.

Extract-based Mixtures make it easy to move between simple and detailed 
cell-free models simply by swapping which Mixture subclass you use. This 
modular approach is ideal for exploring different levels of model 
complexity without changing your Component definitions.

Cell-Based Mixtures
~~~~~~~~~~~~~~~~~~~~~~~

Cell-based Mixtures in BioCRNpyler are designed to model gene 
expression in living cells, capturing features like transcription 
and translation, RNA and protein degradation, and dilution due to 
cell growth and division. These Mixtures configure Mechanisms and 
parameters reflecting typical in vivo environments.

BioCRNpyler includes the following predefined cell-based Mixture 
classes:

- `ExpressionDilutionMixture`: Simplified one-step gene expression 
  with global dilution and degradation.
- `SimpleTxTlDilutionMixture`: Two-step gene expression with explicit 
  mRNA intermediates and dilution.
- `TxTlDilutionMixture`: More detailed mechanistic expression with 
  regulated transcription and enzyme-mediated kinetics, plus dilution.

For example, you can create a simplified in vivo-like Mixture with 
`ExpressionDilutionMixture`::

This Mixture includes the following Mechanisms:

+-------------------------+-----------------------------------------+
| Mechanism Type          | Description                             |
+=========================+=========================================+
| expression              | One-step gene expression combining      |
|                         | transcription and translation           |
+-------------------------+-----------------------------------------+
| rna_degradation         | Degradation of mRNA species             |
+-------------------------+-----------------------------------------+
| protein_degradation     | Degradation of protein species          |
+-------------------------+-----------------------------------------+
| dilution                | Global loss due to cell growth          |
+-------------------------+-----------------------------------------+

During compilation, this Mixture identifies all DNA species and applies 
the expression Mechanism to generate protein production reactions 
directly. It also includes global degradation and dilution to reflect 
loss of molecules over time in growing cells.

You can also use `SimpleTxTlDilutionMixture` to model gene expression 
as two sequential steps with explicit mRNA intermediates::

    simple_txtl_mixture = SimpleTxTlDilutionMixture(
        name='simple_txl_dilution',
        components=[dna_part]
    )

This Mixture includes the following Mechanisms:

+-------------------------+-----------------------------------------+
| Mechanism Type          | Description                             |
+=========================+=========================================+
| transcription           | SimpleTranscription generating mRNA     |
+-------------------------+-----------------------------------------+
| translation             | SimpleTranslation producing protein     |
+-------------------------+-----------------------------------------+
| rna_degradation         | Degradation of mRNA species             |
+-------------------------+-----------------------------------------+
| protein_degradation     | Degradation of protein species          |
+-------------------------+-----------------------------------------+
| dilution                | Global loss due to cell growth          |
+-------------------------+-----------------------------------------+

During compilation, this Mixture finds all DNA Components and 
systematically applies transcription and translation Mechanisms 
to create mRNA and protein species, along with degradation and 
dilution processes for realistic in vivo dynamics.

For more detailed modeling, `TxTlDilutionMixture` includes regulated 
transcription and enzyme-mediated translation::

    txtl_mixture = TxTlDilutionMixture(
        name='detailed_txl_dilution',
        components=[dna_part]
    )

This Mixture includes the following Mechanisms:

+-------------------------+-----------------------------------------+
| Mechanism Type          | Description                             |
+=========================+=========================================+
| transcription           | Regulated transcription (Hill-based)    |
+-------------------------+-----------------------------------------+
| translation             | Translation Mechanism with optional     |
|                         | enzyme dynamics                         |
+-------------------------+-----------------------------------------+
| rna_degradation         | Degradation of mRNA species             |
+-------------------------+-----------------------------------------+
| protein_degradation     | Degradation of protein species          |
+-------------------------+-----------------------------------------+
| dilution                | Global loss due to cell growth          |
+-------------------------+-----------------------------------------+

This Mixture enables modeling of regulatory control, enzyme-mediated 
kinetics, and system-wide dilution, providing a richer, more detailed 
representation of gene expression in cells.


Custom Mixtures
----------------

BioCRNpyler includes built-in Mixture classes for common contexts 
such as cell-free systems and in vivo environments, each pre-configured 
with appropriate Mechanisms and parameters. However, you can also define 
your own custom Mixture classes by subclassing 
:class:`~core.mixture.Mixture` and overriding its setup to include 
your choice of Mechanisms, Components, and parameters.

This extensibility enables modeling of specialized biological systems 
or custom experimental setups.

API Reference
-------------

More information on Mixtures can be found in the following module:

.. autosummary::

   biocrnpyler.core.mixture
