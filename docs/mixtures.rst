.. currentmodule:: biocrnpyler

.. _mixtures_ref:

********
Mixtures
********

A mixture in BioCRNpyler defines the *context* in which components are
compiled into a chemical reaction network (CRN). While
:ref:`components<components_ref>` describe *what* biomolecular parts are
present, and :ref:`mechanisms<mechanisms_ref>` define *how* biological
processes are modeled, the mixture ties these together by specifying
*which* Mechanisms are available, *which* components are present, and
*what* parameters to use.

This separation of design from context is central to BioCRNpyler's
flexibility. The same set of components can be compiled in different
mixtures to produce models with varying levels of detail, different kinetic
assumptions, or different biological environments (e.g. cell-free extract
vs. in vivo).

Defining and Using Mixtures
----------------------------

A mixture is defined by specifying:

- A `name` identifying the modeling context.
- A list of components to include.
- A dictionary of mechanisms to make available during compilation.
- A parameter database shared by all components and mechanisms.

During compilation, each component in the mixture calls the relevant
Mechanisms to generate its species and reactions. The mixture ensures
consistent parameter handling, mechanism availability, and naming
conventions across the entire model.

For example, a simple cell-free system might include transcription and
translation mechanisms, along with components representing DNA assemblies:

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

When compiled, this mixture will automatically generate all species and
reactions needed to model gene expression under the specified assumptions.

Controlling Model Detail
-------------------------

One of the most powerful features of mixtures is that they determine *model
resolution*. By changing which mechanisms are included, you can quickly
shift between simplified and detailed representations without changing your
components.

For example, you can define gene expression in two ways depending on the
mechanisms you supply in the mixture. A very simple model might use a
single combined mechanism:

::

    simple_mixture = Mixture(
        name='simple_expression',
        mechanisms={
            'transcription': OneStepGeneExpression(),
	    'translation': EmptyMechanism()
        }
    )

This would generate a one-step reaction:

.. math::

   DNA \rightarrow DNA + Protein

For more detail, you can use separate mechanisms to model transcription,
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

This approach models RNA explicitly, includes degradation pathways, and
accounts for loss due to cell growth and division.

Parameters and Defaults
------------------------

Mixtures also manage model parameters, providing a central database that
components and mechanisms can query during compilation. This parameter
database stores values such as rate constants, Hill coefficients, or
binding affinities in a structured way.

Parameters in BioCRNpyler are identified by *ParameterKeys*, which specify:

- The mechanism type (e.g. 'transcription', 'translation')
- The parameter name (e.g. 'k', 'K', 'n')
- The component name or part name (optional for specificity)

When a mechanism needs a parameter value during compilation, BioCRNpyler
uses a defaulting hierarchy to search the parameter database. The search
tries to find the most specific match first, falling back to more general
entries if needed. This allows you to define highly specific parameters for
a single component, as well as broad defaults that apply across the entire
model.

For example, suppose you define the following parameters in a mixture:

::

    params = {
        ('transcription', 'k', 'GFP_expression'): 0.2,
        ('transcription', 'k', None): 0.1,
        ('translation', 'k'): 0.5
    }

Here:

- The rate 0.2 will be used for transcription of the 'GFP_expression'
  component specifically.
- The rate 0.1 will be used for any other transcription component
  without its own specific value.
- The translation rate is 0.5 for all translation mechanisms.

You can supply this database when defining the mixture:

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

Additionally, components can have their own local parameter databases,
which override the mixture's parameters for that specific component.  This
design lets you easily manage complex parameter sets while maintaining
clear, reusable defaults across the entire model.

Predefined Mixtures
--------------------

BioCRNpyler includes several predefined mixture classes designed to
represent common experimental contexts. These mixtures come with
appropriate mechanisms, components, and default parameters already
configured, making it easy to set up standard modeling scenarios
quickly. Users can either use these as-is or subclass them to create custom
variations.

Extract-Based Mixtures
~~~~~~~~~~~~~~~~~~~~~~~

Extract-based mixtures are designed to model cell-free transcription-
translation (TX-TL) systems, such as in vitro expression reactions.  These
environments lack the complexity of living cells but are widely used for
prototyping circuits and parts. Extract mixtures typically include gene
expression mechanisms, global dilution or degradation, and suitable
parameter sets.

BioCRNpyler includes several extract-based mixture classes:

- `~mixtures.ExpressionExtract`: Uses `~mechanisms.OneStepGeneExpression`
  for very simple models where transcription and translation are collapsed
  into a single step.

- `~mixtures.SimpleTxTlExtract`: Includes separate transcription and
  translation mechanisms to explicitly model mRNA as an intermediate.

- `~mixtures.EnergyTxTlExtract`: Models transcription and translation with
  explicit representation of RNA polymerase (RNAP), ribosomes, RNases,
  and energy carrier molecules.

For example, you can create a simple extract-based mixture in code as::

    extract_mixture = SimpleExtract(
        name='cell_free_extract',
        components=[dna_part]
    )

This mixture typically includes the following mechanisms:

+----------------+----------------------------------------+
| Mechanism Type | Description                            |
+================+========================================+
| expression     | One-step gene expression combining     |
|                | transcription and translation          |
+----------------+----------------------------------------+
| dilution       | Global mechanism modeling loss due to  |
|                | cell-free extract degradation or decay |
+----------------+----------------------------------------+

During compilation, it will find all DNA species in the model and generate
protein production reactions for them automatically.

This mixture automatically includes the OneStepGeneExpression mechanism for
all DNA components, a global dilution mechanism to model loss over time,
and a parameter database appropriate for cell-free systems. During
compilation, it will find all DNA species in the model and generate protein
production reactions for them automatically.

You can also create an extract-based mixture with
separate transcription and translation mechanisms using
ExpressionExtract::

    extract_mixture = SimpleTxTl(
        name='cell_free_expression',
        components=[dna_part]
	parameter_file='mixtures/extract_parameters.tsv'
    )

This mixture includes the following mechanisms:

+----------------+-------------------------------------------+
| Mechanism Type | Description                               |
+================+===========================================+
| transcription  | Simple transcription mechanism generating |
|                | mRNA from DNA                             |
+----------------+-------------------------------------------+
| translation    | Simple translation mechanism producing    |
|                | protein from mRNA                         |
+----------------+-------------------------------------------+
| dilution       | Global mechanism modeling loss due to     |
|                | extract degradation or decay              |
+----------------+-------------------------------------------+

During compilation, this mixture will create mRNA intermediates and
model gene expression as two sequential steps, allowing more detailed
exploration of transcriptional and translational dynamics.

Extract-based mixtures make it easy to move between simple and detailed
cell-free models simply by swapping which mixture subclass you use. This
modular approach is ideal for exploring different levels of model
complexity without changing your component definitions.

Cell-Based Mixtures
~~~~~~~~~~~~~~~~~~~

Cell-based mixtures in BioCRNpyler are designed to model gene expression in
living cells, capturing features like transcription and translation, RNA
and protein degradation, and dilution due to cell growth and
division. These mixtures configure mechanisms and parameters reflecting
typical in vivo environments.

BioCRNpyler includes the following predefined cell-based mixture classes:

- `~mixtures.ExpressionDilutionMixture`: Simplified one-step gene expression
  with global dilution and degradation.

- `~mixtures.SimpleTxTlDilutionMixture`: Two-step gene expression with
  explicit mRNA intermediates and dilution.

- `~mixtures.TxTlDilutionMixture`: More detailed mechanistic expression
  with regulated transcription and enzyme-mediated kinetics, plus dilution.

For example, you can create a simplified in vivo-like mixture with
`~mixtures.ExpressionDilutionMixture`::

    extract_mixture = ExpressionDilutionMixture(
        name='cell_based_expression',
        components=[dna_part],
	parameter_file='mixtures/cell_parameters.tsv'
    )

This mixture includes the following mechanisms:

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

During compilation, this mixture identifies all DNA species and applies the
expression mechanism to generate protein production reactions directly. It
also includes global degradation and dilution to reflect loss of molecules
over time in growing cells.

You can also use `~mixtures.SimpleTxTlDilutionMixture` to model gene
expression as two sequential steps with explicit mRNA intermediates::

    simple_txtl_mixture = SimpleTxTlDilutionMixture(
        name='simple_txl_dilution',
        components=[dna_part]
    )

This mixture includes the following mechanisms:

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

During compilation, this mixture finds all DNA components and
systematically applies transcription and translation mechanisms to create
mRNA and protein species, along with degradation and dilution processes for
realistic in vivo dynamics.

For more detailed modeling, :class:`~mixtures.TxTlDilutionMixture`
includes regulated transcription and enzyme-mediated translation::

    txtl_mixture = TxTlDilutionMixture(
        name='detailed_txl_dilution',
        components=[dna_part]
    )

This mixture includes the following mechanisms:

+-------------------------+-----------------------------------------+
| Mechanism Type          | Description                             |
+=========================+=========================================+
| transcription           | Regulated transcription (Hill-based)    |
+-------------------------+-----------------------------------------+
| translation             | Translation mechanism with optional     |
|                         | enzyme dynamics                         |
+-------------------------+-----------------------------------------+
| rna_degradation         | Degradation of mRNA species             |
+-------------------------+-----------------------------------------+
| protein_degradation     | Degradation of protein species          |
+-------------------------+-----------------------------------------+
| dilution                | Global loss due to cell growth          |
+-------------------------+-----------------------------------------+

This mixture enables modeling of regulatory control, enzyme-mediated
kinetics, and system-wide dilution, providing a richer, more detailed
representation of gene expression in cells.

PURE Mixtures
~~~~~~~~~~~~~

The PURE (Protein synthesis Using Recombinant Elements) system is a
reconstituted cell-free system, assembled from purified components rather
than from a cell extract.  Because the composition is known, PURE models
can account explicitly for the machinery and the consumable resources that
limit expression.

A single mixture, `~mixtures.PURE`, covers the whole range of PURE models.
Rather than choosing between separate classes, you select the level of
detail with three switches::

    pure_mixture = PURE(
        name='pure',
        components=[dna_part],
        include_machinery=True,
        include_resources=True,
        include_fuel=True,
    )

The switches are nested: each one requires the ones above it.

+---------------------+------------------------------------------------+
| Argument            | Effect when True                               |
+=====================+================================================+
| include_machinery   | Adds RNA polymerase and ribosome, so that      |
|                     | expression competes for a finite pool of       |
|                     | machinery                                      |
+---------------------+------------------------------------------------+
| include_resources   | Adds NTP and amino acid consumption, so that   |
|                     | transcription and translation draw down a      |
|                     | finite supply of building blocks               |
+---------------------+------------------------------------------------+
| include_fuel        | Adds the fuel species (ATP) and energy         |
|                     | utilization, so that expression halts when the |
|                     | system runs out of energy                      |
+---------------------+------------------------------------------------+

Turning all three off gives a minimal model in which expression proceeds
without competition or depletion.  Turning them on one at a time is a
convenient way to see which limitation dominates in a given design: with
machinery alone, two genes compete for ribosomes; adding resources makes
expression taper as amino acids are consumed; adding fuel makes it stop
altogether once ATP is exhausted.

The species names can be changed if the defaults clash with the rest of
your model::

    pure_mixture = PURE(
        name='pure', components=[dna_part], ribosome='Ribosome'
    )

`~mixtures.BasicPURE` remains available and is equivalent to
`~mixtures.PURE` with all three switches set to True.


Custom Mixtures
----------------

BioCRNpyler includes built-in mixture classes for common contexts such as
cell-free systems and in vivo environments, each pre-configured with
appropriate mechanisms and parameters. However, you can also define your
own custom mixture classes by subclassing :class:`~core..Mixture` and
overriding its setup to include your choice of mechanisms, components, and
parameters.

This extensibility enables modeling of specialized biological systems or
custom experimental setups.
