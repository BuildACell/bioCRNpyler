.. currentmodule:: biocrnpyler

.. _components_ref:

**********
Components
**********

Components are the primary building blocks of models in BioCRNpyler.  They
represent biomolecular parts or motifs such as promoters, enzymes,
transcriptional units, or complexes, and serve as an abstraction layer
above the core chemical species and reactions described in the
:ref:`core_ref` section. By encapsulating biological functionality in
modular, reusable objects, Components make it easier to design, share, and
manage complex models.

A key advantage of BioCRNpyler is that components are "context-aware".
They do not directly specify all the species and reactions they
generate; instead, they describe *what* they are (e.g. a promoter,
enzyme, or DNA assembly) and rely on the context provided by the
:ref:`mixtures<mixtures_ref>` to determine *how* they behave. During
compilation, each component calls one or more
:ref:`mechanims<mechanisms_ref>` from the mixture to generate the
appropriate species and reactions for that modeling context. This
design allows the same component to be reused across different models
with varying levels of detail or biological assumptions.

Defining and Using Components
------------------------------

A component typically has a `name` and may include additional attributes
such as:

- :code:`mechanisms`: A dictionary of mechanism objects it can call, which
  can override the defaults provided by the mixture.
- :code:`parameters`: A local parameter database used during compilation.
- :code:`attributes`: User-defined metadata or tags for custom behaviors.

These attributes make components flexible and allow them to adapt to
different contexts without changing their definition.

For example, consider a simple enzyme that catalyzes a conversion of
substrate to product. This can be defined as:

::

    enzyme = bcp.Enzyme(
        name='LacZ',
        substrate='Lactose',
        product='Glucose'
    )

When compiled within a mixture that supplies a catalysis mechanism, this
component will automatically generate the necessary species and reactions
describing enzymatic conversion, applying the appropriate kinetic rate
laws.

Components are usually organized into mixtures, which provide the
Mechanisms and parameter sets needed to interpret the components.
This separation between *design specification* (Components) and
*context* (mixtures) enables modular model construction and easy
exploration of different modeling assumptions. The details of
Mechanisms and mixtures are described in more detail in subsequent
chapters.

Components in BioCRNpyler provide:

- Modular, reusable building blocks for models
- Context-aware behavior via mixtures and mechanisms
- Flexible attributes for customization and parameter management
- Support for multiple mechanisms per process
- Easy extension through subclassing

By building models with components, users can create rich,
biologically grounded designs that are easy to maintain, modify, and
share.

DNA Components
---------------

One of the most common uses of components in BioCRNpyler is modeling
genetic circuits through *DNA components*. These represent transcriptional
units that specify how genes are expressed in a model.

A typical DNA component might include:

- A promoter (which initiates transcription)
- A ribosome binding site (RBS) for translation initiation
- A coding sequence specifying the protein product

For example:

::

    dna_part = bcp.DNAassembly(
        name='GFP_expression',
        promoter='P_lac',
        rbs='RBS_standard',
        protein='GFP'
    )

When compiled in a mixture that supplies transcription and translation
mechanisms, this single component will expand into the full set of species
and reactions required to model gene expression, including mRNA
intermediates and the resulting protein.

The abstraction of DNA components allows modelers to specify genetic
constructs at a high level, while mixtures determine the level of
detail—ranging from simple lumped expression models to detailed mechanistic
representations with explicit enzyme binding and multi-step processes.

Custom Components
-----------------

BioCRNpyler includes many built-in component classes to represent common
biological parts, but users can also define their own components by
subclassing :class:`~core.Component`. This extensibility is essential for
modeling novel motifs, non-standard parts, or organism-specific behaviors.

Custom components implement `~core.Component.update_species` and
`~core.Component.update_reactions` methods, which define how the component
generates its species and reactions during compilation. For example::

    class SimpleEnzyme(Component):
        def __init__(self, name, substrate, product):
            super().__init__(name)
            self.substrate = substrate
            self.product = product

        def update_species(self, ...):
            # Return list of species created
            pass

        def update_reactions(self, ...):
            # Return list of Reaction objects
            pass

By defining custom components, users can encode arbitrary logic or kinetic
forms while still integrating seamlessly with mixtures and Mechanisms.
