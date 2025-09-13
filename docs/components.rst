.. currentmodule:: biocrnpyler

.. _components_ref:

**********
Components
**********

Components are the primary building blocks of models in BioCRNpyler. 
They represent biomolecular parts or motifs such as promoters, 
enzymes, transcriptional units, or complexes, and serve as an abstraction 
layer above the core chemical species and reactions described in the 
:ref:`core_ref` section. By encapsulating biological functionality in 
modular, reusable objects, Components make it easier to design, share, 
and manage complex models.

A key advantage of BioCRNpyler is that Components are "context-aware".
They do not directly specify all the species and reactions they
generate; instead, they describe *what* they are (e.g. a promoter,
enzyme, or DNA assembly) and rely on the context provided by the
:ref:`mixtures_ref` to determine *how* they behave. During
compilation, each Component calls one or more
:ref:`mechanims<mechanisms_ref>` from the Mixture to generate the
appropriate species and reactions for that modeling context. This
design allows the same Component to be reused across different models
with varying levels of detail or biological assumptions.

Defining and Using Components
------------------------------

A Component typically has a `name` and may include additional attributes 
such as:

- `mechanisms`: A dictionary of Mechanism objects it can call, which 
  can override the defaults provided by the Mixture.
- `parameters`: A local parameter database used during compilation.
- `attributes`: User-defined metadata or tags for custom behaviors.

These attributes make Components flexible and allow them to adapt to 
different contexts without changing their definition.

For example, consider a simple enzyme that catalyzes a conversion of 
substrate to product. This can be defined as:

::

    enzyme = bcp.Enzyme(
        name='LacZ',
        substrate='Lactose',
        product='Glucose'
    )

When compiled within a Mixture that supplies a catalysis Mechanism, 
this Component will automatically generate the necessary species and 
reactions describing enzymatic conversion, applying the appropriate 
kinetic rate laws.

Components are usually organized into Mixtures, which provide the
Mechanisms and parameter sets needed to interpret the Components.
This separation between *design specification* (Components) and
*context* (Mixtures) enables modular model construction and easy
exploration of different modeling assumptions. The details of
Mechanisms and Mixtures are described in more detail in subsequent
chapters.

Components in BioCRNpyler provide:

- Modular, reusable building blocks for models
- Context-aware behavior via Mixtures and Mechanisms
- Flexible attributes for customization and parameter management
- Support for multiple Mechanisms per process
- Easy extension through subclassing

By building models with Components, users can create rich,
biologically grounded designs that are easy to maintain, modify, and
share.

DNA Components
---------------

One of the most common uses of Components in BioCRNpyler is modeling 
genetic circuits through *DNA Components*. These represent 
transcriptional units that specify how genes are expressed in a model. 

A typical DNA Component might include:

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

When compiled in a Mixture that supplies transcription and translation 
Mechanisms, this single Component will expand into the full set of 
species and reactions required to model gene expression, including 
mRNA intermediates and the resulting protein.

The abstraction of DNA Components allows modelers to specify genetic 
constructs at a high level, while Mixtures determine the level of 
detail—ranging from simple lumped expression models to detailed 
mechanistic representations with explicit enzyme binding and 
multi-step processes.

Custom Components
-----------------

BioCRNpyler includes many built-in Component classes to represent 
common biological parts, but users can also define their own Components 
by subclassing :class:`~core.component.Component`. This extensibility 
is essential for modeling novel motifs, non-standard parts, or 
organism-specific behaviors.

Custom Components implement `update_species` and `update_reactions` 
methods, which define how the Component generates its species and 
reactions during compilation. For example::

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

By defining custom Components, users can encode arbitrary logic or 
kinetic forms while still integrating seamlessly with Mixtures and 
Mechanisms.

API Reference
-------------

More information on Components can be found in the following module:

.. autosummary::

   biocrnpyler.core.component
