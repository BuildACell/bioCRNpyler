.. currentmodule:: biocrnpyler

.. _mechanisms_ref:

**********
Mechanisms
**********

Mechanisms in BioCRNpyler define "reaction schemas" that describe the
biochemical processes generating species and reactions during model
compilation. They sit between the abstract design of
:ref:`components<components_ref>` and the concrete chemical reactions
and species described in the :ref:`core_ref` section.

A mechanism specifies *how* a given biological process is modeled: for
example, transcription might be represented as a simple lumped
reaction or as a detailed enzyme-mediated mechanism. This flexibility
lets users explore different modeling assumptions by swapping
mechanisms without changing the components themselves.

Mechanisms work in conjunction with :ref:`mixtures<mixtures_ref>`,
which provide the context for compilation. Each mixture holds a set of
available mechanisms and parameters. During compilation, components
call the mechanisms in their mixture to generate the species and
reactions appropriate to that modeling context.

Defining and Using Mechanisms
-----------------------------

Mechanisms are designed to be modular and reusable. Each mechanism has
a `name` and a `mechanism_type` describing the category of process it
implements, such as "transcription", "translation", "binding", or
"degradation". These types let components specify *what* process they
need without knowing *how* it will be modeled.

For example, consider modeling transcription in two ways:

- A simple, lumped reaction:

  .. math::

     \text{DNA} \rightarrow \text{DNA} + \text{RNA}

- A mechanistic Michaelis-Menten form:

  .. math::

     \text{DNA} + \text{RNAP} \leftrightharpoons
         \text{DNA:RNAP} \rightarrow \text{DNA} + \text{RNAP} + \text{RNA}

By changing which mechanism is supplied in the mixture, the same DNA
component will compile into different sets of reactions, enabling easy
switching between modeling abstractions.

Mechanisms are typically created in advance and added to a mixture.
For example::

    tx_simple = bcp.SimpleTranscription()
    tx_mm = bcp.Transcription_MM(rnap=Species('RNAP'))

These mechanisms can then be supplied to different mixtures to set the
desired level of detail.

By defining processes as mechanisms, BioCRNpyler enables models to
remain modular and flexible. This approach allows users to test
different assumptions and contexts without rewriting components,
supports multiple modeling abstractions, and ensures a clear
separation between process type and implementation.

Core Mechanisms for Gene Expression
------------------------------------

BioCRNpyler includes a library of mechanisms designed to model gene
expression and other fundamental processes. These can be easily
swapped to change the level of detail in a model. Key categories
include:

- *Transcription*: Converts DNA to RNA. Includes simple forms without
  enzyme detail, Michaelis-Menten kinetics with explicit RNA
  polymerase binding, and Hill-function-based regulatory repression or
  activation.

- *Translation*: Converts RNA to protein. Supports simple mass-action
  forms and Michaelis-Menten variations modeling ribosome binding.

- *Degradation*: Represents molecular decay or enzymatic
  breakdown. Usually modeled as first-order mass-action reactions, but
  custom mechanisms can specify more complex pathways.

- *Dilution*: Accounts for loss of molecules due to cell growth and
  division. Typically implemented as first-order processes that remove
  species over time at a defined rate.

For example, a mixture might specify both a SimpleTranscription
Mechanism and a SimpleTranslation mechanism to create a streamlined
gene expression model, or replace them with more mechanistic
alternatives to capture enzyme dynamics explicitly.

DNAassembly Components
~~~~~~~~~~~~~~~~~~~~~~

Modeling gene expression in BioCRNpyler typically begins with the
:class:`~components.dna.assembly.DNAassembly` component, which
represents a transcriptional unit specifying how a gene is
expressed. A :class:`~components.dna.assembly.DNAassembly` combines
several biological parts in a single, modular definition:

- *Promoter*: defines where transcription initiates and its regulatory
  properties.
- *Ribosome Binding Site (RBS)*: controls translation initiation
  efficiency.
- *Coding Sequence (CDS)*: specifies the protein product to be produced.
- *Terminator* (optional): marks the end of transcription.

When a :class:`~components.dna.assembly.DNAassembly` is compiled in a
mixture, it uses the available mechanisms to generate the full set of
species and reactions describing gene expression. For example, it will
call the transcription mechanism to produce an mRNA species from the
DNA template and the translation mechanism to produce a protein from
the mRNA.

The species naming follows clear conventions to ensure consistency:

- DNA species: '<gene_name>_DNA'
- mRNA species: '<gene_name>_RNA'
- Protein species: '<protein_name>_protein'

This design makes it easy to trace which DNA parts produce which
outputs in large, multi-gene models while maintaining modularity and
reusability. By swapping different mechanisms in the mixture, the same
:class:`~components.dna.assembly.DNAassembly` component can represent
anything from coarse-grained one-step gene expression to detailed
mechanistic models with enzyme binding and multiple regulatory states.

A typical :class:`~components.dna.assembly.DNAassembly` can be defined
in code as::

    dna_part = bcp.DNAassembly(
        name='GFP_expression',
        promoter='P_lac',
        rbs='RBS_standard',
        protein='GFP'
    )

This definition specifies a promoter, RBS, and protein product, with
default or optional settings for other attributes such as terminators.
When compiled in a mixture with suitable mechanisms, this single
Component will expand into all required species and reactions to model
expression of GFP from the DNA template.

One Step Gene Expression
~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~mechanisms.txtl.OneStepGeneExpression` mechanism
implements a highly simplified model of gene expression by collapsing
transcription and translation into a single reaction. This is useful
for coarse-grained models where intermediate mRNA species are not
represented explicitly.

The reaction schema is:

.. math::

   \text{DNA} \xrightarrow{k} \text{DNA} + \text{Protein}

For example, you can create this mechanism in code as::

    expression_mechanism = bcp.OneStepGeneExpression(
        name='simple_expression'
    )

This mechanism can then be added to a mixture to define the overall
gene expression process.

When used during compilation, this mechanism relies on the mixture to
find all Species with `material_type` set to "DNA". For each such DNA
species, the mixture calls the expression mechanism to generate a
reaction producing the corresponding Protein. The mechanism examines
the DNA component's `protein` attribute to determine which Protein
species to create, using consistent naming
(e.g. '<protein_name>_protein').  If that Protein species does not
exist, it is automatically created.  This design allows the mechanism
to be applied systematically across all DNA parts in the system.

Required Parameters:

+----------------+--------------------------------------------+
| Parameter Name | Description                                |
+================+============================================+
| k              | Overall expression rate constant           |
+----------------+--------------------------------------------+

Simple Transcription and Translation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~mechanisms.txtl.SimpleTranscription` and
:class:`~mechanisms.txtl.SimpleTranslation` mechanisms represent gene
expression as two sequential steps, explicitly including mRNA as an
intermediate species. This approach is common in models that want to
capture transcriptional and translational regulation while retaining a
simple mass-action form.

The reactions are:

Transcription:

.. math::

   DNA \xrightarrow{k_\text{tx}} DNA + mRNA

Translation:

.. math::

   mRNA \xrightarrow{k_\text{tl}} mRNA + Protein

You can create these mechanisms in code as::

    tx_mechanism = bcp.SimpleTranscription(name='simple_tx')
    tl_mechanism = bcp.SimpleTranslation(name='simple_tl')

These mechanisms can then be added to a mixture to model two-step
gene expression with mRNA as an intermediate.

During compilation, these mechanisms use the mixture to apply
transcription and translation to all DNA components. The transcription
mechanism creates an mRNA species named '<gene_name>_RNA' using the
DNA component's `name`.  The translation mechanism produces a Protein
species named '<protein_name>_protein' from the DNA component's
`protein` attribute.  Species are automatically added if they do not
already exist, ensuring consistent, traceable naming across complex
models.

Required Parameters:

+----------------+--------------------------------------------+
| Parameter Name | Description                                |
+================+============================================+
| k_tx           | Transcription rate constant                |
+----------------+--------------------------------------------+
| k_tl           | Translation rate constant                  |
+----------------+--------------------------------------------+

Regulated Transcription
~~~~~~~~~~~~~~~~~~~~~~~~~

The `PositiveHillTranscription` and `NegativeHillTranscription`
Mechanisms model transcriptional regulation using Hill-function-based
rate laws. These are suitable for describing cooperative binding of
activators or repressors.

The general reaction schema is:

.. math::

   DNA \xrightarrow{\rho} DNA + mRNA

With rate law:

- For activation:

  .. math::

     \rho = \frac{k \cdot S^n}{K^n + S^n}

- For repression:

  .. math::

     \rho = \frac{k}{1 + (S/K)^n}

You can create these mechanisms using::

    pos_reg_tx = bcp.PositiveHillTranscription(name='activated_tx')
    neg_reg_tx = bcp.NegativeHillTranscription(name='repressed_tx')

These mechanisms can then be added to a mixture to model transcriptional
activation and repression using Hill functions.

During compilation, the mixture applies these mechanisms to DNA
components and identifies the regulator Species based on user
specification or naming conventions. This ensures the correct
transcription factor or repressor is included in the model. Species
are automatically created as needed.

Required Parameters:

+----------------+---------------------------------------------+
| Parameter Name | Description                                 |
+================+=============================================+
| k              | Maximum transcription rate                  |
+----------------+---------------------------------------------+
| K              | Dissociation constant for regulator binding |
+----------------+---------------------------------------------+
| n              | Hill coefficient (degree of cooperativity)  |
+----------------+---------------------------------------------+

Michaelis-Menten Based Mechanisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Transcription_MM` and `Translation_MM` mechanisms provide more
mechanistic representations using Michaelis-Menten kinetics,
explicitly modeling binding and unbinding of RNA polymerase and
ribosomes.

Transcription_MM:

.. math::

   DNA + RNAP \leftrightharpoons DNA:RNAP \rightarrow DNA + RNAP + mRNA

Translation_MM:

.. math::

   mRNA + Ribo \leftrightharpoons mRNA:Ribo \rightarrow mRNA + Ribo + Protein

You can create these mechanisms using::

    tx_mm = bcp.Transcription_MM(name='tx_mm', rnap=Species('RNAP'))
    tl_mm = bcp.Translation_MM(name='tl_mm', ribosome=Species('Ribo'))

These mechanisms can then be added to a mixture to capture enzyme-mediated
binding and catalysis.

The mixture ensures these mechanisms are applied to all appropriate
DNA components and uses the defined RNAP and Ribo Species. Bound
complexes are automatically named and created during compilation.

Required Parameters:

+----------------+--------------------------------------------+
| Parameter Name | Description                                |
+================+============================================+
| kb             | Binding rate constant                      |
+----------------+--------------------------------------------+
| ku             | Unbinding rate constant                    |
+----------------+--------------------------------------------+
| kc             | Catalytic rate constant                    |
+----------------+--------------------------------------------+

Energy Utilization Mechanisms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `Energy_Transcription_MM` and `Energy_Translation_MM` mechanisms
extend the Michaelis-Menten forms by modeling consumption of an energy
Species during the catalytic step.

Example transcription reaction:

.. math::

   DNA:RNAP + Energy \rightarrow DNA + RNAP + mRNA

Example translation reaction:

.. math::

   mRNA:Ribo + Energy \rightarrow mRNA + Ribo + Protein

You can create these mechanisms using::

    energy_tx_mm = Energy_Transcription_MM(
        name='energy_tx',
        rnap=Species('RNAP'),
        energy=Species('ATP')
    )

    energy_tl_mm = Energy_Translation_MM(
        name='energy_tl',
        ribosome=Species('Ribo'),
        energy=Species('GTP')
    )

These mechanisms can then be added to a mixture to model energy-dependent
transcription and translation.

During compilation, the mixture applies these mechanisms across DNA
Components and ensures the specified energy Species is used in the
reactions. This allows linking gene expression to energy availability
in the system.

Required Parameters:

+----------------+--------------------------------------------+
| Parameter Name | Description                                |
+================+============================================+
| kb             | Binding rate constant                      |
+----------------+--------------------------------------------+
| ku             | Unbinding rate constant                    |
+----------------+--------------------------------------------+
| kc             | Catalytic rate constant                    |
+----------------+--------------------------------------------+

Detailed Models Including Isomerization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The `multi_tx` and `multi_tl` mechanisms offer the most detailed
representations, modeling multiple binding, unbinding, and
isomerization steps in transcription and translation.

Transcription steps:

.. math::

   \begin{aligned}
     {} & \text{DNA} + \text{RNAP} \leftrightharpoons \text{DNA:RNAP}_{closed}
       \rightarrow \text{DNA:RNAP}_{open} \\
     {} & \text{DNA:RNAP}_{open} \rightarrow \text{DNA} + \text{RNAP}
       + \text{mRNA}
   \end{aligned}

Translation steps:

.. math::

   \begin{aligned}
     {} & \text{mRNA} + \text{Ribo} \leftrightharpoons
       \text{mRNA:Ribo}_\text{bound}
       \rightarrow \text{mRNA:Ribo}_\text{active} \\
     {} & \text{mRNA:Ribo}_\text{active} \rightarrow
       \text{mRNA} + \text{Ribo} + \text{Protein}
   \end{aligned}

You can create these mechanisms using::

    multi_tx_mechanism = multi_tx(name='detailed_tx', rnap=Species('RNAP'))
    multi_tl_mechanism = multi_tl(name='detailed_tl', ribosome=Species('Ribo'))

These mechanisms can then be added to a mixture to model multi-step
binding and isomerization events in detail.

During compilation, the mixture applies these mechanisms to all DNA
components, creating intermediate complex Species with clear
naming. This enables modeling of delays, cooperativity, and detailed
enzyme kinetics.

Required Parameters:

+----------------+--------------------------------------------+
| Parameter Name | Description                                |
+================+============================================+
| kb1            | Initial binding rate constant              |
+----------------+--------------------------------------------+
| ku1            | Initial unbinding rate constant            |
+----------------+--------------------------------------------+
| kc1            | Isomerization rate constant                |
+----------------+--------------------------------------------+
| kc2            | Catalytic rate constant                    |
+----------------+--------------------------------------------+

These mechanisms are suited for models requiring a detailed,
biophysically informed representation of gene expression.

Global Mechanisms
------------------

In addition to the mechanisms called by individual components,
mixtures can include *global mechanisms* that operate on all species
in the model.  Global mechanisms are applied after all components have
generated their species and reactions, allowing them to enforce
model-wide processes like dilution, decay, or transformation.

Typical uses for global mechanisms include:

- Modeling *dilution* of all species due to cell growth and division.
- Applying *background degradation* to any unprotected molecules.
- Adding *global conversion* or labeling rules for all species.

When compiling a mixture, global mechanisms receive the entire set of
species generated so far and can add additional species and reactions
as needed. This design ensures consistent, system-wide behavior
without requiring every component to implement these processes
individually.

For example, a mixture might include a dilution mechanism that
automatically adds first-order decay reactions for all cellular
species:

::

    dilution_mechanism = Dilution(rate=0.01)

    detailed_mixture = Mixture(
        name='detailed_expression',
        components=[dna_part],
        mechanisms={
            'transcription': SimpleTranscription(),
            'translation': SimpleTranslation(),
            'dilution': dilution_mechanism
        }
    )

This approach keeps components modular and focused on their own
biological roles while ensuring that important global processes are
included in the compiled CRN.

Global mechanisms can also use the material_type attribute of Species
to apply processes selectively. For example, a degradation mechanism
might only add decay reactions for species with a material type of
"protein" or "RNA", leaving DNA or small molecules untouched. This
filtering allows global degradation processes to accurately model
biological specificity while still being applied automatically across
the entire set of species in the model.


Custom Mechanisms
------------------

Users can also define custom mechanisms by subclassing
:class:`~core.mechanism.Mechanism`. A mechanism must implement
`update_species` and `update_reactions` methods, which specify the
species and reactions to add during compilation.

Example::

    class SimpleDegradation(Mechanism):
        def __init__(self, name='simple_degradation'):
            super().__init__(name, mechanism_type='degradation')

        def update_species(self, ...):
            # Return list of species created
            pass

        def update_reactions(self, ...):
            # Return list of Reaction objects
            pass

This extensibility allows users to implement new kinetic forms,
organism-specific processes, or advanced logic for generating CRNs.
