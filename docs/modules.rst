.. currentmodule:: biocrnpyler

.. _modules_ref:

*******
Modules
*******

A module in BioCRNpyler groups the components that make up a subsystem
together with the mechanisms and parameters that subsystem needs.  Where
:ref:`components<components_ref>` describe individual biomolecular parts and
:ref:`mixtures<mixtures_ref>` describe the context those parts are compiled
in, a module describes a *piece of the design*: a sensing pathway, a
genetic circuit, a reporter.

Modules exist because models grow.  A model with thirty components in one
flat list is hard to read, hard to reuse, and hard to hand to a
collaborator.  Grouping those components into modules lets the code follow
the way the system is described on a whiteboard, and lets a subsystem that
appears twice be written once.

A `~core.Module` is itself a `~core.Component`, so it can be placed
anywhere a component can go, including inside another module.  A module
cannot be compiled on its own, since it describes a design rather than a
context; it has to be placed in a mixture::

    crn = Mixture('test', components=[my_module]).compile_crn()

Defining and Using Modules
--------------------------

A module is defined by specifying:

- A `name` identifying the subsystem.
- The components belonging to it.
- The mechanisms its components need (optional).
- The parameters describing it (optional), as a dictionary or a file.
- Species to add directly to the CRN (optional).

For example, a reporter subsystem might be written as::

    reporter = Module(
        name='reporter',
        components=[
            DNAassembly('GFP', promoter='pGFP', rbs='rGFP')
        ],
        mechanisms={
            'transcription': SimpleTranscription(),
            'translation': SimpleTranslation()
        },
        parameters={'ktx': 0.05, 'ktl': 0.01}
    )

Everything describing the reporter is now in one place.  Note that a module
holds the same kinds of things a mixture does, with one exception: global
mechanisms, which act on every species in a model, remain a property of the
mixture.

Combining Modules and Mixtures
------------------------------

Modules are combined using ``+``.  Adding a module to a mixture returns a
new mixture containing it::

    system = PURE(name='chassis') + reporter + sensor

The original mixture is not modified, so a chassis mixture can be reused
across several systems.  Two modules can also be added to each other,
producing a module that contains both::

    detector = reporter + sensor
    system = PURE(name='chassis') + detector

When two modules are combined this way, each module's mechanisms and
parameters are applied to its own components first, so the combined module
keeps the two subsystems' settings separate.

Wiring Modules Together
-----------------------

Modules are connected by *name*.  Components in different modules that
share a type and a name compile to the same species, and are compiled only
once.  A protein produced in one module and consumed in another is simply
written with the same name in both::

    sensing = Module('sensing', components=[Protein('signal'), ...])
    response = Module('response', components=[Protein('signal'), ...])

Here both modules refer to a single ``protein_signal`` in the compiled CRN.
This is what makes modules composable: there is no separate wiring step,
and no need for one module to hold a reference to another.

The same applies between a module and its mixture.  A module that declares
a component the mixture already provides, such as an RNA polymerase, will
use the mixture's copy rather than creating a second one.

.. note::

   A component declared by more than one module is compiled under the
   context of whichever module is enumerated first.  If two modules
   need different mechanisms or parameters for a component, then they
   do not really share it, and it should be renamed in one of them
   using `~core.Module.instance` (described below).

Mechanisms and Parameters in Modules
------------------------------------

A module's mechanisms and parameters apply to its own components.  This
means two modules in the same mixture can model the same process
differently: one subsystem can use simple transcription while another uses
a Michaelis-Menten model with explicit RNA polymerase, in a single
compiled CRN.

Precedence runs from the inside out:

1. A component's own mechanisms and parameters
2. Its module's
3. The mixture's

So a parameter written directly on a component is used in preference to
one from its module, which is used in preference to one from the mixture.
Mechanisms follow the same order, with a further level below the mixture
for defaults, described in `Default mechanisms`_ below.

This ordering is by *level*, not by how specifically a parameter is keyed.
The :ref:`parameter defaulting<parameters_ref>` rules, which prefer a more
specific `~core.ParameterKey` over a general one, apply *within* a level.
Between levels, the inner level wins even if the outer one names the
parameter more precisely::

    enzyme = Enzyme('betagal', 'XGal', 'indigo', parameters={'kcat': 1})

    readout = Module(
        'readout', components=[enzyme],
        parameters={('michaelis_menten', 'betagal', 'kcat'): 2}
    )

The enzyme's ``kcat`` of 1 is used, even though the module names the
parameter more precisely, because the enzyme is the inner level.

An important consequence is that an outer level cannot override an
inner one.  A parameter set on a mixture will not replace one already set
on a module or a component, however specifically it is written.  To change
a value, set it where it was defined::

    reporter.update_parameters(
        parameters={'ktx': 0.1}, overwrite_parameters=True)

Generally speaking, a parameter should only be specified within a module
if the value of that parameter is fixed in all possible compositions.
Otherwise, leave the parameter value unspecified and set it in the
final mixture.

Note that leaving a parameter unspecified does not oblige anything to
supply it later.  If neither the module nor the mixture gives a value,
the mechanism's own default is used and the model compiles without
complaint, so a module should document the parameters it expects the
mixture to provide.

Default mechanisms
^^^^^^^^^^^^^^^^^^

Mechanisms, unlike parameters, can be given as a *suggestion* rather
than as a requirement.  A component may carry a
`default_mechanism`, which is consulted only after the mixture has been
asked, so it never prevents an outer level from choosing something
else.  A `~components.ChemicalComplex`, for instance, defaults to
`~mechanisms.One_Step_Binding` but will use whatever binding mechanism
its module or mixture specifies.

A module can offer a default in the same way::

    signaling = Module(
        'signaling', components=[...],
        default_mechanism=One_Step_Cooperative_Binding()
    )

For mechanisms, then, the full order is:

1. A component's own mechanisms
2. Its module's mechanisms
3. The mixture's mechanisms
4. Its module's default mechanism
5. The component's own default mechanism

A module's default mechanism replaces the one a component carries from
its own class, on the grounds that the class default is a library
fallback rather than a statement about this particular model.  Where
modules are nested, the innermost module's default is the one offered.

This distinction is worth keeping in mind when deciding what to put in
a module.  A default mechanism is safe to specify, since any outer
level overrides it; a mechanism given in the ``mechanisms`` argument,
like a parameter, is binding on every composition the module appears
in.

Modules and Parameter Files
---------------------------

Because precedence runs inward, the simplest way to manage a large model
is to leave parameters out of the components and modules entirely and put
them in a single parameter file supplied to the mixture.  Modules then
describe the *structure* of the design and the file describes its
*numbers*.

A system defined this way compiles before any parameters are supplied at
all, using the defaults carried by the mechanisms themselves::

    reporter = Module('reporter', components=[DNAassembly('GFP', ...)])
    sensor = Module('sensor', components=[DNAassembly('RFP', ...)])

    system = Mixture('system', mechanisms=..., ) + reporter + sensor
    crn = system.compile_crn()      # compiles using mechanism defaults

Supplying a parameter file to the mixture then sets individual reaction
rates and initial conditions throughout the model, reaching into the
modules and their components::

    mechanism_id            part_id  param_name  param_val
    simple_transcription    pGFP     ktx         2.0
    simple_transcription    pRFP     ktx         7.0
    initial concentration   system   dna_GFP     4.0

Note that reaction parameters are addressed by `part_id`, which for gene
expression is the name of the promoter or RBS rather than the name of the
gene or the module.  Two genes that share a promoter name share its
parameters.

Reusing a Module
----------------

A module can be used more than once in the same mixture by making copies
of it with `~core.Module.instance`.  Each copy is given a name and a
dictionary of the names to change::

    s1 = signaling.instance('s1', rename={'ligand': 'IPTG'})
    s2 = signaling.instance('s2', rename={'ligand': 'aTc'})

    system = PURE(name='chassis') + s1 + s2

Only the names listed are changed.  Anything left alone keeps its name and
is therefore shared between the copies, which is how the copies stay
connected to common parts of the system: if the module also contains a
kinase that is not renamed, both copies use one kinase.

Renaming applies to a component's name, to the species it holds, including
the species inside a complex, and to any parameter keys that refer to it by
name.  All of these matter.  A component's species determine what it
compiles to, and a component's initial concentration is stored under its
own name, so a rename that changed only the component's name would quietly
produce a different model.  Renaming reaches into nested modules as well.

A name in the rename dictionary that matches nothing raises an error rather
than being ignored, since a mistyped name would otherwise leave two copies
silently sharing a component they were meant to separate.

Naming a module's interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since renaming changes every occurrence of a name, it is worth naming the
parts of a module that copies are expected to rebind for their *role*
rather than for whatever the first copy happens to use.  A membrane pore
that transports a different molecule in each copy is clearer, and safer,
written with a placeholder::

    pore_channel = MembraneChannel(
        pore_monomer.product, substrate='cargo', ...)

    pore_AHL = pore.instance('pore_AHL', rename={'cargo': 'AHL', ...})
    pore_XGal = pore.instance('pore_XGal', rename={'cargo': 'XGal', ...})

Had the template used ``'AHL'`` as the substrate name directly, renaming
it to ``'XGal'`` for the second copy would also have renamed every *other*
occurrence of AHL in that module, including the AHL bound inside an
activator complex driving the promoter.

Remember also that parameters are addressed by part name.  If two copies of
a module need different reaction rates from a parameter file, the parts
carrying those rates, such as the promoter and the RBS, have to be renamed
too, not just the gene.

Nested Modules
--------------

Because a module is a component, a module can contain other modules.  A
larger subsystem can therefore be assembled from smaller ones::

    detector_cell = Module(
        'detector_cell',
        components=[qs_sensing, membrane_pore, readout]
    )

Precedence continues to run from the inside out: a component's own
settings, then those of the module holding it, then those of any enclosing
module, and finally the mixture's.
