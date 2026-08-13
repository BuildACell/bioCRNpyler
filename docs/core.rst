.. currentmodule:: biocrnpyler

.. _core_ref:

**************************
Chemical Reaction Networks
**************************

The lowest level representation of biological circuits and systems in
BioCRNpyler is a chemical reaction network (CRN).  While most of the
interaction with the BioCRNpyler package is done at higher levels of
abstract (:ref:`components<components_ref>`,
:ref:`mechanisms<mechanisms_ref>`, and :ref:`mixtures<mixtures_ref>`),
it is useful to understand the lower level representation and there
are occassions when it is easier to work directly at this level.


BioCRNpyler represents chemical reaction networks using
:class:`~core.Species`, :class:`~core.Reaction`, and
:class:`~core.Propensity` objects.

Species
-------

A :class:`~core.Species` represents a molecular type or state
in the model. Species can be simple chemicals, labeled complexes, or
structured entities with attributes such as name, material type, and
optional :class:`~core.Compartment`.

Key features:

- Unique string representations for human-readable output and SBML export
- Support for material types (e.g. DNA, RNA, protein)
- Integration with Compartments for spatial modeling

Example:

.. testcode::

    ATP = bcp.Species('ATP')
    RNAP = bcp.Species('RNAP', material_type='protein')
    GFP_DNA = bcp.Species('GFP', material_type='dna', compartment='cell')

Species objects are used throughout BioCRNpyler to define reactants
and products in reactions. They also support automatic naming
conventions during compilation.

The `material_type` attribute specifies the biological category of the
species, such as DNA, RNA, protein, or small molecules
(e.g. metabolites). This is important for distinguishing different
molecular roles in the model, even when species share the same base
name. For example, 'GFP' as DNA represents the gene encoding GFP,
while 'GFP' as protein is the translated product.  During compilation,
BioCRNpyler uses `material_type` to apply the correct
:ref:`mechanisms<mechanisms_ref>` (e.g. transcription for DNA,
translation for RNA) and to control naming conventions in generated
CRNs and SBML outputs. This ensures that models accurately reflect
biological processes and remain clear and interpretable, even as they
scale in complexity.

The `compartment` attribute assigns a species to a specific spatial or
logical location in the model. By default, species are placed in a
global compartment (named 'default') representing a well-mixed system,
but specifying a `compartment` allows you to distinguish between
locations such as the cytoplasm, nucleus, or vesicles. This enables
more realistic modeling of transport and localization. More details
about defining and using compartments are provided below.

Reactions
---------

A :class:`~core.Reaction` defines a chemical transformation with:

- Inputs (reactants)
- Outputs (products)
- A propensity function (kinetics)

Formally, a chemical reaction can be written as:

.. math::

  A + B \xrightarrow{k} C

where :math:`A` and :math:`B` are the input species, :math:`C` is the
output species, and :math:`k` is the reaction rate.  Formally,
BioCRNpyler represents the reaction using a "propensity function",
which defines the probability that the reaction will occur based on
the concentration of the reactants (and potentially other factors, as
described in more detail below).

Reactions in BioCRNpyler are often generated automatically during
compilation from :ref:`components<components_ref>` and
:ref:`mechanisms<mechanisms_ref>`. However, you can also define them
manually for advanced use cases or for hand-built models.

Example::

    A, B, C = bcp.Species('A'), bcp.Species('B'), bcp.Species('C')
    rxn = bcp.Reaction(
        inputs=[A, B],
        outputs=[C],
        propensity_type=bcp.MassAction(k_forward=1e-2)
    )

Key features:

- Arbitrary stoichiometry and multiple reactants/products
- Flexible propensity functions, including mass-action and Hill kinetics
- SBML export compatibility

The `propensity_type` attribute defines the kinetic rate law governing
the reaction. This specifies how the reaction rate depends on species
concentrations and parameters, such as rate constants or Hill coefficients.
BioCRNpyler includes built-in propensity types like mass-action, positive
and negative Hill functions, and custom symbolic forms. By selecting an
appropriate `propensity_type`, you control the mathematical form of the
reaction rate, enabling both detailed mechanistic models and simplified
phenomenological approximations.

Built-in propensity types:

- :class:`~core.MassAction`: classic law-of-mass-action
  kinetics
- :class:`~core.HillPositive` and
  :class:`~core.HillNegative`: for regulatory logic
- :class:`~core.ProportionalHillPositive` and
  :class:`~core.ProportionalHillNegative`: positive and negative Hill
  functions proportional to a species concentration
- :class:`~core.GeneralPropensity`: user-defined symbolic
  expressions

For example, a Hill repression propensity models regulatory effects with
a rate law of the form:

.. math::

   \rho = \frac{k}{1 + (S/K)^n}

where :math:`k` is the maximum rate, :math:`K` is the dissociation
constant, :math:`n` is the Hill coefficient, and :math:`S` is the
concentration of the repressor species.  This can be implemented in
BioCRNpyler using the following code (for a transcriptional
repressor)::

    DNA = bcp.Species('DNA', material_type='dna')
    mRNA = bcp.Species('RNA', material_type='rna')
    tetR = bcp.Species('tetR', material_type='protein')
    tx_rxn = bcp.Reaction(
        inputs=[DNA],
        outputs=[DNA, mRNA],
        propensity_type=bcp.HillNegative(
            K=10, n=2, k=0.5, s1=tetR
        )
    )

Compartments
------------

A :class:`~core.Compartment` represents a physical or
logical groupings of :class:`~core.Species`. They can model:

- Cellular compartments (e.g. cytoplasm, nucleus)
- Vesicles or other spatial domains
- Logical labels for species organization

Compartments help structure large models and support SBML
location-aware exports.

Example::

    cell = bcp.Compartment('cell')
    ATP = bcp.Species('ATP', compartment=cell)

Key features:

- Named hierarchy for clarity
- Unique SBML-compatible representations
- Seamless integration with Species

By default, all species are placed in a single well-mixed compartment
named 'default', which serves as the global container for reactions
that do not explicitly model spatial separation. Assigning species to
different compartments allows you to represent biological organization
such as cytoplasm, nucleus, or vesicles. Importantly, species with the
same name but in different compartments are treated as distinct and
will not interact with each other unless explicit transport reactions
are defined (as described in more detail in the :ref:`membranes_ref`
section).

Compiling and Simulating CRNs
------------------------------

Once you have defined species and reactions in BioCRNpyler, you can
combine them into a chemical reaction network (CRN) object for export
or simulation. This is the final, fully specified form of your model,
ready for analysis with other tools.

You can create a CRN manually by specifying the list of species and
reactions::

    crn = bcp.ChemicalReactionNetwork(
        species=[s1, s2, ..., sN],
        reactions=[r1, r2, ..., rM],
        initial_concentration_dict=initial_conditions
    )

Here, `initial_concentration_dict` is an optional dictionary mapping
species to their initial values::

    initial_conditions = {
        A: 10,
        B: 1
    }

Any species not given an explicit initial condition will default to 0.
This manual approach is useful for building small models directly or
for custom post-processing.

More typically in BioCRNpyler, you generate a CRN by *compiling* a
mixture. The mixture contains components, mechanisms, and parameters
that automatically generate all the required Species and Reactions::

    crn = mixture.compile_crn()

This method applies the mixture's mechanisms to its components, builds
all Species and Reactions, and returns a fully defined CRN object. The
resulting CRN can then be exported or simulated. More information on
components, mixtures, and mechanisms is given in subsequent chapters.

BioCRNpyler supports exporting CRNs to the SBML standard for broad
compatibility with other modeling tools. The `write_sbml_file` method
saves the CRN in SBML format to a specified file::

    crn.write_sbml_file('my_model.xml')

This file can be opened in simulation software such as COPASI or
MATLAB's SimBiology and shared with collaborators.

For direct simulation in Python, BioCRNpyler also supports integration
with the bioscrape simulator. The `simulate_with_bioscrape_via_sbml` function
runs stochastic or deterministic simulations of the CRN::

    result = crn.simulate_with_bioscrape_via_sbml(
        initial_condition_dict={A: 4, B: 1},
        timepoints=np.linspace(0, 100, 100)
    )
    plt.plot(result['time'], result[['A', 'B', 'C']])
    plt.legend(['A', 'B', 'C'])
    plt.xlabel('Time [min]')
    plt.ylabel('Concentration [nM]')

This function returns simulation results as a NumPy array or pandas
DataFrame (default), making it easy to analyze trajectories, plot
dynamics, or fit parameters.

Together, these tools make it easy to go from a high-level BioCRNpyler
design to a fully specified, exportable, and simulatable chemical
reaction network model.


Plotting Simulation Results
---------------------------

Simulation results can be plotted directly with matplotlib, as above, but
a compiled CRN often contains more species than you wrote down.  A species
that binds to anything else appears both on its own and inside every
complex that contains it, and those complexes have long generated names.
BioCRNpyler provides two plotting routines that take the CRN itself as an
argument and work out the relevant species for you.

`~utils.plotting.plot_all_species_containing` plots every species that
contains a given species, including complexes, along with the total::

    bcp.plot_all_species_containing(results, crn, ribosome)

This is the routine to use when you want to know where a species has gone.
A ribosome in a PURE model, for example, is distributed between its free
form and the translation complexes it participates in; plotting all of
them together shows the free pool being drawn down while the total stays
constant.  Passing a list plots several species on the same axes, and a
`~core.Component` may be given in place of a `~core.Species`.

`~utils.plotting.plot_gene_expression_data` summarizes the expression of
one or more genes as a set of four panels -- DNA, RNA, protein, and
resources::

    bcp.plot_gene_expression_data(results, crn, dna_assembly)

Each panel shows the individual species, any complexes containing them,
and the total.

Both routines let you control how much of each species name appears in
the legend.  `show_material` governs the species itself and the total,
while `show_complex_material` governs the complexes::

    bcp.plot_all_species_containing(
        results, crn, ribosome,
        show_material=False, show_complex_material=False
    )

The distinction matters because expanded complex names grow quickly.  A
ribosome bound to a transcript together with amino acids and ATP is
rendered as::

    complex[metabolite[AAs]:4x_metabolite[ATP]:protein[Ribo]:rna[GFP]]

which will overwhelm a legend.  Turning `show_complex_material` off
abbreviates it to ``[AAs:4x_ATP:Ribo:GFP]``.  It is off by default in
`~utils.plotting.plot_gene_expression_data`, whose four panels leave
little room, and on by default in
`~utils.plotting.plot_all_species_containing`.

`~utils.plotting.plot_all_species_containing` draws on the current axes
and sets the axis labels and legend, so it composes with the usual
matplotlib calls for figures, subplots, titles, and limits.
`~utils.plotting.plot_gene_expression_data` lays out its own two by two
grid of subplots and so takes over the current figure.

Units are handled through the `time_units` and
`concentration_units` arguments.  Time defaults to minutes in both
routines.  Concentration defaults to micromolar in
`~utils.plotting.plot_all_species_containing`, while
`~utils.plotting.plot_gene_expression_data` takes a list giving the
units for the DNA, RNA, protein, and resource panels separately,
defaulting to nanomolar, micromolar, micromolar, and millimolar::

    bcp.plot_all_species_containing(
        results, crn, ribosome,
        time_units=hrs, time_label='hr'
    )

A worked example using both routines is given in
``examples/gfp_expression.py``.
