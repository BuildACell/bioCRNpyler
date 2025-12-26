# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

"""Binding mechanisms for chemical complexes.

The 'binding' mechanisms are used to create the reactions required to
implement binding between two or more chemical species.  Binding occurs
between one or more binding species and a bindee species.  Cooperative binding
involves the use of multimers with a given multiplicity, and binding can be in
either one step (a single reaction for all multimers) or two steps (multimers
form a complex, then bind).

"""

import itertools as it

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction
from ..core.species import Complex, Species, WeightedSpecies


class One_Step_Cooperative_Binding(Mechanism):
    r"""Cooperative binding mechanism for single-step multi-ligand binding.

    A 'binding' mechanism where multiple copies of a binder molecule (A) bind
    cooperatively to a single bindee molecule (B) in one concerted step.  This
    models cooperative binding events where all ligands bind simultaneously
    rather than sequentially.

    The binding reaction is given by

        $$ n 'A' + 'B' <--> 'A'_n:'B' $$

    where $n$ is the cooperativity (number of binders).

    Parameters
    ----------
    name : str, default='one_step_cooperative_binding'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='binding'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/binding_parameters.tsv'
        Path to file containing default parameter values for binding
        mechanisms.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('binding').
    parameter_database : ParameterDatabase
        Database storing default parameters for this mechanism.

    See Also
    --------
    Two_Step_Cooperative_Binding : Sequential cooperative binding mechanism.
    Combinatorial_Cooperative_Binding : Multiple distinct binders binding.
    One_Step_Binding : Simple binding without cooperativity.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is used to model cooperative binding where multiple
    identical ligands bind simultaneously to a receptor. Common examples
    include:

    - Oxygen binding to hemoglobin
    - Transcription factor dimers binding to DNA
    - Cooperative enzyme-substrate interactions

    The mechanism generates a single reversible mass-action reaction with
    forward rate constant 'kb' and reverse rate constant 'ku'.  The
    'cooperativity' parameter determines the stoichiometry of the binding
    reaction. A cooperativity of 2 models dimeric binding, 3 for trimeric,
    etc.

    Required parameters for this mechanism:

    - 'kb' : Forward binding rate constant
    - 'ku' : Reverse unbinding rate constant
    - 'cooperativity' : Number of binder molecules that bind simultaneously

    Examples
    --------
    Create a mechanism for dimeric transcription factor binding:

    >>> promoter = bcp.RegulatedPromoter(
    ...     name='p_dimer',
    ...     regulators='TF_dimer',
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[promoter],
    ...     mechanisms={'binding': bcp.One_Step_Cooperative_Binding()},
    ...     parameters={'cooperativity': 2, 'kb': 0.1, 'ku': 0.01}
    ... )

    """

    def __init__(
        self,
        name='one_step_cooperative_binding',
        mechanism_type='binding',
        parameter_file='mechanisms/binding_parameters.tsv',
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        binder,
        bindee,
        complex_species=None,
        cooperativity=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for cooperative binding.

        Creates the species involved in cooperative binding: the binder,
        bindee, and the resulting complex containing multiple binders bound
        to the bindee.

        Parameters
        ----------
        binder : Species
            The ligand species that binds cooperatively.
        bindee : Species
            The target species being bound to.
        complex_species : Species, optional
            Pre-specified complex species. If None, automatically creates a
            Complex containing $n$ binders and 1 bindee, where $n$ is the
            cooperativity.
        cooperativity : int or float, optional
            Number of binder molecules that bind simultaneously. If None,
            retrieved from component parameters using part_id.
        component : Component, optional
            Component containing parameter values. Required if cooperativity
            is not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            'repr(binder)-repr(bindee)'.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [binder, bindee, complex] where complex is the
            cooperative binding product.

        Raises
        ------
        ValueError
            If neither component nor cooperativity is provided.
        TypeError
            If complex_species is not a Species or None.

        Notes
        -----
        The `cooperativity` parameter determines how many binder molecules
        are incorporated into the complex. For example, cooperativity=2
        creates a complex with 2 binders and 1 bindee.

        """
        if part_id is None:
            part_id = repr(binder) + '-' + repr(bindee)

        if cooperativity is None and component is not None:
            cooperativity = component.get_parameter(
                'cooperativity',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        elif component is None and cooperativity is None:
            raise ValueError(
                "Must pass in a Component or values for cooperativity"
            )

        complexS = None
        if isinstance(complex_species, Species):
            complexS = complex_species
        elif complex_species is not None:
            raise TypeError(
                "complex_species keyword must be a Species, or None."
            )

        if complexS is None:
            complexS = Complex([binder] * int(cooperativity) + [bindee])

        return [binder, bindee, complexS]

    def update_reactions(
        self,
        binder,
        bindee,
        complex_species=None,
        component=None,
        kb=None,
        ku=None,
        part_id=None,
        cooperativity=None,
        **kwargs,
    ):
        """Generate reactions for cooperative binding.

        Creates a single reversible mass-action reaction for the cooperative
        binding of multiple binder molecules to a bindee.

        Parameters
        ----------
        binder : Species
            The ligand species that binds cooperatively.
        bindee : Species
            The target species being bound to.
        complex_species : Species, optional
            Pre-specified complex species. If None, automatically creates a
            Complex containing $n$ binders and 1 bindee.
        component : Component, optional
            Component containing parameter values. Required if kb, ku, or
            cooperativity are not provided directly.
        kb : Parameter or float, optional
            Forward binding rate constant. If None, retrieved from component
            parameters.
        ku : Parameter or float, optional
            Reverse unbinding rate constant. If None, retrieved from
            component parameters.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            'repr(binder)-repr(bindee)'.
        cooperativity : int or float, optional
            Number of binder molecules that bind simultaneously. If None,
            retrieved from component parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single reversible mass-action reaction for
            cooperative binding.

        Raises
        ------
        ValueError
            If component is None and any of kb, ku, or cooperativity is not
            provided.
        TypeError
            If complex_species is not a Species or None.

        Notes
        -----
        The reaction stoichiometry is determined by the `cooperativity`
        parameter:

        - cooperativity * binder + bindee <--> complex

        The forward rate is `kb` and reverse rate is `ku`. The reaction
        follows mass-action kinetics with the appropriate stoichiometric
        coefficients.

        """
        if part_id is None:
            part_id = repr(binder) + '-' + repr(bindee)
        if kb is None and component is not None:
            kb = component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            )
        if ku is None and component is not None:
            ku = component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            )
        if cooperativity is None and component is not None:
            cooperativity = component.get_parameter(
                'cooperativity',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        if component is None and (
            kb is None or ku is None or cooperativity is None
        ):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and "
                "coopertivity."
            )

        complexS = None
        if complex_species is None:
            complex_name = None
        elif isinstance(complex_species, Species):
            complexS = complex_species
        else:
            raise TypeError(
                "complex_species keyword must be a str, Species, or None."
            )

        if complexS is None:
            complexS = Complex(
                [binder] * int(cooperativity) + [bindee], name=complex_name
            )

        inputs = [
            WeightedSpecies(species=binder, stoichiometry=cooperativity),
            WeightedSpecies(species=bindee, stoichiometry=1),
        ]

        rxns = [
            Reaction.from_massaction(
                inputs=inputs, outputs=[complexS], k_forward=kb, k_reverse=ku
            )
        ]
        return rxns


class Two_Step_Cooperative_Binding(Mechanism):
    r"""Sequential cooperative binding mechanism with oligomerization.

    A 'binding' mechanism where multiple binder molecules first oligomerize,
    then the oligomer binds to the bindee in a two-step process. This models
    cooperative binding where ligands must first form a multimeric complex
    before binding to their target.

    The binding process follows two sequential reactions:

    1. $n A <--> A_n$ (oligomerization)
    2. $A_n + B <--> A_n:B$ (binding)

    where $n$ is the cooperativity.

    Parameters
    ----------
    name : str, default='two_step_cooperative_binding'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='binding'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/binding_parameters.tsv'
        Path to file containing default parameter values for binding
        mechanisms.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('binding').
    parameter_database : ParameterDatabase
        Database storing default parameters for this mechanism.

    See Also
    --------
    One_Step_Cooperative_Binding : Single-step cooperative binding.
    Combinatorial_Cooperative_Binding : Multiple distinct binders.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models cooperative binding as a two-step process:

    1. Oligomerization: Binder molecules first associate to form an
       oligomer (dimer, trimer, etc.)
    2. Binding: The oligomer then binds to the target

    This is useful for modeling:

    - Protein dimerization followed by DNA binding
    - Receptor oligomerization before ligand binding
    - Sequential assembly and binding processes

    Required parameters for this mechanism:

    - 'kb1' : Forward rate constant for oligomerization
    - 'ku1' : Reverse rate constant for oligomerization
    - 'kb2' : Forward rate constant for oligomer-bindee binding
    - 'ku2' : Reverse rate constant for oligomer-bindee binding
    - 'cooperativity' : Number of binder molecules in the oligomer

    Examples
    --------
    Model transcription factor dimerization followed by DNA binding:

    >>> mech = bcp.Two_Step_Cooperative_Binding()
    >>> # TF dimerizes (2*TF <-> TF2), then binds DNA (TF2 + DNA <-> TF2:DNA)
    >>> params = {
    ...     'cooperativity': 2,
    ...     'kb1': 0.1, 'ku1': 0.01,  # Dimerization rates
    ...     'kb2': 1.0, 'ku2': 0.001   # DNA binding rates
    ... }

    Model trimeric receptor assembly and activation:

    >>> mech = bcp.Two_Step_Cooperative_Binding()
    >>> params = {
    ...     'cooperativity': 3,  # Trimeric receptor
    ...     'kb1': 0.05, 'ku1': 0.1,   # Trimerization
    ...     'kb2': 10.0, 'ku2': 0.01   # Ligand binding
    ... }

    """

    def __init__(
        self,
        name='two_step_cooperative_binding',
        mechanism_type='binding',
        parameter_file='mechanisms/binding_parameters.tsv',
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        binder,
        bindee,
        component=None,
        complex_species=None,
        n_mer_species=None,
        cooperativity=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for two-step cooperative binding.

        Creates the species involved in sequential cooperative binding:
        binder, bindee, oligomer (n-mer), and final complex.

        Parameters
        ----------
        binder : Species
            The ligand species that oligomerizes then binds.
        bindee : Species
            The target species that the oligomer binds to.
        component : Component, optional
            Component containing parameter values. Required if cooperativity
            is not provided directly.
        complex_species : Species, optional
            Pre-specified final complex species. If None, automatically
            creates a Complex containing the n-mer and bindee.
        n_mer_species : Species, optional
            Pre-specified oligomer species. If None, automatically creates a
            Complex containing $n$ binders.
        cooperativity : int or float, optional
            Number of binders in the oligomer. If None, retrieved from
            component parameters.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            'repr(binder)-repr(bindee)'.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [binder, bindee, complex, n_mer] where n_mer is
            the oligomer and complex is the final bound product.

        Raises
        ------
        ValueError
            If neither component nor cooperativity is provided.
        TypeError
            If n_mer_species or complex_species is not a Species or None.

        Notes
        -----
        The n_mer represents the oligomerized form of the binder (e.g., a
        dimer for cooperativity=2). The complex represents the n_mer bound
        to the bindee.

        """
        if part_id is None:
            part_id = repr(binder) + '-' + repr(bindee)

        if cooperativity is None and component is not None:
            cooperativity = component.get_parameter(
                'cooperativity',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        elif component is None and cooperativity is None:
            raise ValueError(
                "Must pass in a Component or values for cooperativity"
            )

        n_mer = None
        if isinstance(n_mer_species, Species):
            n_mer = n_mer_species
        elif n_mer_species is not None:
            raise TypeError(
                "n_mer_species keyword nust be a Species, or None. Not "
                + str(n_mer_species)
            )

        if n_mer is None:
            n_mer = Complex(cooperativity * [binder])

        complexS = None
        if isinstance(complex_species, Species):
            complexS = complex_species
        elif complex_species is not None:
            raise TypeError(
                "complex_species keyword must be a Species, or None. Not "
                + str(complex_species)
            )

        if complexS is None:
            complexS = Complex([n_mer, bindee])
        return [binder, bindee, complexS, n_mer]

    def update_reactions(
        self,
        binder,
        bindee,
        kb=None,
        ku=None,
        component=None,
        part_id=None,
        cooperativity=None,
        complex_species=None,
        n_mer_species=None,
        **kwargs,
    ):
        """Generate reactions for two-step cooperative binding.

        Creates two sequential reactions: oligomerization of binders followed
        by oligomer binding to the bindee.

        Parameters
        ----------
        binder : Species
            The ligand species that oligomerizes then binds.
        bindee : Species
            The target species that the oligomer binds to.
        kb : list of float or Parameter, optional
            Forward rate constants [kb1, kb2] for oligomerization and binding.
            If None, retrieved from component parameters.
        ku : list of float or Parameter, optional
            Reverse rate constants [ku1, ku2] for oligomerization and binding.
            If None, retrieved from component parameters.
        component : Component, optional
            Component containing parameter values. Required if kb, ku, or
            cooperativity are not provided.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            'repr(binder)-repr(bindee)'.
        cooperativity : int or float, optional
            Number of binders in the oligomer. If None, retrieved from
            component parameters.
        complex_species : Species, optional
            Pre-specified final complex species.
        n_mer_species : Species, optional
            Pre-specified oligomer species.
        **kwargs
            Additional keyword arguments passed to update_species.

        Returns
        -------
        list of Reaction
            List containing two reactions:

            1. Oligomerization: cooperativity*binder <--> n_mer
            2. Binding: n_mer + bindee <--> complex

        Raises
        ------
        ValueError
            If component is None and `kb`, `ku`, or `cooperativity` is not
            provided, or if `kb` and `ku` do not contain exactly 2 values
            each.

        Notes
        -----
        The two-step process uses separate rate constants:

        - 'kb1', 'ku1': Control oligomerization kinetics
        - 'kb2', 'ku2': Control oligomer-bindee binding kinetics

        This separation allows modeling of processes where oligomerization
        and binding have different kinetic properties.

        """
        if part_id is None:
            repr(binder) + '-' + repr(bindee)
        if (
            kb is None or ku is None or cooperativity is None
        ) and component is not None:
            kb1 = component.get_parameter(
                'kb1', part_id=part_id, mechanism=self
            )
            kb2 = component.get_parameter(
                'kb2', part_id=part_id, mechanism=self
            )
            ku1 = component.get_parameter(
                'ku1', part_id=part_id, mechanism=self
            )
            ku2 = component.get_parameter(
                'ku2', part_id=part_id, mechanism=self
            )
            cooperativity = component.get_parameter(
                'cooperativity',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        elif component is None and (
            kb is None or ku is None or cooperativity is None
        ):
            raise ValueError(
                "Must pass in a Component or values for kb, ku, and "
                "cooperativity"
            )
        elif len(kb) != len(ku) != 2:
            raise ValueError(
                "kb and ku must contain 2 values each for 'two-step binding'"
            )
        else:
            kb1, kb2 = kb
            ku1, ku2 = ku

        binder, bindee, complexS, n_mer = self.update_species(
            binder,
            bindee,
            component=component,
            complex_species=complex_species,
            n_mer_species=n_mer_species,
            cooperativity=cooperativity,
            part_id=part_id,
            **kwargs,
        )

        inputs_for_rxn1 = [
            WeightedSpecies(species=binder, stoichiometry=cooperativity)
        ]
        rxns = [
            Reaction.from_massaction(
                inputs=inputs_for_rxn1,
                outputs=[n_mer],
                k_forward=kb1,
                k_reverse=ku1,
            ),
            Reaction.from_massaction(
                inputs=[n_mer, bindee],
                outputs=[complexS],
                k_forward=kb2,
                k_reverse=ku2,
            ),
        ]

        return rxns


class Combinatorial_Cooperative_Binding(Mechanism):
    """Combinatorial binding mechanism for multiple distinct ligands.

    A 'binding' mechanism where different types of binder molecules can bind
    to a bindee in various combinations, each with its own cooperativity. This
    models complex regulatory scenarios where multiple transcription factors
    or ligands can bind to the same target in different combinations, each
    producing a distinct complex.

    The mechanism generates all possible binding combinations and the
    reactions between them, considering individual binding affinities and
    cooperativities for each binder type.

    Parameters
    ----------
    name : str, default='Combinatorial_Cooperative_binding'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='binding'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/binding_parameters.tsv'
        Path to file containing default parameter values for binding
        mechanisms.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('binding').
    parameter_database : ParameterDatabase
        Database storing default parameters for this mechanism.

    See Also
    --------
    One_Step_Cooperative_Binding : Single binder type cooperative binding.
    CombinatorialPromoter : Component that uses this mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is designed for modeling complex regulatory logic where:

    - Multiple different regulators can bind to the same target
    - Each regulator can have its own cooperativity (e.g., some bind as
      dimers, others as monomers)
    - All possible combinations of bound states are generated
    - Each transition between states has specific rate constants

    The mechanism generates a complete reaction network connecting all
    possible bound states. For $n$ different binders, this creates $2^n - 1$
    different complexes (excluding the unbound state).

    Required parameters for this mechanism (per binder):

    - 'kb': Forward binding rate
    - 'ku': Reverse unbinding rate
    - 'cooperativity': Number of molecules binding together

    This is commonly used for:

    - Complex promoter regulation with multiple transcription factors
    - Multi-ligand receptor systems
    - Combinatorial protein complex assembly

    Examples
    --------
    Model a promoter with two different transcription factors:

    >>> A, B = bcp.Species('A'), bcp.Species('B')
    >>> AND_promoter = bcp.CombinatorialPromoter(
    ...     'AND_promoter', [A, B], tx_capable_list=[[A, B]], leak=False)
    ... AND_assembly = bcp.DNAassembly(
    ...     'AND', promoter=AND_promoter, rbs='medium', protein='GFP')
    ... mixture = bcp.ExpressionExtract(
    ...     name='AND_mixture', components=[AND_assembly],
    ...     parameter_file=[
    ...         'mechanisms/binding_parameters.tsv',
    ...         'mixtures/extract_parameters.tsv',
    ...     ]
    ... )
    ... crn = mixture.compile_crn()

    """

    def __init__(
        self,
        name='Combinatorial_Cooperative_binding',
        mechanism_type='binding',
        parameter_file='mechanisms/binding_parameters.tsv',
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def make_cooperative_complex(self, combo, bindee, cooperativity):
        """Create a complex with multiple cooperative binders.

        Constructs a complex species containing the specified combination of
        binders (each repeated according to its cooperativity) bound to the
        bindee.

        Parameters
        ----------
        combo : tuple or list of Species
            Combination of binder species to include in the complex.
        bindee : Species
            The target species being bound to.
        cooperativity : dict
            Dictionary mapping binder names (str) to their cooperativity
            values (int). Determines how many copies of each binder are
            included.

        Returns
        -------
        Species or Complex
            If only bindee is present (empty combo), returns bindee alone.
            Otherwise returns a Complex containing all binders (repeated per
            cooperativity) and the bindee.

        Notes
        -----
        For each binder in combo, the method adds cooperativity[binder.name]
        copies to the complex. For example, if binder A has cooperativity 2
        and binder B has cooperativity 1, the complex for combo=[A, B] would
        contain [A, A, B, bindee].

        """
        complexed_species_list = []
        for binder in combo:
            binder_cooperativity = int(cooperativity[binder.name])
            # I hope that cooperativity is an int! what if it isn't
            if binder_cooperativity > 1:
                complexed_species_list += [binder] * binder_cooperativity
            else:
                complexed_species_list += [binder]
        complexed_species_list += [bindee]
        if len(complexed_species_list) == 1:
            myspecies = complexed_species_list[0]
        else:
            myspecies = Complex(complexed_species_list)
        return myspecies

    def update_species(
        self,
        binders,
        bindee,
        cooperativity=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate all species for combinatorial binding.

        Creates all possible complexes from combinations of binders bound to
        the bindee, considering each binder's cooperativity.

        Parameters
        ----------
        binders : list of Species
            List of different binder species that can bind in combinations.
        bindee : Species
            The target species being bound to.
        cooperativity : dict, optional
            Dictionary mapping binder names to cooperativity values. If None
            for any binder, retrieved from component parameters.
        component : Component, optional
            Component containing parameter values. Required if cooperativity
            values are not provided.
        part_id : str, optional
            Base identifier for parameter lookup. Individual binder parameters
            are looked up as 'part_id_bindername'.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List of all possible complexes from binding combinations. For n
            binders, returns $2^n - 1$ complexes (all combinations except
            unbound bindee).

        Raises
        ------
        ValueError
            If component is None and cooperativity is not provided for all
            binders.

        Notes
        -----
        This method generates all possible combinations of binders:
        - Single binders: A:bindee, B:bindee, etc.
        - Pairs: A:B:bindee, A:C:bindee, etc.
        - Higher combinations up to all binders bound simultaneously

        Each complex respects the individual cooperativity of its binders.

        """
        cooperativity_dict = {}
        out_species = []

        # Figure out part_id (if list, assume first element is primary ID)
        if part_id is None:
            prefix = ''
            suffix_list = []
        elif isinstance(part_id, list):
            prefix = part_id[0] + '_'
            suffix_list = part_id[1:]
        else:
            prefix = part_id + '_'
            suffix_list = []

        for binder in binders:
            binder_partid = prefix + binder.name
            if (
                (cooperativity is None)
                or (
                    isinstance(cooperativity, dict)
                    and binder_partid not in cooperativity
                )
                and (component is not None)
            ):
                # here we are extracting the relevant cooperativity value from
                # the dictionary which should be passed in as the
                # cooperativity argument
                coop_val = component.get_parameter(
                    'cooperativity',
                    part_id=[binder_partid] + suffix_list,
                    mechanism=self,
                    return_numerical=True,
                )
            elif (
                isinstance(cooperativity, dict)
                and binder_partid in cooperativity
            ):
                coop_val = cooperativity[binder_partid]
            if component is None and (cooperativity is None):
                raise ValueError(
                    "Must pass in a Component or values for kb, ku, and "
                    "coopertivity."
                )
            cooperativity_dict[binder.name] = coop_val

        for i in range(1, len(binders) + 1):
            for combo in it.combinations(binders, i):
                # go through every possible combination of reactants and dna
                # and make all the complexes
                out_species += [
                    self.make_cooperative_complex(
                        combo, bindee, cooperativity_dict
                    )
                ]
        return out_species

    def update_reactions(
        self,
        binders,
        bindee,
        component=None,
        kbs=None,
        kus=None,
        part_id=None,
        cooperativity=None,
        **kwargs,
    ):
        """Generate reactions for all combinatorial binding transitions.

        Creates reactions connecting all possible binding states, where each
        reaction represents one binder type associating or dissociating while
        other binders remain bound.

        Parameters
        ----------
        binders : list of Species
            List of different binder species that can bind in combinations.
        bindee : Species
            The target species being bound to.
        component : Component, optional
            Component containing parameter values. Required if rate constants
            or cooperativities are not provided.
        kbs : dict, optional
            Dictionary mapping binder names to forward rate constants (kb).
            If None for any binder, retrieved from component parameters.
        kus : dict, optional
            Dictionary mapping binder names to reverse rate constants (ku).
            If None for any binder, retrieved from component parameters.
        part_id : str, optional
            Base identifier for parameter lookup. Individual binder parameters
            are looked up as 'part_id_bindername'.
        cooperativity : dict, optional
            Dictionary mapping binder names to cooperativity values. If None
            for any binder, retrieved from component parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of all reactions connecting binding states. Each reaction
            represents adding or removing one binder type to/from an existing
            complex.

        Raises
        ------
        ValueError
            If component is None and any required parameters (kb, ku,
            cooperativity) are not provided.

        Notes
        -----
        The reaction network connects all possible binding states such that:

        - Each reaction adds or removes exactly one binder type
        - Rate constants are specific to each binder
        - Cooperativity determines the stoichiometry of each binder

        For $n$ binders, this generates approximately $n 2^(n-1)$ reactions,
        connecting the $2^n$ possible states (including unbound).

        Each binder requires three parameters:

        - 'kb': Forward binding rate constant
        - 'ku': Reverse unbinding rate constant
        - 'cooperativity': Stoichiometry of that binder

        The algorithm avoids generating duplicate reactions by tracking which
        transitions have been created.

        """
        binder_params = {}

        # Figure out part_id (if list, assume first element is primary ID)
        if part_id is None:
            prefix = ''
            suffix_list = []
        elif isinstance(part_id, list):
            prefix = part_id[0] + '_'
            suffix_list = part_id[1:]
        else:
            prefix = part_id + '_'
            suffix_list = []

        for binder in binders:
            binder_partid = prefix + binder.name
            if (isinstance(kbs, dict) and binder not in kbs) or (
                not isinstance(kbs, dict) and component is not None
            ):
                kb = component.get_parameter(
                    'kb',
                    part_id=[binder_partid] + suffix_list,
                    mechanism=self,
                )
            elif isinstance(kbs, dict) and binder in kbs:
                kb = kbs[binder.name]
            elif not isinstance(kbs, dict) and component is None:
                raise ValueError(
                    "Must pass in a Component or values for kb, ku, and "
                    "cooperativity."
                )
            if (isinstance(kus, dict) and binder not in kus) or (
                kus is None and component is not None
            ):
                ku = component.get_parameter(
                    'ku',
                    part_id=[binder_partid] + suffix_list,
                    mechanism=self,
                )
            elif isinstance(kus, dict) and binder in kus:
                ku = kus[binder.name]
            elif not isinstance(kus, dict) and component is None:
                raise ValueError(
                    "Must pass in a Component or values for kb, ku, and "
                    "cooperativity."
                )
            if (
                (cooperativity is None)
                or (
                    isinstance(cooperativity, dict)
                    and binder.name not in cooperativity
                )
                and component is not None
            ):
                coop_val = component.get_parameter(
                    'cooperativity',
                    part_id=[binder_partid] + suffix_list,
                    mechanism=self,
                    return_numerical=True,
                )
            elif (
                isinstance(cooperativity, dict)
                and binder.name in cooperativity
            ):
                coop_val = cooperativity[binder.name]
            if component is None and (
                kb is None or ku is None or cooperativity is None
            ):
                raise ValueError(
                    "Must pass in a Component or values for kb, ku, and "
                    "cooperativity."
                )
            binder_params[binder] = {
                'kb': kb,
                'ku': ku,
                'cooperativity': coop_val,
            }
        # out_rxns = []
        rxndict = {}
        coop_dict = {
            a.name: binder_params[a]['cooperativity'] for a in binder_params
        }
        for i in range(1, len(binders) + 1):
            for combo in it.combinations(binders, i):
                # come up with all combinations of binders

                product = self.make_cooperative_complex(
                    combo, bindee, coop_dict
                )
                # this is the complex which becomes the product
                for binder in combo:
                    # then in each case, one binder can be added to make
                    # this combination
                    reactant = tuple(set(combo) - set([binder]))
                    rxn_prototype = (binder, reactant)

                    # this part makes a describer of the reaction; which
                    # reactants are combining?
                    # print(rxn_prototype)
                    # print(rxndict)
                    if rxn_prototype in rxndict:
                        # if we already did this reaction then forget about it
                        continue
                    else:
                        reactant_complex = self.make_cooperative_complex(
                            reactant, bindee, coop_dict
                        )

                        inputs = [
                            WeightedSpecies(
                                species=binder,
                                stoichiometry=binder_params[binder][
                                    'cooperativity'
                                ],
                            ),
                            WeightedSpecies(
                                species=reactant_complex, stoichiometry=1
                            ),
                        ]

                        reaction = Reaction.from_massaction(
                            inputs=inputs,
                            outputs=[product],
                            k_forward=binder_params[binder]['kb'],
                            k_reverse=binder_params[binder]['ku'],
                        )
                        rxndict[rxn_prototype] = reaction
        return [rxndict[a] for a in rxndict]


class One_Step_Binding(Mechanism):
    r"""Simple binding mechanism for multiple species without cooperativity.

    A 'binding' mechanism to model the simultaneous binding of multiple
    species into a single complex in one concerted step. Unlike cooperative
    binding mechanisms, this treats all species equally without cooperativity
    factors - each species contributes exactly one molecule to the complex.

    The binding reaction follows:

      $$ 'S'_1 + 'S'_2 + ... + 'S'_n <--> 'S'_1:'S'_2:...:'S'_n $$

    Parameters
    ----------
    name : str, default='one_step_binding'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='binding'
        Type classification of this mechanism.
    parameter_file : str, default='mechanisms/binding_parameters.tsv'
        Path to file containing default parameter values for binding
        mechanisms.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('binding').
    parameter_database : ParameterDatabase
        Database storing default parameters for this mechanism.

    See Also
    --------
    One_Step_Cooperative_Binding : Binding with cooperativity.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is the simplest binding model where multiple distinct
    species come together to form a complex. Each species contributes
    exactly one molecule (stoichiometry of 1) to the complex.

    Common applications include:

    - Protein complex formation from distinct subunits
    - Multi-component enzyme assembly
    - Simple receptor-ligand binding
    - Formation of heterodimers or heterotrimers

    The mechanism generates a single reversible mass-action reaction with
    rate constants 'kb' (forward) and 'ku' (reverse).

    Key differences from cooperative binding:

    - No 'cooperativity' parameter - all stoichiometries are 1
    - Can handle arbitrary lists of different species
    - Simpler parameter structure (single 'kb', 'ku' for the entire reaction)

    Required parameters for this mechanism:

    - 'kb' : Forward binding rate constant
    - 'ku' : Reverse unbinding rate constant

    Examples
    --------
    Model receptor-ligand binding:

    >>> mech = bcp.One_Step_Binding()
    >>> ligand, receptor = bcp.Species('L'), bcp.Species('R')
    >>> mech.update_species(
    ...     binder=ligand,
    ...     bindee=receptor,
    ...     kb=1.0, ku=0.001
    ... )
    [L, R, complex_L_R_]

    Model formation of a three-protein complex:

    >>> proteins = [
    ...     bcp.Species(s, material_type='protein') for s in ['A', 'B', 'C']]
    >>> rxns = mech.update_reactions(
    ...     binder=proteins[0],
    ...     bindee=proteins[1:],
    ...     kb=0.1, ku=0.01
    ... )
    >>> # Generates: A + B + C <--> A:B:C

    """

    def __init__(
        self,
        name='one_step_binding',
        mechanism_type='binding',
        parameter_file='mechanisms/binding_parameters.tsv',
    ):
        Mechanism.__init__(
            self, name, mechanism_type, parameter_file=parameter_file
        )

    def update_species(
        self,
        binder,
        bindee,
        component=None,
        complex_species=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for simple multi-species binding.

        Creates the list of species involved in the binding reaction: all
        input species plus the resulting complex.

        Parameters
        ----------
        binder : Species or list of Species
            The first species or list of species to bind. Automatically
            converted to list if single species provided.
        bindee : Species or list of Species
            The second species or list of species to bind. Automatically
            converted to list if single species provided.
        component : Component, optional
            Component containing this mechanism (unused but kept for API
            consistency).
        complex_species : Species, optional
            Pre-specified complex species. If None, automatically creates a
            Complex containing all binders and bindees.
        part_id : str, optional
            Identifier for parameter lookup. If None, automatically generated
            from species names as `name1_name2_..._nameN`.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing all input species plus the complex. Format:
            [binder1, ..., bindee1, ..., complex].

        Notes
        -----
        The binder/bindee distinction is primarily for API consistency with
        other binding mechanisms. Functionally, all species are treated
        equally in the binding reaction.

        The complex is created as a Complex object containing all input
        species in the order: binders + bindees.

        """
        if not isinstance(binder, list):
            binder = [binder]
        if not isinstance(bindee, list):
            bindee = [bindee]
        species = binder + bindee
        if part_id is None:
            part_id = ''
            for s in species:
                part_id += s.name + '_'
            part_id = part_id[:-1]

        if complex_species is None:
            complex_species = Complex(species)

        return species + [complex_species]

    def update_reactions(
        self,
        binder,
        bindee,
        component=None,
        complex_species=None,
        part_id=None,
        kb=None,
        ku=None,
        **kwargs,
    ):
        """Generate reaction for simple multi-species binding.

        Creates a single reversible mass-action reaction for the binding of
        all species into a complex.

        Parameters
        ----------
        binder : Species or list of Species
            The first species or list of species to bind. Automatically
            converted to list if single species provided.
        bindee : Species or list of Species
            The second species or list of species to bind. Automatically
            converted to list if single species provided.
        component : Component, optional
            Component containing parameter values. Required if kb or ku are
            not provided directly.
        complex_species : Species, optional
            Pre-specified complex species. If None, automatically creates a
            Complex containing all binders and bindees.
        part_id : str, optional
            Identifier for parameter lookup. If None, automatically generated
            from species names as `name1_name2_..._nameN`.
        kb : Parameter or float, optional
            Forward binding rate constant. If None, retrieved from component
            parameters.
        ku : Parameter or float, optional
            Reverse unbinding rate constant. If None, retrieved from
            component parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single reversible mass-action reaction for the
            binding of all species.

        Raises
        ------
        ValueError
            If component is None and kb or ku is not provided.

        Notes
        -----
        The reaction has equal stoichiometry (1) for all species:
        species1 + species2 + ... + speciesN <--> complex

        This is simpler than cooperative binding mechanisms which can have
        varying stoichiometries for different species.

        """
        if not isinstance(binder, list):
            binder = [binder]
        if not isinstance(bindee, list):
            bindee = [bindee]
        species = binder + bindee
        if part_id is None:
            part_id = ''
            for s in species:
                part_id += s.name + '_'
            part_id = part_id[:-1]

        if (kb is None or ku is None) and component is not None:
            kb = component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            )
            ku = component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            )
        elif component is None and (kb is None or ku is None):
            raise ValueError(
                "Must pass in a Component or values for kb and ku"
            )

        if complex_species is None:
            complex_species = Complex(species)

        return [
            Reaction.from_massaction(
                inputs=species,
                outputs=[complex_species],
                k_forward=kb,
                k_reverse=ku,
            )
        ]
