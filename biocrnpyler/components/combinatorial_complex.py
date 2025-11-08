#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


from itertools import permutations

from ..core.component import Component
from ..core.species import Complex, ComplexSpecies


class CombinatorialComplex(Component):
    """Complex formed through combinatorial binding of multiple species.

    A `CombinatorialComplex` component represents a complex that can form
    through multiple combinatorial binding pathways. The component
    enumerates all possible intermediate complexes and generates binding
    reactions between initial states and final states, optionally
    constrained by intermediate states and excluded states. Uses a
    'binding' mechanism to generate combinatorial binding reactions.

    Parameters
    ----------
    final_states : ComplexSpecies or list of ComplexSpecies
        The final complex(es) to be formed. All binding reactions
        ultimately lead to these states.
    initial_states : list of Species or ComplexSpecies, optional
        Starting species that bind together to form final_states. If None,
        defaults to all individual species contained within final_states.
    intermediate_states : list of ComplexSpecies, optional
        Allowed intermediate complexes formed during binding. Restricts
        the binding pathway. If None, all possible intermediates are
        enumerated.
    excluded_states : list of Species or ComplexSpecies, optional
        Species or complexes that are NOT allowed to form. If None, no
        complexes are excluded.
    name : str, optional
        Name of the component. If None, automatically generated from
        final_states names.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    final_states : list of ComplexSpecies
        List of final complex states.
    initial_states : list of Species or ComplexSpecies
        List of initial binding species.
    intermediate_states : list of ComplexSpecies or None
        List of allowed intermediate complexes, or None if unrestricted.
    excluded_states : list
        List of excluded species/complexes.
    sub_species : list of Species
        All individual species contained in final_states.
    combination_dict : dict
        Dictionary storing computed binding combinations.

    See Also
    --------
    ChemicalComplex : Simple complex of two or more molecules.
    Component : Base class for biomolecular components.
    ComplexSpecies : Species subclass for molecular complexes.

    Notes
    -----
    The combinatorial binding process generates reactions based on the
    provided constraints:

    Case 1 - only `final_states` given: all species in final_states bind
    combinatorially:

        individual_species <--> all_intermediates <--> final_states

    Case 2 - `final_states` + `initial_states`: binding starts from
    specified initial states directly to final states:

        initial_states <--> final_states

    Case 3 - `final_states` + `intermediate_states`: binding restricted to
    specified intermediates:

        individual_species <--> intermediate_states <--> final_states

    Case 4 - `final_states` + `initial_states` + `intermediate_states`: both
    initial and intermediate constraints applied:

        initial_states <--> intermediate_states <--> final_states

    The component name is automatically generated as a concatenation of
    final_states names separated by underscores if not provided.

    Examples
    --------
    Example 1: Full combinatorial binding

    >>> A = bcp.Species('A')
    >>> B = bcp.Species('B')
    >>> C = bcp.Species('C')
    >>> final = bcp.Complex([A, B, C])
    >>> cc = bcp.CombinatorialComplex(
    ...     final_states=final,
    ...     mechanisms={'binding': bcp.One_Step_Binding()},
    ...     parameters={'kb': 1e-1, 'ku': 1e-1})

    Initial states default to [A, B, C]. All intermediates
    [A, B], [A, C], [B, C] are enumerated, resulting in 6 reversible
    reactions:

    1. A + B <--> Complex([A, B])
    2. A + C <--> Complex([A, C])
    3. B + C <--> Complex([B, C])
    4. Complex([A, B]) + C <--> Complex([A, B, C])
    5. Complex([A, C]) + B <--> Complex([A, B, C])
    6. Complex([B, C]) + A <--> Complex([A, B, C])

    Example 2: Constrained initial states

    >>> initial = [bcp.Complex([A, B]), bcp.Complex([A, C])]
    >>> cc = bcp.CombinatorialComplex(
    ...     final_states=final, initial_states=initial,
    ...     mechanisms={'binding': bcp.One_Step_Binding()},
    ...     parameters={'kb': 1e-1, 'ku': 1e-1})

    Results in 2 reactions:

    1. Complex([A, B]) + C <--> Complex([A, B, C])
    2. Complex([A, C]) + B <--> Complex([A, B, C])

    Example 3: Restricted intermediate states

    >>> inter = [bcp.Complex([A, B]), bcp.Complex([A, C])]
    >>> cc = bcp.CombinatorialComplex(
    ...     final_states=final, intermediate_states=inter,
    ...     mechanisms={'binding': bcp.One_Step_Binding()},
    ...     parameters={'kb': 1e-1, 'ku': 1e-1})

    Results in 4 reactions:

    1. A + B <--> Complex([A, B])
    2. A + C <--> Complex([A, C])
    3. Complex([A, B]) + C <--> Complex([A, B, C])
    4. Complex([A, C]) + B <--> Complex([A, B, C])

    Example 4: Multiple final states with homodimers

    >>> final = [bcp.Complex([A, A, B]), bcp.Complex([A, B, B])]
    >>> cc = bcp.CombinatorialComplex(
    ...     final_states=final,
    ...     mechanisms={'binding': bcp.One_Step_Binding()},
    ...     parameters={'kb': 1e-1, 'ku': 1e-1})

    Results in 7 reactions including homodimer formation:

    1. A + A <--> Complex([A, A])
    2. Complex([A, A]) + B <--> Complex([A, A, B])
    3. B + B <--> Complex([B, B])
    4. Complex([B, B]) + A <--> Complex([A, B, B])
    5. A + B <--> Complex([A, B])
    6. Complex([A, B]) + A <--> Complex([A, A, B])
    7. Complex([A, B]) + B <--> Complex([A, B, B])

    """

    def __init__(
        self,
        final_states,
        initial_states=None,
        intermediate_states=None,
        excluded_states=None,
        name=None,
        **kwargs,
    ):
        # The order these run in is important! (TODO: why?)
        self.final_states = final_states
        self.initial_states = initial_states
        self.intermediate_states = intermediate_states
        self.excluded_states = excluded_states

        # used to store combinations of species during update
        self.combination_dict = {}

        # Call super
        if name is None:
            name = ''
            for s in self.final_states:
                name += s.name + '_'
            name = name[:-1]
        super().__init__(name, **kwargs)

    # Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        """List of final complex states to be formed.

        Returns
        -------
        list of ComplexSpecies

        """
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        """Set the final complex states.

        Parameters
        ----------
        final_states : ComplexSpecies or list of ComplexSpecies
            Final complex(es) to be formed through combinatorial binding.

        Raises
        ------
        ValueError
            If any element in final_states is not a ComplexSpecies.

        Notes
        -----
        Also creates a list of all sub-species (individual species
        contained in the complexes) stored in `self.sub_species`.

        """
        final_states = self.set_species(final_states)
        if not isinstance(final_states, list):
            final_states = [final_states]

        # all final_states must be ComplexSpecies
        if not all([isinstance(s, ComplexSpecies) for s in final_states]):
            raise ValueError(
                f"final_states must be a list of {ComplexSpecies} (or "
                f"subclasses thereof). Recieved: {final_states}."
            )

        self._final_states = final_states

        # Then create a list of all sub-species included in
        # final_states Complexes
        self.sub_species = []
        for s in self.final_states:
            self.sub_species += s.species_set

    # Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        """List of initial states for binding.

        Returns
        -------
        list of Species or ComplexSpecies
        """
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        """Set the initial binding states.

        Parameters
        ----------
        initial_states : list of Species or ComplexSpecies, optional
            Starting species for combinatorial binding. If None, defaults
            to all individual species in final_states (sub_species).

        Raises
        ------
        ValueError
            If any initial state is not contained in sub_species or is not
            a ComplexSpecies made from sub_species.

        Notes
        -----
        Initial states must either be individual species from the
        final_states or ComplexSpecies composed of those species.

        """
        # set initial states
        if initial_states is None:
            self._initial_states = self.sub_species
        else:
            initial_states = self.set_species(initial_states)
            if not isinstance(initial_states, list):
                initial_states = [initial_states]

            for s in initial_states:
                if not (
                    s in self.sub_species
                    or (
                        isinstance(s, ComplexSpecies)
                        and all(
                            [ss in self.sub_species for ss in s.species_set]
                        )
                    )
                ):
                    raise ValueError(
                        f"Invalid initial species {s}; initial_states must "
                        "either be contained in the final_states or a "
                        "ComplexSpecies made of Species in the final_states."
                    )
            self._initial_states = initial_states

    # Intermediate states allows the user to restrict the complexes
    # formed between the intial state and final state
    @property
    def intermediate_states(self):
        """List of allowed intermediate complexes.

        Returns
        -------
        list of ComplexSpecies or None
        """
        return self._intermediate_states

    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        """Set the allowed intermediate complex states.

        Parameters
        ----------
        intermediate_states : list of ComplexSpecies, optional
            Allowed intermediate complexes formed between initial and
            final states. If None, all possible intermediates are
            enumerated.

        Raises
        ------
        ValueError
            If any intermediate state is not a ComplexSpecies, or if any
            contains species not in sub_species.

        Notes
        -----
        Restricting intermediate states limits the binding pathways and
        can reduce the number of reactions generated.

        """
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            intermediate_states = self.set_species(intermediate_states)
            if not isinstance(intermediate_states, list):
                intermediate_states = [intermediate_states]

            # All intermediate_states must be ComplexSpecies or
            # OrderdedComplexSpecies
            if not all(
                [isinstance(s, ComplexSpecies) for s in intermediate_states]
            ):
                raise ValueError(
                    f"intermediate must be a list of {ComplexSpecies} "
                    "(or subclasses thereof). "
                    f"Recieved: {intermediate_states}."
                )
            # All intermediate_states must be made of sub_species
            for s in intermediate_states:
                intermediate_sub_species = s.species_set
                if not all(
                    [
                        ss in self.sub_species
                        for ss in intermediate_sub_species
                    ]
                ):
                    raise ValueError(
                        f"intermediate species {s} contains subspecies not "
                        "in the final_states."
                    )

            self._intermediate_states = intermediate_states

    # Excluded states allows the user to exclude specific Species from
    # being enumerated
    @property
    def excluded_states(self):
        """list: Species or complexes excluded from enumeration."""
        return self._excluded_states

    @excluded_states.setter
    def excluded_states(self, excluded_states):
        """Set the excluded species and complexes.

        Parameters
        ----------
        excluded_states : list of Species or ComplexSpecies, optional
            Species or complexes that are NOT allowed to form. If None,
            no exclusions are applied (empty list).

        Notes
        -----
        Excluded states will not appear as reactants or products in
        generated reactions. Useful for preventing unwanted binding
        pathways.

        """
        if excluded_states is None:
            self._excluded_states = []
        else:
            excluded_states = self.set_species(excluded_states)
            if not isinstance(excluded_states, list):
                excluded_states = [excluded_states]
            self._excluded_states = excluded_states

    def compute_species_to_add(self, s0, sf):
        """Compute species needed to convert s0 into complex sf.

        Parameters
        ----------
        s0 : Species or ComplexSpecies
            Starting species or complex.
        sf : ComplexSpecies
            Target final complex.

        Returns
        -------
        list of Species or None
            List of species that need to be added to s0 to form sf. Returns
            None if s0 contains species not in sf or more copies of any
            species than sf contains.

        Raises
        ------
        ValueError
            If sf is not a ComplexSpecies.

        Notes
        -----
        This method compares the stoichiometry of species in s0 and sf to
        determine what needs to be added. If s0 contains more of any
        species than sf, or contains species not in sf, None is returned.

        """
        if not isinstance(sf, ComplexSpecies):
            raise ValueError(f"sf must be a ComplexSpecies. Recieved {sf}")

        species_to_add = []
        for s in sf.species_set:
            if s0.monomer_eq(
                s
            ):  # this is used instead of == to deal with the
                # potential for different parents
                s0_count = 1
            elif isinstance(s0, ComplexSpecies):
                s0_count = s0.monomer_count(
                    s
                )  # This is used instead of .count to deal with
                # species with different parents
            else:
                s0_count = 0

            # Add the correct number of s to the list
            sf_count = sf.monomer_count(s)
            if s0_count < sf_count:
                species_to_add += (sf_count - s0_count) * [s]
            elif s0_count > sf_count:
                species_to_add = None  # In this case, s0 contains
                # more stuff than sf, so
                # nothing should be returned
                break
            else:
                pass  # if they have the same number, do not add it

        # s0 contains more or different species than sf, return None
        if (not isinstance(s0, ComplexSpecies)) and (
            sf.monomer_count(s0) == 0
        ):
            species_to_add = None
        elif (
            isinstance(s0, ComplexSpecies)
            and any(
                [
                    s0.monomer_count(s) > sf.monomer_count(s)
                    for s in s0.species_set
                ]
            )
            and sf.monomer_count(s0) == 0
        ):
            species_to_add = None
        elif len(species_to_add) == 0:
            species_to_add = None

        return species_to_add

    def get_combinations_between(self, s0, sf):
        """Get all binding combinations to form complex sf from s0.

        Enumerates all possible binding orders (permutations) to construct
        the final complex sf starting from s0, generating tuples of
        (binder, bindee, complex_species) for each binding step.

        Parameters
        ----------
        s0 : Species or ComplexSpecies
            Starting species or complex.
        sf : ComplexSpecies
            Target final complex.

        Returns
        -------
        list of tuple
            List of (binder, bindee, complex_species) tuples representing
            all possible binding combinations. Each tuple represents one
            binding step. Returns empty list if no combinations are
            possible.

        Notes
        -----
        The method:

        1. Computes which species need to be added to s0 to form sf
        2. Generates all permutations of these species (different binding
           orders)
        3. For each permutation, creates binding steps: binder + bindee -->
           complex
        4. Filters out any combinations involving excluded_states

        Examples
        --------
        If s0 = A and sf = Complex([A, B, C]), and no exclusions:

        Returns combinations for different binding orders:

        - Order 1: (B, A, [A,B]), (C, [A,B], [A,B,C])
        - Order 2: (C, A, [A,C]), (B, [A,C], [A,B,C])

        """
        species_to_add = self.compute_species_to_add(s0, sf)

        if species_to_add is None or len(species_to_add) == 0:
            return []
        else:
            # combinations (binder, bindee, complex_species) to be returned
            combinations = []
            # get all permutations of the species for different binding orders
            perms = permutations(species_to_add, len(species_to_add))

            # Iterate through permutations
            for perm in perms:
                bindee = s0
                for s in perm:
                    binder = s
                    if isinstance(bindee, ComplexSpecies):
                        s_list = list(bindee.species)
                    else:
                        s_list = [bindee]

                    s_list.append(s)
                    cs = Complex(s_list)

                    # append the new combination if it isn't excluded
                    if (
                        binder not in self.excluded_states
                        and bindee not in self.excluded_states
                        and cs not in self.excluded_states
                    ):
                        combinations.append((binder, bindee, cs))
                    # update bindee
                    bindee = cs

            return combinations

    def update_species(self):
        """Use 'binding' mechanism to generate combinatorial species.

        Uses the 'binding' mechanism to generate species for all possible
        binding combinations between initial_states and final_states,
        optionally constrained by intermediate_states and excluding
        excluded_states.

        Returns
        -------
        list of Species
            List of all unique species generated, including initial states,
            final states, and all intermediate complexes along binding
            pathways.

        Notes
        -----
        The method handles two cases:

        With intermediate_states:

            1. Generate combinations: initial_states --> intermediate_states
            2. Generate combinations: intermediate_states --> final_states

        Without intermediate_states:

            Generate combinations: initial_states --> final_states directly

        Duplicate species are automatically removed from the final list.
        The combination_dict is populated during this process for use by
        `update_reactions`.

        """
        mech_b = self.get_mechanism('binding')
        species = []
        species_added_dict = {}  # save which combinations have
        # already been added
        self.combination_dict = {}  # this should be recomputed every
        # updated species

        # If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    # Get combinatorial species between s0 and si
                    if (s0, si) not in self.combination_dict:
                        self.combination_dict[s0, si] = (
                            self.get_combinations_between(s0, si)
                        )

                    # iterate through combinations of species between
                    # s0 and si
                    for binder, bindee, cs in self.combination_dict[s0, si]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

            for si in self.intermediate_states:
                for sf in self.final_states:
                    # Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = (
                            self.get_combinations_between(si, sf)
                        )

                    # iterate through combinations of species between
                    # si and sf
                    for binder, bindee, cs in self.combination_dict[si, sf]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

        # If there are no intermediate restrictions, compute
        # combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    # Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = (
                            self.get_combinations_between(s0, sf)
                        )

                    # iterate through combinations of species between
                    # s0 and sf
                    for binder, bindee, cs in self.combination_dict[s0, sf]:
                        if (binder, bindee, cs) not in species_added_dict:
                            species_added_dict[binder, bindee, cs] = True
                            species += mech_b.update_species(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

        return list(set(species))

    def update_reactions(self):
        """Use 'binding' mechanism to generate combinatorial reactions.

        Uses the 'binding' mechanism to generate reactions for all possible
        binding combinations between initial_states and final_states,
        optionally constrained by intermediate_states and excluding
        excluded_states.

        Returns
        -------
        list of Reaction
            List of all binding reactions (forward and reverse) along all
            enumerated pathways.

        Notes
        -----
        The method handles two cases:

        With intermediate_states:

            1. Generate reactions: initial_states <--> intermediate_states
            2. Generate reactions: intermediate_states <--> final_states

        Without intermediate_states:
            Generate reactions: initial_states <--> final_states directly

        Duplicate reactions are automatically filtered out. The method uses
        combination_dict computed by update_species() or computes it if
        needed. Reactions are symmetric, so (binder, bindee, complex) and
        (bindee, binder, complex) are treated as duplicates.

        """
        mech_b = self.get_mechanism('binding')
        reactions = []
        reactions_added_dict = {}  # save which combinations have
        # already been added in order to
        # not add duplicates

        # If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    # Get combinatorial species between s0 and si
                    if (s0, si) not in self.combination_dict:
                        self.combination_dict[s0, si] = (
                            self.get_combinations_between(s0, si)
                        )

                    # iterate through combinations of species between
                    # s0 and si
                    for binder, bindee, cs in self.combination_dict[s0, si]:
                        if (
                            binder,
                            bindee,
                            cs,
                        ) not in reactions_added_dict and (
                            bindee,
                            binder,
                            cs,
                        ) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

            for si in self.intermediate_states:
                for sf in self.final_states:
                    # Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = (
                            self.get_combinations_between(si, sf)
                        )

                    # iterate through combinations of species between
                    # si and sf
                    for binder, bindee, cs in self.combination_dict[si, sf]:
                        if (
                            binder,
                            bindee,
                            cs,
                        ) not in reactions_added_dict and (
                            bindee,
                            binder,
                            cs,
                        ) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

        # If there are no intermediate restrictions, compute
        # combinations in one step
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    # Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = (
                            self.get_combinations_between(s0, sf)
                        )

                    # iterate through combinations of species between
                    # s0 and sf
                    for binder, bindee, cs in self.combination_dict[s0, sf]:
                        if (
                            binder,
                            bindee,
                            cs,
                        ) not in reactions_added_dict and (
                            bindee,
                            binder,
                            cs,
                        ) not in reactions_added_dict:
                            reactions_added_dict[binder, bindee, cs] = True
                            reactions += mech_b.update_reactions(
                                binder=binder,
                                bindee=bindee,
                                complex_species=cs,
                                component=self,
                                part_id=self.name
                                + '_'
                                + str(binder)
                                + '_'
                                + str(bindee),
                            )

        return reactions
