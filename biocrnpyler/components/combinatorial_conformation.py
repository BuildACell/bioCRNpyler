#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import warnings
from itertools import permutations

from ..components.dna.promoter import Promoter
from ..core.component import Component
from ..core.species import Complex, ComplexSpecies, PolymerConformation
from ..mechanisms.conformation import One_Step_Reversible_Conformation_Change


class CombinatorialConformation(Component):
    """Polymer conformation with combinatorial internal binding complexes.

    A `CombinatorialConformation` component represents a polymer
    conformation (made of one unique OrderedPolymerSpecies) with multiple
    internal complexes that can bind and unbind in many different ways.
    Unlike `CombinatorialComplex` where individual species are added one at
    a time, this component adds groups of species in single steps to form
    the appropriate complexes. Uses a 'conformation_change' mechanism.

    Parameters
    ----------
    final_states : PolymerConformation or list of PolymerConformation
        One or more final polymer conformations to be formed. All must
        contain the same unique OrderedPolymerSpecies.
    initial_states : list of PolymerConformation, optional
        Initial polymer conformations that can bind/unbind to become
        final_states. If None or empty, defaults to the bare polymer
        without complexes.
    intermediate_states : list of PolymerConformation, optional
        Allowed intermediate conformations formed when converting
        initial_states to final_states. If None, all possible
        intermediate conformations are enumerated.
    excluded_states : list of PolymerConformation, optional
        Polymer conformations that will NOT be formed during enumeration.
        If None, no conformations are excluded.
    state_part_ids : dict, optional
        Dictionary mapping PolymerConformation to string, used to generate
        shorter part-ids for conformations.
    name : str, optional
        Name of the component. If None, uses the internal polymer name.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    final_states : list of PolymerConformation
        List of final conformation states.
    initial_states : list of PolymerConformation
        List of initial conformation states.
    intermediate_states : list of PolymerConformation or None
        List of allowed intermediate conformations, or None if
        unrestricted.
    excluded_states : list of PolymerConformation
        List of excluded conformations.
    internal_polymer : OrderedPolymerSpecies
        The unique polymer species common to all conformations.
    state_part_ids : dict
        Dictionary for custom part-id naming.
    combination_dict : dict
        Dictionary storing computed conformation changes.

    See Also
    --------
    CombinatorialComplex : Combinatorial binding of simple complexes.
    PolymerConformation : Species subclass for polymer conformations.
    Component : Base class for biomolecular components.

    Notes
    -----
    Key differences from `CombinatorialComplex`:

    - Operates on `PolymerConformation` objects instead of simple `Species`
    - All conformations must share the same `OrderedPolymerSpecies`
    - Adds groups of species simultaneously to form complexes
    - Uses 'conformation_change' mechanism instead of 'binding'

    Reaction generation: The component generates conformation change
    reactions based on constraints:

    - Without intermediate_states:
      initial_states <--> final_states

    - With intermediate_states:
      initial_states <--> intermediate_states <--> final_states

    Validation requirements: All conformations must:

    1. Be PolymerConformation objects
    2. Contain exactly one unique OrderedPolymerSpecies
    3. Have the same internal polymer

    Examples
    --------
    Create a simple conformational change system:

    >>> A, B, C, S = (bcp.Species(s) for s in ['A', 'B', 'C', 'S'])
    >>> pc = bcp.PolymerConformation(polymer=[A, A, B, C])
    >>> # Form a complex A:B by binding positions 0 and 2
    >>> c1 = bcp.Complex([pc.polymers[0][0], pc.polymers[0][2]])
    >>> pc1 = c1.parent
    >>> # Form two complexes: A:B and A:C:S (S is external)
    >>> c2 = bcp.Complex([pc1.polymers[0][1], pc1.polymers[0][3], S])
    >>> pc2 = c2.parent
    >>> # Create component to enumerate reactions
    >>> cc = bcp.CombinatorialConformation(
    ...     final_states=pc2,
    ...     parameters={'kf': 1, 'kr': 0.01})

    Using a Mixture to generate species and reactions:

    >>> mixture = bcp.Mixture(components=[cc])
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        final_states,
        initial_states=None,
        intermediate_states=None,
        excluded_states=None,
        state_part_ids=None,
        name=None,
        **kwargs,
    ):
        if state_part_ids is None:
            self.state_part_ids = {}
        else:
            self.state_part_ids = state_part_ids
        self.internal_polymer = None  # set inside final_states setter
        self.final_states = final_states
        self.initial_states = initial_states
        self.intermediate_states = intermediate_states
        self.excluded_states = excluded_states

        if name is None:
            name = str(self.internal_polymer.name)
        Component.__init__(
            self,
            name=name,
            default_mechanism=One_Step_Reversible_Conformation_Change(),
            **kwargs,
        )

    # Helper function to assert the correct class type
    def _assert_conformation(self, states, input_name='states'):
        """Validate that states are proper PolymerConformations.

        Parameters
        ----------
        states : list
            List of states to validate.
        input_name : str, default='states'
            Name of the parameter being validated (for error messages).

        Raises
        ------
        ValueError
            If states are not PolymerConformations, do not contain exactly
            one polymer, or do not share the same OrderedPolymerSpecies.

        Notes
        -----
        Sets self.internal_polymer on first call if not already set.

        """
        if not all([isinstance(s, PolymerConformation) for s in states]):
            raise ValueError(
                f"{input_name} must be a list of PolymerConformations. "
                f"Recieved: {states}."
            )
        if not all([len(s.polymers) == 1 for s in states]):
            raise ValueError(
                f"All PolymerConformations in {input_name} must contain "
                f"a single unique internal OrderedPolymerSpecies. "
                f"Recieved: {states}."
            )

        if self.internal_polymer is None:
            self.internal_polymer = copy.deepcopy(states[0].polymers[0])
            self.internal_polymer.parent = None

        if not all(
            [
                str(p) == str(self.internal_polymer)
                for s in states
                for p in s.polymers
            ]
        ):
            raise ValueError(
                f"All PolymerConformations in {input_name} must contain "
                f"a single unique internal OrderedPolymerSpecies. "
                f"Recieved: {states}."
            )

    def get_species(self):
        """Get the bare polymer conformation.

        Returns
        -------
        PolymerConformation
            The internal polymer without any complexes.

        """
        return PolymerConformation(polymer=self.internal_polymer)

    # Getters and setters
    # Final States stores the end complexes that will be formed
    @property
    def final_states(self):
        """List of final conformation states.

        Returns
        -------
        list of PolymerConformation

        """
        return self._final_states

    @final_states.setter
    def final_states(self, final_states):
        """Set the final conformation states.

        Parameters
        ----------
        final_states : PolymerConformation or list of PolymerConformation
            Final conformation(s) to be formed.

        Raises
        ------
        ValueError
            If validation fails (see _assert_conformation).

        """
        final_states = self.set_species(final_states)
        if not isinstance(final_states, list):
            final_states = [final_states]

        self._assert_conformation(final_states, 'final_states')
        self._final_states = final_states

    # Initial states stores the starting states used in binding reactions
    @property
    def initial_states(self):
        """List of initial conformation states.

        Returns
        -------
        list of PolymerConformation
        """
        return self._initial_states

    @initial_states.setter
    def initial_states(self, initial_states):
        """Set the initial conformation states.

        Parameters
        ----------
        initial_states : list of PolymerConformation, optional
            Initial conformations. If None or empty, defaults to bare
            polymer conformation.

        Raises
        ------
        ValueError
            If validation fails (see _assert_conformation).

        """
        # set initial states
        if initial_states is None or len(initial_states) == 0:
            self._initial_states = [
                PolymerConformation(polymer=self.internal_polymer)
            ]
        else:
            initial_states = self.set_species(initial_states)
            if not isinstance(initial_states, list):
                initial_states = [initial_states]

            # all initial_states must be PolymerConformation
            self._assert_conformation(initial_states, 'initial_states')

            self._initial_states = initial_states

    # Intermediate states allows the user to restrict the complexes formed
    # between the intial state and final state
    @property
    def intermediate_states(self):
        """List of allowed intermediates.

        Returns
        -------
        list of PolymerConformation or None

        """
        return self._intermediate_states

    @intermediate_states.setter
    def intermediate_states(self, intermediate_states):
        """Set the allowed intermediate conformation states.

        Parameters
        ----------
        intermediate_states : list of PolymerConformation, optional
            Allowed intermediate conformations. If None, all possible
            intermediates are enumerated.

        Raises
        ------
        ValueError
            If validation fails (see _assert_conformation).

        """
        if intermediate_states is None:
            self._intermediate_states = None
        else:
            intermediate_states = self.set_species(intermediate_states)
            if not isinstance(intermediate_states, list):
                intermediate_states = [intermediate_states]

            # All intermediate_states must be PolymerConformations
            self._assert_conformation(
                intermediate_states, 'intermediate_states'
            )

            self._intermediate_states = intermediate_states

    # excluded_states are PolymerConformations which are not allowed to form
    @property
    def excluded_states(self):
        """List of excluded conformations.

        Returns
        -------
        list of PolymerConformation

        """
        return self._excluded_states

    @excluded_states.setter
    def excluded_states(self, excluded_states):
        """Set the excluded conformation states.

        Parameters
        ----------
        excluded_states : list of PolymerConformation, optional
            Conformations that are NOT allowed to form. If None, no
            exclusions (empty list).

        Raises
        ------
        ValueError
            If validation fails (see _assert_conformation).

        """
        if excluded_states is None:
            self._excluded_states = []
        else:
            # All excluded states must be PolymerConformations
            excluded_states = self.set_species(excluded_states)
            if not isinstance(excluded_states, list):
                excluded_states = [excluded_states]
            self._assert_conformation(excluded_states, 'excluded_states')
            self._excluded_states = excluded_states

    def compute_species_changes(self, s0, sf):
        """Compute changes needed to convert conformation s0 into sf.

        Analyzes what species need to be added and which complexes need to
        be merged to transform the initial conformation s0 into the final
        conformation sf. Assumes both conformations share the same
        underlying polymer.

        Parameters
        ----------
        s0 : PolymerConformation
            Starting conformation.
        sf : PolymerConformation
            Target final conformation.

        Returns
        -------
        tuple of (dict, dict) or False
            Returns False if s0 cannot be additively transformed into sf.
            Otherwise returns (species_changes, merged_complexes) where:

            - species_changes: dict mapping (complex, positions) to list of
              external species to add
            - merged_complexes: dict mapping (complex, positions) to list of
              complexes from s0 that merge to form sf

        Notes
        -----
        Returns False if:

        - s0 has more complexes at any position than sf
        - Any complex in sf cannot be formed additively from s0

        """
        # print("computing species changes between", s0, "and", sf)
        if any(
            [
                len(s0.get_complexes_at(0, i))
                > len(sf.get_complexes_at(0, i))
                for i in range(len(self.internal_polymer))
            ]
        ):
            species_changes = False
        else:
            species_changes = {}
            merged_complexes = {}

            # Figure out what happens to create each complex in sf
            for cf in sf.complexes:
                pf_inds = sf.get_polymer_positions(cf, 0)

                # initially, assume all external Species should be
                # added species are stored as (s, index) in order to
                # differentiate species in different parts of the
                # polymer
                cf_species_and_inds = [
                    (cf.species[i], pf_inds[i])
                    for i in range(len(cf.species))
                ]
                # string version used to check for equality
                cf_species_and_inds_str = [
                    str(i) for i in cf_species_and_inds
                ]
                sub_complex_species_and_inds = []
                sub_complex_species_and_inds_str = []
                merged_complexes[cf, pf_inds] = []

                # Iterate through speces in the polymer sf
                for ind in pf_inds:
                    # get complex of the polymer at this location
                    found_complexes = s0.get_complexes_at(0, ind)
                    for c0 in found_complexes:
                        p0_inds = s0.get_polymer_positions(c0, 0)
                        # Ensure each subcomplex c0 is only evaluated
                        # once then add all the Species in the complex
                        # to the sub_complex_species_and_inds list
                        if (
                            (cf, pf_inds) not in merged_complexes
                            or c0 not in merged_complexes[cf, pf_inds]
                        ):
                            merged_complexes[cf, pf_inds].append(c0)

                            sub_complex_species_and_inds += [
                                (c0.species[i], p0_inds[i])
                                for i in range(len(c0.species))
                            ]
                            sub_complex_species_and_inds_str += [
                                str(i) for i in sub_complex_species_and_inds
                            ]

                # cf cannot be created additively from the complexes in s0
                sub_complex_exclusive_species = [
                    i
                    for i in sub_complex_species_and_inds_str
                    if i not in cf_species_and_inds_str
                ]
                cf_complex_exclusive_species = [
                    i
                    for i in cf_species_and_inds_str
                    if i not in sub_complex_species_and_inds_str
                ]
                if len(sub_complex_exclusive_species) > 0:
                    species_changes[cf, pf_inds] = False
                # keep track of all the external species added to
                # produce the new complex
                elif len(cf_complex_exclusive_species) > 0:
                    species_changes[cf, pf_inds] = [
                        cf_species_and_inds[i][0]
                        for i in range(len(cf_species_and_inds_str))
                        if (
                            cf_species_and_inds_str[i]
                            not in sub_complex_species_and_inds_str
                            and cf_species_and_inds[i][1] is None
                        )
                    ]
        # If any complex cannot be created, return False
        if (
            species_changes is False
            or any([species_changes[k] is False for k in species_changes])
            or (len(species_changes) == 0 and len(merged_complexes) == 0)
        ):
            # print("cannot convert.")
            return False
        else:
            # print("species_changes =", species_changes)
            # print("merged_complexes =", merged_complexes)
            return species_changes, merged_complexes

    def get_combinations_between(self, s0, sf):
        """Get all conformation change combinations from s0 to sf.

        Enumerates all possible orders of complex formation to transform
        conformation s0 into sf, generating tuples representing each step.

        Parameters
        ----------
        s0 : PolymerConformation
            Starting conformation.
        sf : PolymerConformation
            Target final conformation.

        Returns
        -------
        list of tuple
            List of (old_state, species_to_add, new_state) tuples
            representing all possible transformation pathways. Each tuple
            represents one conformation change step. Returns empty list if
            no valid pathways exist.

        Notes
        -----
        The method:

        1. Computes which species/complexes change between s0 and sf
        2. Generates all permutations (different formation orders)
        3. For each permutation, creates conformational change steps
        4. Filters out any combinations involving excluded_states

        Unlike `CombinatorialComplex`, this method adds groups of species
        simultaneously to form complete complexes at polymer positions.

        """
        # print("geting combinations between", s0, "and", sf)
        X = self.compute_species_changes(s0, sf)

        if X is False:
            return []
        else:
            species_changes, merged_complexes = X
            changes_list = []
            for cf in sf.complexes:
                inds = sf.get_polymer_positions(cf, 0)
                # print(
                #     "     checking ", cf, "inds =", inds,
                #     "in species_changes", (cf, inds) in species_changes,
                #     "in merged_complexes", (cf, inds) in merged_complexes)
                if (
                    (cf, inds) in merged_complexes
                    and len(merged_complexes[(cf, inds)]) == 1
                    and merged_complexes[(cf, inds)][0].monomer_eq(cf)
                ) or (cf, inds) not in merged_complexes:
                    complexes_to_merge = None
                else:
                    complexes_to_merge = merged_complexes.get((cf, inds))

                if (cf, inds) not in species_changes or len(
                    species_changes[cf, inds]
                ) == 0:
                    species_to_add = None
                else:
                    species_to_add = species_changes.get((cf, inds))

                # print(
                #     "(cf, inds, species_to_add, complexes_to_merge)=",
                #     (cf, inds, species_to_add, complexes_to_merge))
                if (
                    species_to_add is not None
                    or complexes_to_merge is not None
                ):
                    changes_list.append(
                        (cf, inds, species_to_add, complexes_to_merge)
                    )

            # get all permutations of the species for different binding orders
            perms = permutations(changes_list, len(changes_list))

            # combinations (binder, bindee, complex_species) to be returned
            combinations = []

            # Iterate through permutations
            for perm in perms:
                old_state = s0
                for cf, inds, species_to_add, complexes_to_merge in perm:
                    # Determine all the species that are bound together
                    species_list = []
                    for i, p_loc in enumerate(inds):
                        if p_loc is None:
                            species_list.append(cf.species[i])
                        else:
                            species_list.append(old_state.polymers[0][p_loc])

                    new_complex = Complex(species_list)

                    # Determine which complexes to remove from conformation
                    if complexes_to_merge is None:
                        merged_complexes = None
                    else:
                        merged_complexes = [
                            new_complex.parent.get_complex(c)
                            for c in complexes_to_merge
                        ]

                    # create a new polymer conformation by removing
                    # the merged complexes
                    if merged_complexes not in [[], None]:
                        new_state = new_complex.parent.copy_remove_complexes(
                            merged_complexes
                        )
                    else:
                        new_state = new_complex.parent

                    # append the new combination if it isn't excluded
                    if (
                        old_state not in self.excluded_states
                        and new_state not in self.excluded_states
                    ):
                        combinations.append(
                            (old_state, species_to_add, new_state)
                        )

                    # update bindee
                    old_state = new_state
            # print("combinations = ", combinations)
            return combinations

    def _get_part_id(self, state):
        """Get part ID for a conformation state.

        Parameters
        ----------
        state : PolymerConformation
            The conformation state.

        Returns
        -------
        str
            Custom part ID if state is in state_part_ids, otherwise string
            representation of the state.

        """
        if state in self.state_part_ids:
            return self.state_part_ids[state]
        else:
            return str(state)

    def update_species(self):
        """Use 'conformation_change' mechanism to generate species.

        Uses the 'conformation_change' mechanism to generate species for
        all possible conformation transformations between `initial_states`
        and `final_states`, optionally constrained by `intermediate_states`
        and excluding `excluded_states`.

        Returns
        -------
        list of Species
            List of all unique species generated, including
            polymer conformations and any additional species involved in
            conformation changes.

        Notes
        -----
        Duplicate species are automatically removed from the final list.
        The `combination_dict` is populated during this process for use by
        `update_reactions`.

        """
        mech_c = self.get_mechanism('conformation_change')
        species = []
        self.combination_dict = {}  # should recompute every updated species

        # If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    if si != s0:
                        # Get combinatorial species between s0 and si
                        if (s0, si) not in self.combination_dict:
                            self.combination_dict[s0, si] = (
                                self.get_combinations_between(s0, si)
                            )

                        # iterate through combinations of species
                        # between s0 and si
                        for (
                            s00,
                            additional_species,
                            sff,
                        ) in self.combination_dict[s0, si]:
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            species += mech_c.update_species(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

            for si in self.intermediate_states:
                for sf in self.final_states:
                    if si != sf:
                        # Get combinatorial species between si and sf
                        if (si, sf) not in self.combination_dict:
                            self.combination_dict[si, sf] = (
                                self.get_combinations_between(si, sf)
                            )

                        # iterate through combinations of species
                        # between si and sf
                        for (
                            s00,
                            additional_species,
                            sff,
                        ) in self.combination_dict[si, sf]:
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            species += mech_c.update_species(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

        # If there are no intermediate restrictions, compute
        # combinations in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    if s0 != sf:
                        # Get combinatorial species between s0 and sf
                        if (s0, sf) not in self.combination_dict:
                            self.combination_dict[s0, sf] = (
                                self.get_combinations_between(s0, sf)
                            )

                        # iterate through combinations of species
                        # between s0 and sf
                        for (
                            s00,
                            additional_species,
                            sff,
                        ) in self.combination_dict[s0, sf]:
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            species += mech_c.update_species(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

        return list(set(species))

    def update_reactions(self):
        """Use 'conformation_change' mechanism to generate reactions.

        Uses the 'conformation_change' mechanism to generate reactions for
        all possible conformation transformations between initial_states
        and final_states, optionally constrained by intermediate_states
        and excluding excluded_states.

        Returns
        -------
        list of Reaction
            List of all conformation change reactions (forward and reverse)
            along all enumerated pathways.

        Notes
        -----
        The method handles two cases:

        With intermediate_states:

            1. Generate reactions: initial_states <--> intermediate_states
            2. Generate reactions: intermediate_states <--> final_states

        Without intermediate_states: generate reactions:
        initial_states <--> final_states directly.

        Duplicate reactions are automatically filtered out using
        `reactions_added_dict`. The method uses `combination_dict` computed by
        `update_species` or computes it if needed.

        """
        mech_c = self.get_mechanism('conformation_change')
        reactions = []
        # save which combinations have already been added in order to
        # not add duplicates
        reactions_added_dict = {}

        # If there are intermediates, compute combinations in two steps
        if self.intermediate_states is not None:
            for s0 in self.initial_states:
                for si in self.intermediate_states:
                    # Get combinatorial species between s0 and si
                    if (s0, si) not in self.combination_dict:
                        self.combination_dict[s0, si] = (
                            self.get_combinations_between(s0, si)
                        )

                    # Iterate thru combinations of species between s0 and si
                    for s00, additional_species, sff in self.combination_dict[
                        s0, si
                    ]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            reactions += mech_c.update_reactions(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

            for si in self.intermediate_states:
                for sf in self.final_states:
                    # Get combinatorial species between si and sf
                    if (si, sf) not in self.combination_dict:
                        self.combination_dict[si, sf] = (
                            self.get_combinations_between(si, sf)
                        )

                    # Iterate thru combinations of species between si and sf
                    for s00, additional_species, sff in self.combination_dict[
                        si, sf
                    ]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            reactions += mech_c.update_reactions(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

        # If there are no intermediate restrictions, compute combinations
        # in onestep
        else:
            for s0 in self.initial_states:
                for sf in self.final_states:
                    # Get combinatorial species between s0 and sf
                    if (s0, sf) not in self.combination_dict:
                        self.combination_dict[s0, sf] = (
                            self.get_combinations_between(s0, sf)
                        )

                    # Iterate thru combinations of species between s0 and sf
                    for s00, additional_species, sff in self.combination_dict[
                        s0, sf
                    ]:
                        if (s00, sff) not in reactions_added_dict:
                            reactions_added_dict[s00, sff] = True
                            part_id = (
                                f"{self._get_part_id(s00)}"
                                + f"-{self._get_part_id(sff)}"
                            )
                            reactions += mech_c.update_reactions(
                                s0=s00,
                                sf=sff,
                                additional_species=additional_species,
                                component=self,
                                part_id=part_id,
                            )

        return reactions


class CombinatorialConformationPromoter(CombinatorialConformation, Promoter):
    """Combinatorial conformation with transcriptionally active states.

    A `CombinatorialConformationPromoter` combines `CombinatorialConformation`
    and `Promoter` functionality, creating a polymer with combinatorial
    conformations where certain conformations can transcribe/express
    RNA/protein products. Specific conformations can be designated as
    transcriptionally active ('on') or inactive ('off').

    Parameters
    ----------
    promoter_states : list of PolymerConformation
        Polymer conformations used by the promoter. These states are
        designated as either 'on' or 'off' based on promoter_states_on.
    promoter_location : int
        Index of the monomer in the polymer conformation that represents
        the promoter location on the polymer.
    promoter_states_on : bool, default=True
        If True, all `promoter_states` are transcribable. If False, all
        states except `promoter_states` are transcribable.
    activating_complexes : list of ComplexSpecies, optional
        Complexes that activate polymer conformations for transcription
        regardless of promoter_states.
    inactivating_complexes : list of ComplexSpecies, optional
        Complexes that inactivate polymer conformations, preventing
        transcription even if otherwise active.
    intermediate_states : list of PolymerConformation, optional
        Allowed intermediate conformations (see `CombinatorialConformation`).
    final_states : list of PolymerConformation, optional
        Final conformations (see `CombinatorialConformation`).
    name : str, default='CombinatorialConformationPromoter'
        Name of the component.
    **kwargs
        Additional keyword arguments passed to the parent class constructors.

    Attributes
    ----------
    promoter_states : list of PolymerConformation
        List of designated promoter states.
    promoter_states_on : bool
        Whether promoter_states are active or inactive.
    promoter_location : int
        Polymer position of the promoter.
    activating_complexes : list of ComplexSpecies
        Complexes that activate transcription.
    inactivating_complexes : list of ComplexSpecies
        Complexes that prevent transcription.
    conformation_species : list
        All conformation species (populated by update_species).

    See Also
    --------
    CombinatorialConformation : Base class for conformational changes.
    Promoter : Base class for transcription initiation.
    PolymerConformation : Polymer with internal complexes.

    Notes
    -----
    A conformation is transcriptionally active if:

    1. (conformation in promoter_states AND promoter_states_on=True) OR
       (conformation not in promoter_states AND promoter_states_on=False)
    2. OR any activating_complex is present in the conformation
    3. AND no inactivating_complex is present

    If inactivating_complex conflicts with active_state or
    active_complex, a warning is issued and transcription is prevented.

    The promoter location determines which DNA species from the polymer is
    used for transcription initiation.

    Examples
    --------
    Create a promoter (operon) with conformational regulation:

    >>> A, B, F = (bcp.Species(s) for s in ['A', 'B', 'F'])
    >>> op = bcp.PolymerConformation(polymer=[B, A, B])
    >>> OF0 = bcp.Complex([op.polymers[0][0], F]).parent  # F bound at pos'n 0
    >>> OF1 = bcp.Complex([op.polymers[0][2], F]).parent  # F bound at pos'n 2
    >>> OF2 = tbp.Complex([OF1.polymers[0][2], F]).parent # F bound to both
    >>> # Looped conformations
    >>> L0 = Complex([op.polymers[0][0], op.polymers[0][1], F]).parent
    >>> L1 = Complex([op.polymers[0][2], op.polymers[0][1], F]).parent
    >>> # Define fully bound looped states
    >>> L0F1 = bcp.Complex(
    ...     [OF1.polymers[0][0], OF1.polymers[0][1], F]).parent
    >>> L1F0 = bcp. Complex(
    ...     [OF0.polymers[0][2], OF0.polymers[0][1], F]).parent
    >>> # Create promoter with specific active states
    >>> ccp = bcp.CombinatorialConformationPromoter(
    ...     name="CCP",
    ...     intermediate_states=[OF0, OF1],
    ...     final_states=[OF2, L0F1, L1F0],
    ...     promoter_states=[L0F1, L1F0, L0, L1],  # transcribed states
    ...     promoter_states_on=True,
    ...     promoter_location=1
    ... )

    With repression (toggle `promoter_states_on`):

    >>> # Same setup as above, but with promoter_states_on=False
    >>> # Now only states NOT in promoter_states will transcribe
    >>> ccp = bcp.CombinatorialConformationPromoter(
    ...     name="CCP",
    ...     intermediate_states=[OF0, OF1],
    ...     final_states=[OF2, L0F1, L1F0],
    ...     promoter_states=[L0F1, L1F0, L0, L1],
    ...     promoter_states_on=False,
    ...     promoter_location=1
    ... )
    >>> # Use in a DNAassembly for transcription
    >>> assy = bcp.DNAassembly(
    ...     name="X", dna=op, promoter=ccp, rbs="rbs", protein="X")
    >>> mixture = bcp.Mixture(
    ...     components=[assy],
    ...     mechanisms=[bcp.SimpleTranscription(), bcp.SimpleTranslation()],
    ...     parameters={'kf': 1, 'kr': 0.01, 'ktx': 1, 'ktl': 1}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self,
        promoter_states,
        promoter_location,
        promoter_states_on=True,
        activating_complexes=None,
        inactivating_complexes=None,
        intermediate_states=None,
        final_states=None,
        name='CombinatorialConformationPromoter',
        **kwargs,
    ):
        Promoter.__init__(self, name, **kwargs)
        CombinatorialConformation.__init__(
            self,
            name=name,
            intermediate_states=intermediate_states,
            final_states=final_states,
            **kwargs,
        )

        self.promoter_states = promoter_states
        self.promoter_states_on = promoter_states_on

        if activating_complexes is None:
            self.activating_complexes = []
        else:
            self.activating_complexes = activating_complexes

        if inactivating_complexes is None:
            self.inactivating_complexes = []
        else:
            self.inactivating_complexes = inactivating_complexes

        assert all(
            [isinstance(c, ComplexSpecies) for c in self.activating_complexes]
        )
        assert all(
            [
                isinstance(c, ComplexSpecies)
                for c in self.inactivating_complexes
            ]
        )

        if promoter_location not in range(len(self.internal_polymer)):
            raise ValueError(
                f"promoter_location must be an index of the polymer "
                f"{self.internal_polymer}. Recieved {promoter_location}."
            )
        else:
            self.promoter_location = promoter_location

    # promoter_states are PolymerConformations which are used by the
    # Promoter class.  These can be ON or OFF
    @property
    def promoter_states(self):
        """List of designated promoter states.

        Returns
        -------
        list of PolymerConformation

        """
        return self._promoter_states

    @promoter_states.setter
    def promoter_states(self, promoter_states):
        """Set the promoter conformational states.

        Parameters
        ----------
        promoter_states : list of PolymerConformation, optional
            Conformations designated as promoter states. If None, empty
            list.

        Raises
        ------
        ValueError
            If validation fails (see _assert_conformation).

        """
        if promoter_states is None:
            self._promoter_states = []
        else:
            # All excluded states must be PolymerConformations
            promoter_states = self.set_species(promoter_states)
            if not isinstance(promoter_states, list):
                promoter_states = [promoter_states]
            self._assert_conformation(promoter_states, 'promoter_states')
            self._promoter_states = promoter_states

    def update_species(self):
        """Generate species for conformation changes and transcription.

        Generates species from both conformational changes (via
        CombinatorialConformation) and transcription (via Promoter) for
        conformations that are transcriptionally active.

        Returns
        -------
        list of Species
            List of all unique species including conformation states and
            transcription-related species (RNAP complexes, transcripts,
            etc.) for active conformations.

        Notes
        -----
        For each conformation, determines if it is transcriptionally active
        based on:

        - Whether it is in promoter_states (and promoter_states_on setting)
        - Presence of activating_complexes
        - Absence of inactivating_complexes

        Only active conformations generate transcription species via
        Promoter.update_species().

        """
        self.conformation_species = CombinatorialConformation.update_species(
            self
        )
        promoter_species = []
        for s in self.conformation_species:
            if isinstance(s, PolymerConformation):
                active_state = False
                if (
                    s in self.promoter_states and self.promoter_states_on
                ) or (
                    s not in self.promoter_states
                    and not self.promoter_states_on
                ):
                    active_state = True

                active_complex = False
                if any(
                    [
                        s.get_complex(c) is not None
                        for c in self.activating_complexes
                    ]
                ):
                    active_complex = True

                innactive_complex = False
                if any(
                    [
                        s.get_complex(c) is not None
                        for c in self.inactivating_complexes
                    ]
                ):
                    innactive_complex = True
                    if (
                        active_state and len(self.promoter_states) > 0
                    ) or active_complex:
                        warnings.warn(
                            "Inactive_complex conflicts with active_complex "
                            "or active_state in "
                            "CombinatorialConformationPromoter. "
                            "Defaulting to inactive."
                        )

                if (active_state or active_complex) and not innactive_complex:
                    self.dna_to_bind = s.polymers[0][self.promoter_location]
                    # Reset name for unique addressable promoter states
                    self.name = self._get_part_id(s)
                    promoter_species += Promoter.update_species(self)

        return list(set(promoter_species + self.conformation_species))

    def update_reactions(self):
        """Generate reactions for conformation changes and transcription.

        Generates reactions from both conformational changes (via
        CombinatorialConformation) and transcription (via Promoter) for
        conformations that are transcriptionally active.

        Returns
        -------
        list of Reaction
            List of all reactions including conformation change reactions
            and transcription reactions (RNAP binding, elongation, etc.)
            for active conformations.

        Notes
        -----
        For each conformation, determines if it is transcriptionally active
        using the same logic as update_species(). Only active
        conformations generate transcription reactions via
        `Promoter.update_reactions`.

        The component name is temporarily changed to a state-specific name
        for each conformation to ensure unique reaction identifiers.

        """
        if not hasattr(self, 'conformation_species'):
            self.update_species

        conformation_reactions = CombinatorialConformation.update_reactions(
            self
        )
        promoter_reactions = []
        old_name = self.name
        for s in self.conformation_species:
            if isinstance(s, PolymerConformation):
                active_state = False
                if (
                    s in self.promoter_states and self.promoter_states_on
                ) or (
                    s not in self.promoter_states
                    and not self.promoter_states_on
                ):
                    active_state = True

                active_complex = False
                if any(
                    [
                        s.get_complex(c) is not None
                        for c in self.activating_complexes
                    ]
                ):
                    active_complex = True

                innactive_complex = False
                if any(
                    [
                        s.get_complex(c) is not None
                        for c in self.inactivating_complexes
                    ]
                ):
                    innactive_complex = True
                    if (
                        active_state and len(self.promoter_states) > 0
                    ) or active_complex:
                        warnings.warn(
                            "inactive_complex conflicts with active_complex "
                            "or active_state in "
                            "CombinatorialConformationPromoter. "
                            "Defaulting to innactive."
                        )

                if (active_state or active_complex) and not innactive_complex:
                    self.dna_to_bind = s.polymers[0][self.promoter_location]
                    # Reset name for unique addressable promoter states
                    self.name = self._get_part_id(s)
                    promoter_reactions += Promoter.update_reactions(self)

        self.name = old_name
        return promoter_reactions + conformation_reactions
