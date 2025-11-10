################################################################
# DNA_construct: a higher level for construct compilation
# Author: Andrey Shur
# Latest update: 12/21/2020
#
################################################################


# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from ...core.component import Component
from ...core.polymer import OrderedMonomer, OrderedPolymer
from ...core.species import ComplexSpecies, OrderedPolymerSpecies, Species
from ...utils import all_comb
from ..basic import DNA, RNA


class Construct(Component, OrderedPolymer):
    """Base class for ordered genetic constructs with multiple parts.

    A Construct represents an ordered arrangement of genetic parts (promoters,
    RBS, coding sequences, terminators, etc.) that form a functional unit.
    This class provides the infrastructure for handling complex genetic
    constructs with multiple components, their enumeration, and generation
    of combinatorial variants. Constructs can be linear or circular and
    support both forward and reverse orientations of their constituent parts.

    Parameters
    ----------
    parts_list : list of list
        List of parts in the format [[part, direction], [part, direction],
        ...]  where each part must be an OrderedMonomer and direction is
        'forward' or 'reverse'.
    name : str, optional
        Name of the construct. If None, automatically generated from parts.
    circular : bool, default=False
        If True, the construct is circular (e.g., plasmid). If False, linear.
    mechanisms : dict or list, optional
        Custom mechanisms for this construct, overriding mixture defaults.
    parameters : dict, optional
        Parameter values specific to this construct.
    attributes : list of str, optional
        List of attribute tags for the construct.
    initial_concentration : float, optional
        Initial concentration of the construct species.
    component_enumerators : list, optional
        List of enumerator objects that generate construct variants.
    make_dirless_hash : bool, default=True
        If True, generates direction-independent hash for construct
        comparison.
    **kwargs
        Additional keyword arguments passed to Component constructor.

    Attributes
    ----------
    parts_list : list
        Ordered list of parts in the construct.
    circular : bool
        Whether the construct is circular.
    component_enumerators : list
        Enumerators for generating construct variants.
    out_components : list or None
        Cached list of output components from enumeration.
    predicted_rnas : list or None
        Cached list of predicted RNA species.
    predicted_proteins : list or None
        Cached list of predicted protein species.

    See Also
    --------
    DNA_construct : DNA-specific construct implementation.
    RNA_construct : RNA-specific construct implementation.
    DNA_part : Base class for individual DNA parts.
    OrderedPolymer : Base class for ordered polymer structures.

    Notes
    -----
    Constructs support several advanced features:

    - Part enumeration: Automatically generates all functional variants
      based on the parts present (e.g., all promoter-RBS combinations)
    - Combinatorial complexes: Generates all possible binding states
      when regulatory proteins bind to different parts
    - Direction-free comparison: Can identify equivalent constructs
      regardless of orientation or circular permutation

    The class maintains caches for enumerated components, RNA products, and
    protein products to avoid redundant computation.

    Examples
    --------
    Create a simple linear construct:

    >>> promoter = bcp.Promoter('ptet')
    >>> rbs = bcp.RBS('RBS_standard')
    >>> cds = bcp.CDS('GFP')
    >>> parts = [[promoter, 'forward'], [rbs, 'forward'], [cds, 'forward']]
    >>> construct = bcp.Construct(
    ...     parts_list=parts,
    ...     name='gene_circuit',
    ...     circular=False
    ... )

    Create a circular plasmid:

    >>> ori = bcp.DNA_part('p15A')
    >>> terminator = bcp.Terminator('BBa_B0022')
    >>> parts = [
    >>>     [ori, 'forward'], [promoter, 'forward'], [rbs, 'forward'],
    >>>     [cds, 'forward'], [terminator, 'forward']
    >>> ]
    >>> plasmid = bcp.Construct(
    ...     parts_list=parts,
    ...     name='pExpression',
    ...     circular=True
    ... )

    """

    def __init__(
        self,
        parts_list,
        name=None,
        circular=False,
        mechanisms=None,  # custom mechanisms
        parameters=None,  # customized parameters
        attributes=None,
        initial_concentration=None,
        component_enumerators=None,
        make_dirless_hash=True,
        **kwargs,
    ):
        if component_enumerators is None:
            component_enumerators = []
        self.component_enumerators = component_enumerators

        OrderedPolymer.__init__(self, parts_list, default_direction='forward')
        self.circular = circular
        if name is None:
            name = self.make_name()  # automatic naming
        # self.length = len(parts_list)  # RMM, 14 Sep '25: was not being set
        material_type = kwargs.pop('material_type', None)  # noqa: F841
        Component.__init__(
            self=self,
            name=name,
            mechanisms=mechanisms,
            parameters=parameters,
            attributes=attributes,
            initial_concentration=initial_concentration,
            **kwargs,
        )
        # TODO: material_type is not used; setting it here generates
        # an error self.material_type = material_type # RMM, 14 Sep
        # '25: was not being set
        self.update_parameters()
        self.transcripts = []
        self.out_components = None
        self.predicted_rnas = None
        self.predicted_proteins = None
        self.directionless_hash = None
        if make_dirless_hash:
            self.update_permutation_hash()

    @property
    def parts_list(self):
        return self.polymer

    def make_name(self):
        """Generate a systematic name for the construct based on its parts.

        Creates a name by concatenating all part names with underscores.
        Parts in reverse orientation are suffixed with '_r', and circular
        constructs are suffixed with '_o'.

        Returns
        -------
        str
            The generated construct name.

        Examples
        --------
        Linear construct with forward parts:

        >>> promoter = bcp.Promoter('pLac')
        >>> rbs = bcp.CDS('RBS1')
        >>> cds = bcp.CDS('GFP')
        >>> construct = bcp.DNA_construct(
        ...     [[promoter, 'forward'], [rbs, 'forward'], [cds, 'forward']]
        ... )
        >>> construct.make_name()
        'pLac_RBS1_GFP'

        Linear construct with reversed part:

        >>> construct = bcp.DNA_construct(
        ...     [[promoter, 'forward'], [rbs, 'reverse'], [cds, 'forward']]
        ... )
        >>> construct.make_name()
        'pLac_RBS1_r_GFP'

        Circular construct:

        >>> construct = bcp.DNA_construct(
        ...     [[promoter, 'forward'], [rbs, 'forward'], [cds, 'forward']],
        ...     circular=True
        ... )
        >>> construct.make_name()
        'pLac_RBS1_GFP_o'

        """
        output = ''
        outlst = []
        for part in self.parts_list:
            pname = part.name
            if part.direction == 'reverse':
                pname += '_r'
            outlst += [pname]
        output = '_'.join(outlst)
        if self.circular:
            output += '_o'
        return output

    def get_part(self, part=None, part_type=None, name=None, index=None):
        """Find and return a part from the construct by various criteria.

        Searches the construct's parts list for a part matching the given
        criteria. Only one search criterion should be provided at a time.

        Parameters
        ----------
        part : DNA_part, optional
            A specific DNA_part instance to find. Matches both type and name.
        part_type : type, optional
            Find all parts that are instances of this type (e.g., Promoter).
        name : str, optional
            Find part(s) with this exact name.
        index : int, optional
            Return the part at this position in parts_list.

        Returns
        -------
        DNA_part, list of DNA_part, or None
            Single matching part if exactly one match found. List of parts if
            multiple matches found. None if no matches found.

        Raises
        ------
        ValueError
            If multiple search criteria are provided simultaneously, or if
            invalid types are provided for parameters.

        Warns
        -----
        UserWarning
            If multiple matching parts are found (returns list).

        Notes
        -----
        The search is performed with the following priority:

        1. Index (direct lookup)
        2. Part instance (type and name must match)
        3. Name (string match)
        4. Type (isinstance check)

        Only one search criterion should be provided at a time.

        Examples
        --------
        Find a part by name:

        >>> promoter = bcp.Promoter('ptet')
        >>> rbs = bcp.RBS('RBS_standard')
        >>> cds = bcp.CDS('GFP')
        >>> parts = [
        ...     [promoter, 'forward'], [rbs, 'forward'], [cds, 'forward']
        ... ]
        >>> construct = bcp.Construct(
        ...     parts_list=parts,
        ...     name='gene_circuit',
        ...     circular=False
        ... )
        >>> promoter = construct.get_part(name='ptet')

        Get the third part in the construct:

        >>> third_part = construct.get_part(index=2)

        Find all RBS parts:

        >>> all_rbs = construct.get_part(part_type=bcp.RBS)

        """
        if [part, name, index, part_type].count(None) != 3:
            raise ValueError(
                f"get_component requires a single keyword. "
                f"Recieved component={part}, name={name}, index={index}."
            )
        if not (isinstance(part, DNA_part) or part is None):
            raise ValueError(
                f"component must be of type DNA_part. Recieved {part}."
            )
        if not (isinstance(part_type, type) or part_type is None):
            raise ValueError(
                f"part_type must be a type. Recieved {part_type}."
            )
        if not (isinstance(name, str) or name is None):
            raise ValueError(f"name must be of type str. Recieved {name}.")
        if not (isinstance(index, int) or index is None):
            raise ValueError(f"index must be of type int. Recieved {index}.")

        matches = []
        if index is not None:
            matches.append(self.parts_list[index])
        else:
            for comp in self.parts_list:
                if part is not None:
                    if type(part) is type(comp) and comp.name == part.name:
                        matches.append(comp)
                elif name is not None:
                    if comp.name == name:
                        matches.append(comp)
                elif part_type is not None:
                    if isinstance(comp, part_type):
                        matches.append(comp)
        if len(matches) == 0:
            return None
        elif len(matches) == 1:
            return matches[0]
        else:
            warn(
                "get_part found multiple matching components. "
                "A list has been returned."
            )
            return matches

    def reverse(self):
        """Reverse the construct without modifying the underlying DNA.

        Creates a reversed representation of the construct where all parts
        are in reverse order and each part's direction is flipped. This is
        useful for generating the reverse complement of a construct.

        Returns
        -------
        Construct
            Returns self after reversing the parts list and flipping
            directions.

        Notes
        -----
        This method modifies the construct in place by:

        1. Reversing the order of parts in the parts_list
        2. Flipping each part's direction (forward <--> reverse)

        The underlying DNA sequence is not modified, only the representation
        changes.

        Examples
        --------
        Reverse a simple construct:

        >>> promoter = bcp.Promoter('ptet')
        >>> gene = bcp.CDS('GFP')
        >>> construct = bcp.Construct(
        ...    [[promoter, 'forward'], [gene, 'forward']])
        >>> construct.reverse()
        Construct = GFP_r_pLac_r

        """
        # Reverses everything, without actually changing the DNA.
        # Also updates the name and stored, since this is now a different
        # Construct.
        OrderedPolymer.reverse(self)
        self.reset_stored_data()
        self.name = self.make_name()
        # self.update_base_species()
        return self

    def get_reversed(self):
        """Create a deep copy of the construct with reversed orientation.

        Returns a new construct that is the reverse complement of this
        construct, with all parts in reverse order and flipped directions.
        The original construct is not modified.

        Returns
        -------
        Construct
            A new construct object with reversed parts and directions.

        See Also
        --------
        reverse : Reverse the construct in place.

        Examples
        --------
        Get reversed version without modifying original:

        >>> promoter = bcp.Promoter('ptet')
        >>> cds = bcp.CDS('GFP')
        >>> original = bcp.Construct(
        ...     [[promoter, 'forward'], [cds, 'forward']])
        >>> original.get_reversed()
        Construct = GFP_r_ptet_r

        """
        outcon = copy.deepcopy(self)
        outcon.reverse()
        return outcon

    def get_circularly_permuted(self, new_first_position):
        """Create a circularly permuted version of this construct.

        Returns a new construct where the circular ordering of parts starts
        at a different position. Only valid for circular constructs.

        Parameters
        ----------
        new_first_position : int
            The index of the part that should become the first position in
            the new construct.

        Returns
        -------
        DNA_construct
            A new circular DNA_construct with parts reordered starting from
            the specified position.

        Raises
        ------
        ValueError
            If the construct is linear (circular permutation only applies
            to circular constructs).

        Notes
        -----
        The parts list is rotated so that `parts_list[new_first_position]`
        becomes the new first element, maintaining the circular structure.

        Examples
        --------
        Permute a circular plasmid:

        >>> ori = bcp.DNA_part('p15A')
        >>> promoter = bcp.Promoter('ptet')
        >>> cds = bcp.CDS('GFP')
        >>> plasmid = bcp.DNA_construct(
        ...     [[ori, 'forward'], [promoter, 'forward'], [cds, 'forward']],
        ...     circular=True
        ... )
        >>> plasmid
        DNA_construct = p15A_ptet_GFP_o
        >>> plasmid.get_circularly_permuted(1)
        DNA_construct = ptet_GFP_p15A_o

        """
        if not self.circular:
            return ValueError("cannot circularly permute linear construct")
        else:
            return DNA_construct(
                self.parts_list[new_first_position:]
                + self.parts_list[:new_first_position],
                circular=True,
                component_enumerators=self.component_enumerators,
                attributes=self.attributes,
                mechanisms=self.mechanisms,
                mixture=self.mixture,
            )

    def set_mixture(self, mixture):
        """Set the mixture containing this construct and all its parts.

        Propagates the mixture reference to all parts in the construct,
        ensuring they share access to the same parameter database and
        mechanisms.

        Parameters
        ----------
        mixture : Mixture
            The mixture object that contains this construct.

        Notes
        -----
        This method ensures that all parts in the construct have access to
        the same mixture-level parameters and mechanisms, maintaining
        consistency across the entire construct.

        """
        self.mixture = mixture
        for part in self.parts_list:
            part.set_mixture(mixture)

    def update_permutation_hash(self):
        """Update the direction-independent hash for this construct.

        Generates a unique string representation of the construct that is
        invariant to direction (forward/reverse) and circular permutations
        (for circular constructs). This enables comparison of functionally
        equivalent constructs regardless of their representation.

        Notes
        -----
        The hash is stored in the `directionless_hash` attribute and is used
        for identifying equivalent constructs that differ only in orientation
        or circular starting position.

        The hash is computed using the `omnihash` class method, which finds
        the most alphabetically ordered representation of the construct.

        See Also
        --------
        omnihash : Class method that computes the direction and rotation-free
                   hash.

        """
        self.directionless_hash, _, _ = Construct.omnihash(self)

    def update_base_species(self, base_name=None, attributes=None):
        """Update the base species representation of this construct.

        Sets the `base_species` attribute to a Species object representing
        the construct's primary chemical species in the CRN.

        Parameters
        ----------
        base_name : str, optional
            Name for the base species. If None, uses the construct's name.
        attributes : list of str, optional
            Attribute tags to add to the species.

        Notes
        -----
        The base species serves as the chemical representation of the
        construct in CRN compilation, with material_type matching the
        construct's material_type (e.g., 'dna' for DNA_constructs).

        """
        if base_name is None:
            self.base_species = self.set_species(
                self.name,
                material_type=self.material_type,
                attributes=attributes,
            )
        else:
            self.base_species = self.set_species(
                base_name,
                material_type=self.material_type,
                attributes=attributes,
            )

    def update_parameters(self, overwrite_parameters=True):
        """Update parameters for the construct and all its parts.

        Propagates parameter updates from the construct's parameter database
        to all parts in the parts list, ensuring consistent parameters
        throughout the construct.

        Parameters
        ----------
        overwrite_parameters : bool, default=True
            If True, new parameter values overwrite existing ones in the
            parts. If False, existing parameters in parts are preserved.

        Notes
        -----
        This method:

        1. Updates the construct's own parameters via the parent Component
           class
        2. Propagates these parameters to each part in the construct

        This ensures that all parts have access to the same parameter values
        unless explicitly overridden at the part level.

        """
        Component.update_parameters(
            self=self, parameter_database=self.parameter_database
        )
        for part in self.parts_list:
            part.update_parameters(
                parameter_database=self.parameter_database,
                overwrite_parameters=overwrite_parameters,
            )

    def add_mechanism(
        self,
        mechanism,
        mech_type=None,
        overwrite=False,
        optional_mechanism=False,
    ):
        """Add a mechanism to the construct and all its parts.

        Adds the mechanism to the construct's mechanism dictionary and
        propagates it to all parts in the construct.

        Parameters
        ----------
        mechanism : Mechanism
            The mechanism object to add.
        mech_type : str, optional
            The mechanism type key. If None, uses the mechanism's
            `mechanism_type` attribute.
        overwrite : bool, default=False
            If True, overwrites existing mechanisms with the same type. If
            False, raises ValueError for duplicate types.
        optional_mechanism : bool, default=False
            If True, suppresses ValueError when a mechanism key conflict
            occurs and `overwrite` is False.

        Notes
        -----
        This method:

        1. Adds the mechanism to the construct via the parent Component
           class
        2. Propagates the mechanism to each part in the construct

        This ensures mechanism consistency across the entire construct.

        """
        Component.add_mechanism(
            self,
            mechanism,
            mech_type=mech_type,
            overwrite=overwrite,
            optional_mechanism=optional_mechanism,
        )
        for part in self.parts_list:
            part.add_mechanism(
                mechanism,
                mech_type=mech_type,
                overwrite=overwrite,
                optional_mechanism=optional_mechanism,
            )

    def __repr__(self):
        return 'Construct = ' + self.make_name()

    def __contains__(self, obj2):
        """Check if the construct contains a specific part.

        Tests whether a DNA_part or copy of a DNA_part exists in this
        construct's parts list.

        Parameters
        ----------
        obj2 : DNA_part
            The part to search for in the construct.

        Returns
        -------
        bool
            True if the part (or a copy with the same name and type) is in
            the construct, False otherwise.

        Notes
        -----
        This method supports two types of containment checks:

        1. Direct membership: The exact part object is in the construct
           (checked via `obj2.parent == self`)
        2. Copy membership: A different part object with the same type and
           name exists in the construct

        Examples
        --------
        Check if construct contains a part:

        >>> promoter = bcp.Promoter('ptet')
        >>> rbs = bcp.RBS('RBS_standard')
        >>> cds = bcp.CDS('GFP')
        >>> parts = [
        ...    [promoter, 'forward'], [rbs, 'forward'], [cds, 'forward']
        ... ]
        >>> construct = bcp.Construct(
        ...     parts_list=parts,
        ...     name='gene_circuit'
        ... )
        >>> promoter in construct
        True

        >>> unknown_part = bcp.Promoter('plac')
        >>> unknown_part in construct
        False

        """
        if isinstance(obj2, DNA_part):
            # if we got a DNA part it could mean one of two things:
            # 1 we want to know if a dna part is anywhere
            # 2 we want to know if a specific DNA part is in here
            # this is complicated by the fact that we want to have the
            # same DNA part be reusable in many locations
            if obj2.parent == self:
                # the object should already know if it's a part of me
                return True
            elif obj2.parent is None:
                # this object has been orphaned.  that means we are
                # looking for matching objects in any position
                new_obj2 = copy.copy(obj2).unclone()
                uncloned_list = [
                    copy.copy(a).unclone() for a in self.parts_list
                ]
                return new_obj2 in uncloned_list
            else:
                return False
        elif isinstance(obj2, str):
            # if we get a string, that means we want to know if the
            # name exists anywhere
            return obj2 in str(self)

    def get_species(self):
        """Returns species of DNA construct, using OrderedPolymerSpecies."""
        ocomplx = []
        for part in self.parts_list:
            partspec = copy.copy(part.dna_species)
            ocomplx += [partspec.set_dir(part.direction)]
        out_species = OrderedPolymerSpecies(
            ocomplx, circular=self.circular, material_type=self.material_type
        )

        return out_species

    def located_allcomb(self, spec_list):
        """Generate all combinatorial placement dictionaries for species.

        Creates all possible combinations of species placements when multiple
        species can bind at different positions in the construct.

        Parameters
        ----------
        spec_list : list of Species
            List of species that have position attributes, typically
            ComplexSpecies that can bind at specific locations.

        Returns
        -------
        list of dict
            List of dictionaries where each dict maps positions (int) to
            [species, direction] pairs representing one possible combinatorial
            binding state.

        Notes
        -----
        This method handles the combinatorics of placing multiple binding
        species at different positions. For example:

        - Species A binds at position 0
        - Species B binds at position 3
        - Species C also binds at position 0

        The method generates all valid combinations: `{0: [A, direction]}`,
        `{0: [C, direction]}`, `{0: [A, direction], 3: [B, direction]}`,
        `{0: [C, direction], 3: [B, direction]}`, etc.

        The algorithm:

        1. Extracts positions from spec_list
        2. Groups species by position into prototype_list
        3. Generates all position combinations via all_comb
        4. For each position combination, generates all possible species
           selections at those positions

        """
        # first we have to construct the list we are tracing paths through

        spec_indexes = [a.position for a in spec_list]  # extract all
        # indexes the following takes apart the lists because i don't
        # yet know how to deal with multiple binders at the same time
        compacted_indexes = sorted(list(set(spec_indexes)))

        prototype_list = [None] * len(compacted_indexes)
        for spec in spec_list:
            # now, spec is a list which contains all the elements
            # which are defined for each variant.
            #
            # go through every element and put it in the right place
            proto_ind = compacted_indexes.index(
                spec.position
            )  # where to put it?
            if prototype_list[proto_ind] is None:
                # if nothing's been placed here, then create a list
                prototype_list[proto_ind] = [spec]
            else:
                # if something is already here, then add to the list
                prototype_list[proto_ind] += [spec]
        # at this point we have a list that looks like this:
        # [[[part1,0],[part2,0]],[[part2,3]],[[part3,12],[part5,12]]
        # next step is to pick one of the first list (either [part1,0]
        # or [part2,0])
        # one of the second list (only [part2,3] is our option), one
        # of the third list, etc for all possible choices made this way
        comb_list = all_comb(compacted_indexes)

        def recursive_path(in_list):
            """Takes every possible "path" through a list of lists.

            For example:
            input:  [[A,B],[C],[D,E]]

            we have two options for the first position, one option for the
            second, and two options for the third.

            so, output:[[A,C,D],[B,C,D],[A,C,E],[B,C,E]]

            """
            if len(in_list) == 1:
                out_list = []
                for a in in_list[0]:
                    out_list += [[a]]
                return out_list
            elif len(in_list) == 0:
                return []
            else:
                out_list = []
                for a in in_list[0]:
                    out_list += [[a] + z for z in recursive_path(in_list[1:])]
                return out_list

        outlist = []
        for combo in comb_list:
            combo_sublists = []
            for combo_index in combo:
                combo_sublists += [
                    prototype_list[compacted_indexes.index(combo_index)]
                ]
            outlist += recursive_path(combo_sublists)
        outdict_list = []
        for combo in outlist:
            replacedict = {}
            for spec in combo:
                replacedict.update({spec.position: [spec, spec.direction]})
            outdict_list += [replacedict]
        return outdict_list

    def make_polymers(self, species_lists, backbone):
        """Create polymer species from combinatorial binding combinations.

        Generates OrderedPolymerSpecies by replacing specific monomers in a
        backbone polymer with bound versions according to the replacement
        dictionaries.

        Parameters
        ----------
        species_lists : list of dict
            List of replacement dictionaries where each dict maps positions
            (int) to [species, direction] pairs indicating which monomers
            to replace with bound versions.
        backbone : OrderedPolymerSpecies
            The base polymer species serving as the template for creating
            bound variants.

        Returns
        -------
        list of OrderedPolymerSpecies
            List of polymer species representing all specified combinatorial
            binding states.

        Notes
        -----
        This method takes a backbone polymer (typically the unbound construct)
        and creates variants where specific positions contain bound complexes.

        For example, given backbone `<A,B,C>` and replacements:

        - `{0: [A:RNAP, forward]}` creates `<[A:RNAP],B,C>`
        - `{0: [A:RNAP, forward], 1: [B:RNAP, forward]}` creates
          `<[A:RNAP],[B:RNAP],C>`

        """
        polymers = []
        for combo in species_lists:
            # to make combinatorially bound versions of a species, we take
            # the unbound species (backbone) and replace monomers of it
            # with bound versions (combo)
            #
            # combo is a dictionary of the form: {number:OrderedMonomer,...}
            # where the number is which element should be replaced and the
            # OrderedMonomer is what will be replacing that element
            new_backbone = OrderedPolymerSpecies.from_polymer_species(
                backbone, combo
            )
            polymers += [new_backbone]  # we make a new OrderedComplexSpecies
        return polymers

    def update_combinatorial_complexes(self, active_components):
        """Generate all combinatorial binding state species for the construct.

        Given components that can bind to different positions in the
        construct, this method generates all possible combinations of binding
        states by mixing and matching bound species at different locations.

        Parameters
        ----------
        active_components : list of Component
            Components that generate binding complexes with the construct.
            Each is assumed to bind at only one position.

        Returns
        -------
        list of OrderedPolymerSpecies
            All possible combinatorial binding states of the construct,
            including the unbound state and all single and multiple binding
            combinations.

        Notes
        -----
        The method:

        1. Collects all binary complex species from each active component
        2. Identifies unique positioned complexes
        3. Generates all combinatorial placements using `located_allcomb`
        4. Creates polymer species for each combination using `make_polymers`

        For example, given construct `<A,B,C>` with two components that
        create `<[A:RNAP],B,C>` and `<A,[B:RNAP],C>`, this method generates:

        - `<A,B,C>` (unbound)
        - `<[A:RNAP],B,C>` (A bound)
        - `<A,[B:RNAP],C>` (B bound)
        - `<[A:RNAP],[B:RNAP],C>` (both bound)

        assuming A and B act independently.

        """
        species = []
        for part in active_components:
            # first we make binary complexes
            sp_list = part.update_species()
            species += sp_list
        unique_complexes = []
        # species need to be uniqueified
        unique_species = list(set(species))
        for specie in unique_species:
            # in this list we extract all the variants of the complexes
            # from possible_backbones that exist in our species list.
            if (
                isinstance(specie, ComplexSpecies)
                and specie.position is not None
            ):
                unique_complexes += [specie]
        # unique_complexes now has a list of all the non-combinatorial
        # complexes we can make
        combinatorial_complexes = unique_complexes

        # all possible combinations of binders are made here
        allcomb = self.located_allcomb(combinatorial_complexes)
        allcomb += [{}]  # unbound dna should also be used

        # now, all possibilities have been enumerated.
        # we construct the OrderedPolymerSpecies
        out_polymers = self.make_polymers(allcomb, self.get_species())
        return out_polymers

    # Overwrite Component.enumerate_components
    def enumerate_constructs(self, previously_enumerated=None):
        """Run all enumerators to generate new construct variants.

        Applies all component enumerators to this construct to generate
        derived constructs (e.g., RNA_constructs from transcription).

        Parameters
        ----------
        previously_enumerated : set or list, optional
            Collection of constructs that have already been enumerated, used
            to prevent infinite recursion and duplicate enumeration.

        Returns
        -------
        list of Construct
            New constructs generated by all enumerators. For DNA_constructs
            with the default TxExplorer, this includes all possible
            RNA_construct transcripts.

        Notes
        -----
        Each enumerator's `enumerate_components` method is called with this
        construct and the previously enumerated set, allowing enumerators to
        explore transcriptional units, translational products, or other
        construct-derived components.

        See Also
        --------
        TxExplorer : Default enumerator for DNA transcription exploration.
        TlExplorer : Default enumerator for RNA translation exploration.

        """
        new_constructs = []
        for enumerator in self.component_enumerators:
            new_comp = enumerator.enumerate_components(
                component=self, previously_enumerated=previously_enumerated
            )
            new_constructs += new_comp
        return new_constructs

    def combinatorial_enumeration(self):
        """Generate components for all combinatorial binding states.

        Creates copies of parts that can react with different combinatorial
        binding states of the construct, ensuring reactions are generated for
        all possible binding configurations.

        Returns
        -------
        list of Component
            Components configured to react with different combinatorial
            binding states of the construct.

        Notes
        -----
        This method handles the generation of components that account for
        multiple simultaneous binding events. For example, given construct
        `<A,B,C>` where both A and B can bind RNAP:

        - Binary complexes: `<[A:RNAP],B,C>` and `<A,[B:RNAP],C>`
        - Combinatorial complex: `<[A:RNAP],[B:RNAP],C>`

        The method returns multiple versions of components A and B, each
        configured to bind to different pre-existing binding states:

        - A component binding to `<A,B,C>` --> `<[A:RNAP],B,C>`
        - A component binding to `<A,[B:RNAP],C>` --> `<[A:RNAP],[B:RNAP],C>`
        - B component binding to `<A,B,C>` --> `<A,[B:RNAP],C>`
        - B component binding to `<[A:RNAP],B,C>` --> `<[A:RNAP],[B:RNAP],C>`

        This ensures proper reaction enumeration for all binding
        combinations.

        """
        # Looks at combinatorial states of constructs to generate DNA_parts
        # my_polymer = self.get_species()
        self.update_parameters()

        # Go through parts
        active_components = []
        for part in self.parts_list:
            if hasattr(part, 'update_component'):
                dummy_species = self.get_species()[part.position]
                updated_components = part.update_component(
                    dummy_species, practice_run=True
                )
                if updated_components is not None:
                    active_components += [updated_components]
        # this next part creates "combinatorial" bound complexes, given
        # singly bound complexes generated above.  for example, let's say
        # we got <[A:RNAP],B,C> and <A,[B:RNAP],C> from A and B
        # individually, then a combinatorial construct would be
        # <[A:RNAP],[B:RNAP],C>
        combinatorial_complexes = self.update_combinatorial_complexes(
            active_components
        )

        combinatorial_components = []
        for comb_specie in combinatorial_complexes:
            # after making all combinatorially bound parts, we must seed
            # components with the right species so that the right reactions
            # are made. For example, A cannot react with
            # <[A:RNAP],[B:RNAP],C> because that's the species you get
            # after A has already reacted. So here we are finding species
            # with the proper positions unbound, so they can be fed to the
            # proper components
            if isinstance(comb_specie, OrderedPolymerSpecies):
                for part in active_components:
                    part_pos = part.position
                    if not isinstance(comb_specie[part_pos], ComplexSpecies):
                        # in this case the position of interest is not
                        # complexed. Thus, we need to update!
                        combinatorial_components += [
                            part.update_component(comb_specie[part_pos])
                        ]
        return combinatorial_components

    def enumerate_components(self, previously_enumerated=None):
        """Generate all derived components and constructs from this construct.

        Combines both construct enumeration (e.g., transcripts) and
        combinatorial component enumeration (for binding states) to produce
        a complete set of derived components.

        Parameters
        ----------
        previously_enumerated : set or list, optional
            Collection of components already enumerated, used to prevent
            infinite recursion.

        Returns
        -------
        list of Component
            All new components and constructs generated, including:

            - New constructs from enumerators (e.g., RNA_constructs)
            - Components for combinatorial binding states

        Notes
        -----
        This method generates new components in two scenarios:

        1. Binding-induced species: When a component creates a species
           that binds to part of the construct. For example, `<A,B,C>`
           --> `<[A:RNAP],B,C>` would return component A configured
           for this binding.

        2. Combinatorial binding states: When multiple components can
           bind simultaneously. For construct `<A,B,C>` where both A
           and B can bind RNAP:

           - Binary species: `<[A:RNAP],B,C>` and `<A,[B:RNAP],C>`
           - Combinatorial: `<[A:RNAP],[B:RNAP],C>`

           Returns components A and B configured to bind to various pre-bound
           states, ensuring all binding combinations are enumerated.

        3. New constructs: Generated by `enumerate_constructs()`. For example,
           a DNA_construct with promoter A generates an RNA_construct
           containing `<B,C>`.

        """
        # Runs component enumerator to generate new constructs
        new_constructs = self.enumerate_constructs(
            previously_enumerated=previously_enumerated
        )

        # Looks at combinatorial states of constructs to generate DNA_parts
        combinatorial_components = self.combinatorial_enumeration()

        return combinatorial_components + new_constructs

    @classmethod
    def get_partstring(cls, part):
        """Generate a string identifier for a part including its direction.

        Creates a unique string representation of a part that includes its
        name and direction but not its position, useful for construct
        comparison.

        Parameters
        ----------
        part : DNA_part or OrderedMonomer
            The part to generate a string identifier for.

        Returns
        -------
        str
            String combining the part's name and direction (e.g.,
            'promoter_pLacforward' or 'gene_GFPreverse').

        Notes
        -----
        This method creates an "orphan" copy of the part (without position
        or direction) to get the base name, then appends the direction to
        create a direction-aware identifier.

        """
        orphan = part.get_orphan()
        orphan.direction = None
        orphan.position = None
        curname = str(orphan)
        curdir = part.direction
        return curname + curdir

    @classmethod
    def get_partlist_hash(cls, partlist):
        """Generate a hash string for an ordered list of parts.

        Creates a unique string identifier for a parts list by concatenating
        part names with position indices, capturing the order and content of
        parts (but not their absolute positions).

        Parameters
        ----------
        partlist : list
            List of (part_string, part) tuples where part_string is typically
            generated by `get_partstring`.

        Returns
        -------
        str
            Hash string representing the ordered parts list.

        Notes
        -----
        The hash format concatenates each part's string representation with
        its index in the list, separated by underscores. This creates a
        unique identifier for the parts sequence.

        """
        partlist_str = '_'.join(
            [
                str(a[0]) + str(b)
                for a, b in zip(partlist, range(len(partlist)))
            ]
        )
        return partlist_str

    @classmethod
    def create_hashless_reverse(cls, construct):
        """Create a reversed construct without computing its hash.

        Generates a reversed version of the construct with parts in reverse
        order and flipped directions, but skips hash computation to avoid
        infinite recursion during hash calculations.

        Parameters
        ----------
        construct : Construct
            The construct to reverse.

        Returns
        -------
        Construct
            A new construct with reversed parts order and flipped directions,
            with `make_dirless_hash=False` to prevent hash computation.

        Notes
        -----
        This method is used internally by hash computation routines that need
        to compare forward and reverse orientations. Setting
        `make_dirless_hash=False` prevents infinite loops where hash
        computation would trigger reverse computation, which would trigger
        hash computation, etc.

        The circularity status of the construct is preserved.

        """
        rev_con = [a.get_orphan() for a in construct]
        for rev_part, origpart in zip(rev_con, construct):
            rev_part.direction = {'forward': 'reverse', 'reverse': 'forward'}[
                origpart.direction
            ]
        rev_con = rev_con[::-1]
        rev_con = Construct(
            rev_con, make_dirless_hash=False, circular=construct.circular
        )
        return rev_con

    @classmethod
    def rotation_free_hash(cls, construct):
        """Compute the most alphabetically ordered circular permutation hash.

        Finds the circular permutation of the construct that produces the
        most alphabetically ordered sequence of parts, providing a canonical
        representation for circular constructs regardless of starting
        position.

        Parameters
        ----------
        construct : Construct
            The circular construct to hash. Should have `circular=True`.

        Returns
        -------
        hash : str
            String hash of the most alphabetically ordered permutation.
        direction : int
            Always 1 (forward direction, since only rotation is considered).
        first_position : int
            The position that should be used as the first position to
            recreate this canonical permutation.

        Notes
        -----
        This method evaluates every possible starting position for a circular
        construct and selects the permutation that produces the most
        alphabetically ordered sequence when parts are compared
        lexicographically.

        To recreate the canonical form from the original construct:

        1. Rotate to start at `first_position`

        The direction is always 1 because this method only considers
        rotations, not reversals.

        Examples
        --------
        For a circular construct with parts A, B, C, if starting at C gives
        the most alphabetically ordered sequence (C, A, B), then
        `first_position` would be 2.

        """

        def circular_next(part, construct):
            if isinstance(part, int):
                part_pos = part
            else:
                part_pos = part.position + 1
            part_pos = part_pos % (len(construct))
            return construct[part_pos]

        best_partlist = [(None, None)]
        for part in construct:
            parthash = Construct.get_partstring(part)
            if best_partlist[0][0] is None:
                best_partlist = [(parthash, part)]
            elif parthash > best_partlist[0][0]:
                best_partlist = [(parthash, part)]
            elif parthash == best_partlist[0][0]:
                test_partlist = [(parthash, part)]
                testpos = part.position

                while (
                    test_partlist[-1][0]
                    == best_partlist[len(test_partlist) - 1][0]
                ):
                    if testpos - part.position >= len(construct):
                        # this means we tested every part and they would
                        # come in at the same alphabetical position
                        break
                    testpos += 1
                    next_testpart = circular_next(testpos, construct)
                    test_partlist += [
                        (
                            Construct.get_partstring(next_testpart),
                            next_testpart,
                        )
                    ]

                    if len(best_partlist) < len(test_partlist):
                        # this means we haven't made the hash string for
                        # the next part down the line
                        next_part = circular_next(
                            best_partlist[-1][1], construct
                        )
                        next_parthash = Construct.get_partstring(next_part)
                        best_partlist += [(next_parthash, next_part)]

                if (
                    test_partlist[-1][0]
                    > best_partlist[len(test_partlist) - 1][0]
                ):
                    best_partlist = test_partlist
        while len(best_partlist) < len(construct):
            nextpart = circular_next(best_partlist[-1][1], construct)
            best_partlist += [(Construct.get_partstring(nextpart), nextpart)]
        return (
            Construct.get_partlist_hash(best_partlist),
            1,
            best_partlist[0][1].position,
        )

    @classmethod
    def direction_rotation_free_hash(cls, construct):
        """Compute the best hash considering both rotation and direction.

        Finds the most alphabetically ordered representation of a circular
        construct by evaluating all rotations in both forward and reverse
        orientations.

        Parameters
        ----------
        construct : Construct
            The circular construct to hash.

        Returns
        -------
        hash : str
            String hash of the most alphabetically ordered permutation in
            either direction.
        direction : int
            Direction of the best ordering: 1 for forward, -1 for reverse.
        first_position : int
            The position that should be used as the first position in the
            best permutation.

        Notes
        -----
        This method:

        1. Computes the best forward rotation using `rotation_free_hash`
        2. Creates a reversed construct and computes its best rotation
        3. Returns whichever produces the more alphabetically ordered hash

        To recreate the canonical form:

        1. If direction is -1, reverse the construct
        2. Rotate to start at `first_position`

        This provides complete normalization for circular constructs,
        accounting for both rotation and direction symmetries.

        """
        rev_con = Construct.create_hashless_reverse(construct)
        forward_hash, _, pos = Construct.rotation_free_hash(construct)
        reverse_hash, _, posrev = Construct.rotation_free_hash(rev_con)
        if forward_hash > reverse_hash:
            return forward_hash, 1, pos
        else:
            return reverse_hash, -1, posrev

    @classmethod
    def linear_direction_free_hash(cls, construct):
        """Compute the best hash for a linear construct in either direction.

        Determines which orientation (forward or reverse) produces the most
        alphabetically ordered sequence for a linear construct.

        Parameters
        ----------
        construct : Construct
            The linear construct to hash.

        Returns
        -------
        hash : str
            String hash of the most alphabetically ordered orientation.
        direction : int
            Direction of the best ordering: 1 for forward, -1 for reverse.
        first_position : int
            Always 0 for linear constructs (no circular permutation).

        Notes
        -----
        This method compares forward and reverse orientations part-by-part,
        stopping as soon as one orientation is determined to be more
        alphabetically ordered than the other.

        To recreate the canonical form:

        1. If direction is -1, reverse the construct
        2. Start at position 0 (always the first position for linear
           constructs)

        Unlike circular constructs, linear constructs have a defined start
        and end, so the first_position is always 0.

        """
        rev_con = Construct.create_hashless_reverse(construct)

        test_partlist = []
        test_reverse_partlist = []
        winner = None
        for part, revpart in zip(construct, rev_con):
            # test both forward and reverse at the same time
            if winner is None:
                parthash = Construct.get_partstring(part)
                rev_parthash = Construct.get_partstring(revpart)
                if parthash > rev_parthash:
                    winner = 'forward'
                elif parthash < rev_parthash:
                    winner = 'reverse'
                test_partlist += [(parthash, part)]
                test_reverse_partlist += [(rev_parthash, revpart)]
            elif winner == 'forward':
                # as soon as one of the two directions comes out on top
                # just go with that one
                test_partlist += [(Construct.get_partstring(part), part)]
            elif winner == 'reverse':
                # as soon as one of the two directions comes out on top
                # just go with that one
                test_reverse_partlist += [
                    (Construct.get_partstring(revpart), revpart)
                ]
        if winner == 'forward':
            rethash = Construct.get_partlist_hash(test_partlist)
            return rethash, 1, 0
        elif winner == 'reverse':
            rethash = Construct.get_partlist_hash(test_reverse_partlist)
            return rethash, -1, 0

    @classmethod
    def omnihash(cls, construct):
        """Compute a canonical hash for the construct.

        Creates the most alphabetically ordered representation of a construct,
        accounting for direction (forward/reverse) and, for circular
        constructs, rotation. This provides a unique canonical identifier
        for functionally equivalent constructs.

        Parameters
        ----------
        construct : Construct
            The construct to hash (can be linear or circular).

        Returns
        -------
        hash : str
            Canonical hash string with 'circular' or 'linear' suffix
            indicating construct topology.
        direction : int
            Direction of the canonical form: 1 for forward, -1 for reverse.
        first_position : int
            Starting position for the canonical form. For circular constructs,
            this is the optimal rotation point. For linear constructs, always
            0.

        Notes
        -----
        For circular constructs:

        1. Evaluates all rotations in both orientations
        2. Selects the most alphabetically ordered permutation
        3. Appends 'circular' to the hash

        For linear constructs:

        1. Compares forward and reverse orientations
        2. Selects the most alphabetically ordered orientation
        3. Appends 'linear' to the hash

        To recreate the canonical form:

        1. If direction is -1, reverse the construct
        2. For circular constructs, rotate to start at `first_position`
        3. For linear constructs, `first_position` is always 0

        This hash enables identification of equivalent constructs that differ
        only in representation (rotation or direction).

        Examples
        --------
        Two circular constructs with the same parts in different rotations
        will produce the same omnihash:

        >>> A, B, C = (bcp.DNA_part(s) for s in ['A', 'B', 'C'])
        >>> construct1 = bcp.DNA_construct(
        ...     [[A, 'forward'], [B, 'forward'], [C, 'forward']],
        ...     circular=True)
        >>> construct2 = bcp.DNA_construct(
        ...     [[B, 'forward'], [C, 'forward'], [A, 'forward']],
        ...     circular=True)
        >>> hash1, _, _ = bcp.Construct.omnihash(construct1)
        >>> hash2, _, _ = bcp.Construct.omnihash(construct2)
        >>> hash1 == hash2
        True

        """
        if construct.circular:
            rhash, flip, posit = Construct.direction_rotation_free_hash(
                construct
            )
            return rhash + 'circular', flip, posit
        else:
            rhash, flip, posit = Construct.linear_direction_free_hash(
                construct
            )
            return rhash + 'linear', flip, posit

    def __hash__(self):
        return OrderedPolymer.__hash__(self)

    def __eq__(self, construct2):
        """Test equality between two constructs.

        Two constructs are considered equal if they have the same string
        representation and the same name.

        Parameters
        ----------
        construct2 : Construct
            The other construct to compare with.

        Returns
        -------
        bool
            True if constructs are equal, False otherwise.

        Notes
        -----
        This is a simple equality test based on string representation. It
        does not use deep comparison of parts or the direction-independent
        hash. For more sophisticated equivalence testing that accounts for
        rotations and reversals, use the `omnihash` method.

        See Also
        --------
        omnihash : Compute canonical hash accounting for rotation and
                   direction.

        """
        # TODO: make this be a python object comparison
        if self.__repr__() == construct2.__repr__() and (
            self.name == construct2.name
        ):
            return True
        else:
            return False

    def update_species(self):
        """Generate species for the construct.

        Returns
        -------
        list of Species
            List containing the construct's primary species representation.

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn()` to generate the species associated with this
        construct. For most constructs, this returns a single species
        representing the entire construct.

        """
        species = [self.get_species()]
        return species

    def reset_stored_data(self):
        """Clear all cached enumeration and prediction data.

        Resets the cached results from component enumeration, RNA prediction,
        and protein prediction, forcing these to be recomputed on the next
        access.

        Notes
        -----
        This method should be called whenever the construct is modified in a
        way that would invalidate cached data, such as:

        - Reversing the construct
        - Changing parts
        - Updating mechanisms or parameters

        The cached attributes that are reset:

        - `out_components`: Results from `enumerate_components`
        - `predicted_rnas`: List of predicted RNA products
        - `predicted_proteins`: List of predicted protein products

        """
        self.out_components = None
        self.predicted_rnas = None
        self.predicted_proteins = None

    def changed(self):
        """Handle construct changes by resetting caches and updating name.

        Called when the construct has been modified, this method resets all
        cached data and regenerates the construct's name to reflect its
        current state.

        Notes
        -----
        This method performs two operations:

        1. Resets all cached enumeration and prediction data via
           `reset_stored_data`
        2. Regenerates the construct name via `make_name` to ensure it
           reflects the current parts configuration

        This should be called after any structural modification to the
        construct.

        """
        self.reset_stored_data()
        self.name = self.make_name()

    def update_reactions(self, norna=False):
        """Generate reactions for the construct.

        Returns
        -------
        list of Reaction
            Empty list. Base `Construct` class does not generate reactions
            directly. Subclasses override this method to provide specific
            reaction generation.

        Parameters
        ----------
        norna : bool, default=False
            If True, RNA-related reactions are excluded (used in some
            subclass implementations).

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn()`. The base implementation returns an empty
        list because the base Construct class does not generate reactions
        directly - reactions are generated by the parts within the construct
        through their associated mechanisms.

        Subclasses like `DNA_construct` and `RNA_construct` may override this
        to provide construct-specific reaction generation.

        """
        return []


class DNA_construct(Construct, DNA):
    """DNA construct representing a functional genetic circuit.

    A DNA_construct is a specialized Construct for DNA sequences that can
    contain promoters, RBS sites, coding sequences, terminators, and other
    genetic elements. It supports transcription to generate RNA constructs and
    provides DNA-specific functionality. The class uses the 'transcription'
    mechanism to generate RNA products and related species/reactions during
    CRN compilation.

    Parameters
    ----------
    parts_list : list of list
        List of parts in format [[part, direction], ...] where each part
        must be a DNA_part or OrderedMonomer.
    name : str, optional
        Name of the DNA construct. If None, automatically generated.
    circular : bool, default=False
        If True, represents a circular DNA molecule (e.g., plasmid).
    mechanisms : dict or list, optional
        Custom mechanisms for this construct, overriding mixture defaults.
    parameters : dict, optional
        Parameter values specific to this construct.
    attributes : list of str, optional
        List of attribute tags for the construct.
    initial_concentration : float, optional
        Initial concentration of the DNA construct.
    copy_parts : bool, default=True
        If True, makes deep copies of parts when adding to construct.
    component_enumerators : list, optional
        List of enumerators for generating construct variants. Defaults to
        [TxExplorer()] which explores transcriptional variants.
    **kwargs
        Additional keyword arguments passed to parent constructors.

    Attributes
    ----------
    material_type : str
        Always 'dna' for DNA constructs.
    predicted_rnas : list or None
        Cached list of RNA_construct objects that can be transcribed.
    predicted_proteins : list or None
        Cached list of protein species that can be produced.

    See Also
    --------
    Construct : Base class for all constructs.
    RNA_construct : RNA version of constructs.
    DNA_part : Base class for DNA parts within constructs.
    TxExplorer : Default enumerator for transcriptional exploration.

    Notes
    -----
    DNA_constructs support several key features:

    - Transcription enumeration: Automatically identifies all possible
      transcripts based on promoter positions and orientations
    - Protein prediction: Predicts protein products from transcripts
      containing RBS sites
    - Circular DNA: Special handling for plasmids and other circular
      DNA molecules
    - Component enumeration: Generates functional variants based on
      the genetic parts present

    The default TxExplorer enumerator automatically explores all possible
    transcriptional units in the construct.

    Examples
    --------
    Create a simple gene expression construct:

    >>> promoter = bcp.Promoter('ptet')
    >>> rbs = bcp.RBS('RBS_standard')
    >>> cds = bcp.CDS('GFP')
    >>> terminator = bcp.Terminator('BBa_B0022')
    >>> parts = [
    ...     [promoter, 'forward'], [rbs, 'forward'],
    ...     [cds, 'forward'], [terminator, 'forward']
    ... ]
    >>> gene = bcp.DNA_construct(
    ...     parts_list=parts,
    ...     name='expression_cassette'
    ... )

    Create a circular plasmid:

    >>> ori = bcp.DNA_part('p15A')
    >>> plasmid_parts = [
    ...     [ori, 'forward'], [promoter, 'forward'], [rbs, 'forward'],
    ...     [gene, 'forward'], [terminator, 'forward']
    ... ]
    >>> plasmid = bcp.DNA_construct(
    ...     parts_list=plasmid_parts,
    ...     name='pUC19_GFP',
    ...     circular=True,
    ...     initial_concentration=10
    ... )

    """

    def __init__(
        self,
        parts_list,
        name=None,
        circular=False,
        mechanisms=None,  # custom mechanisms
        parameters=None,  # customized parameters
        attributes=None,
        initial_concentration=None,
        copy_parts=True,
        component_enumerators=None,
        **kwargs,
    ):
        self.material_type = 'dna'
        if component_enumerators is None:
            from ..construct_explorer import TxExplorer

            component_enumerators = [TxExplorer()]

        Construct.__init__(
            self=self,
            parts_list=parts_list,
            name=name,
            circular=circular,
            mechanisms=mechanisms,
            parameters=parameters,
            attributes=attributes,
            initial_concentration=initial_concentration,
            component_enumerators=component_enumerators,
            **kwargs,
        )
        DNA.__init__(self=self, name=self.name)

    def __repr__(self):
        return 'DNA_construct = ' + self.make_name()


class RNA_construct(Construct, RNA):
    """RNA construct representing a functional transcript.

    An RNA construct represents an RNA molecule that can be translated into
    proteins. Unlike DNA constructs, RNA constructs can only be linear (not
    circular) and primarily support translation rather than transcription.
    This class uses the 'translation' mechanism to generate protein products
    and related species/reactions during CRN compilation.

    Parameters
    ----------
    parts_list : list of list
        List of parts in format [[part, direction], ...] where parts
        represent functional RNA elements (RBS sites, coding sequences, etc.).
    name : str, optional
        Name of the RNA construct. If None, automatically generated.
    promoter : Promoter, optional
        Reference to the promoter that produced this RNA transcript.
    component_enumerators : list, optional
        List of enumerators for generating construct variants. Defaults to
        [TlExplorer()] which explores translational variants.
    length : int, default=0
        Length of the RNA molecule in nucleotides.
    **kwargs
        Additional keyword arguments passed to parent constructors.

    Attributes
    ----------
    material_type : str
        Always 'rna' for RNA constructs.
    promoter : Promoter or None
        The promoter that controls transcription of this RNA.
    length : int
        Length of the RNA transcript.
    circular : bool
        Always False for RNA (RNA cannot be circular).

    See Also
    --------
    Construct : Base class for all constructs.
    DNA_construct : DNA version of constructs.
    TlExplorer : Default enumerator for translational exploration.
    RNA : Base class for RNA components.

    Notes
    -----
    RNA_constructs have several key characteristics:

    - Linear only: RNA molecules cannot be circular
    - Translation focus: Primarily generates protein products through
      translation mechanisms
    - RBS enumeration: Automatically identifies all ribosome binding
      sites and potential translation products
    - No transcription: RNA cannot be transcribed to produce other RNA

    The default TlExplorer enumerator automatically explores all possible
    translational units in the RNA construct based on RBS positions.

    Examples
    --------
    Create an mRNA with RBS and coding sequence:

    >>> rbs1 = bcp.RBS('RBS1')
    >>> cds1 = bcp.CDS('GFP')
    >>> parts = [[rbs1, 'forward'], [cds1, 'forward']]
    >>> mrna = bcp.RNA_construct(
    ...     parts_list=parts,
    ...     name='mRNA_GFP'
    ... )

    Create a polycistronic mRNA with multiple RBS-CDS pairs:

    >>> rbs2 = bcp.RBS('RBS2')
    >>> cds2 = bcp.CDS('RFP')
    >>> strong_promoter = bcp.Promoter('pstrong')
    >>> parts = [
    ...     [rbs1, 'forward'], [cds1, 'forward'],
    ...     [rbs2, 'forward'], [cds2, 'forward']
    ... ]
    >>> polycistronic = bcp.RNA_construct(
    ...     parts_list=parts,
    ...     name='mRNA_operon',
    ...     promoter=strong_promoter
    ... )

    """

    def __init__(
        self,
        parts_list,
        name=None,
        promoter=None,
        component_enumerators=None,
        length=0,
        **kwargs,
    ):
        self.material_type = 'rna'
        self.promoter = promoter
        self.length = length
        if component_enumerators is None:
            from ..construct_explorer import TlExplorer

            component_enumerators = [TlExplorer()]

        Construct.__init__(
            self,
            parts_list=parts_list,
            circular=False,
            name=name,
            component_enumerators=component_enumerators,
            **kwargs,
        )
        RNA.__init__(self=self, name=self.name)

    def __repr__(self):
        # The name of an RNA should be different from DNA, right?
        return 'RNA_construct = ' + self.name


# DNA_part: a component-like intermediate class necessary for DNA_construct
# Author: Andrey Shur
# Latest update: 6/4/2020
#
class DNA_part(Component, OrderedMonomer):
    """Base class for individual DNA parts in constructs.

    A DNA_part represents a single functional genetic element (promoter, RBS,
    coding sequence, terminator, etc.) that can be assembled into larger
    DNA_constructs. Parts have position and direction within constructs and
    serve as the modular building blocks for synthetic genetic circuits.
    Unlike full Components, DNA_parts do not have initial concentrations -
    these must be set on the containing construct or assembly.

    Parameters
    ----------
    name : str
        Name of the DNA part.
    assembly : DNAassembly or OrderedPolymer, optional
        The assembly or construct containing this part.
    direction : str, optional
        Orientation of the part: 'forward' or 'reverse'.
    pos : int, optional
        Position of this part within its parent construct.
    sequence : str, optional
        DNA sequence of the part (for future sequence-level modeling).
    no_stop_codons : list, optional
        List of reading frames without stop codons. Used for identifying
        potential coding sequences. Default is empty list.
    material_type : str, default='part'
        Material classification for the part.
    **kwargs
        Additional keyword arguments passed to Component constructor.
        Note: 'initial_concentration' is not allowed for DNA_parts.

    Attributes
    ----------
    name : str
        Name of the part.
    sequence : str or None
        DNA sequence of the part.
    material_type : str
        Material classification ('part').
    no_stop_codons : list
        Reading frames without stop codons.
    assembly : DNAassembly or None
        Reference to containing assembly.
    position : int or None
        Position within parent construct.
    direction : str
        Orientation ('forward' or 'reverse').

    See Also
    --------
    Promoter : DNA_part for transcriptional control.
    RBS : DNA_part for translational control.
    DNA_construct : Container for multiple DNA_parts.
    OrderedMonomer : Base class for positioned elements.

    Raises
    ------
    AttributeError
        If 'initial_concentration' is provided (not allowed for DNA_parts).

    Notes
    -----
    DNA_parts are the fundamental building blocks of genetic constructs:

    - Modular: Can be reused in different constructs
    - Directional: Support forward and reverse orientations
    - Positional: Track their location within constructs
    - No concentration: Cannot have initial concentrations (only
      constructs/assemblies can)

    The no_stop_codons attribute is used to identify potential open reading
    frames for translation.

    Examples
    --------
    Create a generic DNA part:

    >>> part = bcp.DNA_part(
    ...     name='regulatory_element',
    ...     direction='forward'
    ... )

    Create a part with sequence information:

    >>> promoter_part = bcp.DNA_part(
    ...     name='pLac',
    ...     sequence='ATGCGATCG...',
    ...     direction='forward'
    ... )

    Use within a construct:

    >>> gene_part = bcp.DNA_part(
    ...    name='GFP',
    ...    sequence='TGAGTAAAGGAGAAGAA...',
    ...     direction='forward'
    ... )
    >>> parts = [[promoter_part, 'forward'], [gene_part, 'forward']]
    >>> construct = bcp.DNA_construct(parts_list=parts)

    """

    def __init__(
        self,
        name,
        assembly=None,
        direction=None,
        pos=None,
        sequence=None,
        no_stop_codons=[],
        material_type='part',
        **kwargs,
    ):
        # Modular component sequence.
        #
        # These get compiled into working components
        if 'initial_concentration' in kwargs:
            raise AttributeError(
                "DNA_part should not recieve initial_concentration keyword. "
                "Pass this into the DNAassembly or DNA_construct instead."
            )

        # Store/Process DNA part keywords
        self.sequence = sequence
        self.material_type = material_type

        # Most parts have stop codons. If you want your part to not have
        # stop codons, put "forward" and/or "reverse". Some parts will set
        # this, others won't. Default is that it has stop codons
        self.no_stop_codons = no_stop_codons
        if self.no_stop_codons is None:
            self.no_stop_codons = []

        Component.__init__(self=self, name=name, **kwargs)
        if isinstance(assembly, OrderedPolymer):
            OrderedMonomer.__init__(
                self, position=pos, parent=assembly, direction=direction
            )
        else:
            self.assembly = assembly
            OrderedMonomer.__init__(self, position=pos, direction=direction)

    @property
    def dna_species(self):
        """Species: The chemical species representation of this DNA part.

        Returns a Species object with material_type='part' representing this
        DNA_part as a chemical species in the CRN.

        """
        return Species(self.name, material_type='part')

    def __repr__(self):
        myname = self.name
        if self.position is not None:
            myname += '_' + str(self.position)
        if self.direction == 'reverse':
            myname += '_r'
        return myname

    def __hash__(self):
        return OrderedMonomer.__hash__(self) + hash(self.name)

    def __eq__(self, other):
        """Test equality between two DNA_parts.

        Parts are equal if they have the same type, name, parent assembly/
        construct, direction, and position.

        Parameters
        ----------
        other : DNA_part
            The other part to compare with.

        Returns
        -------
        bool
            True if parts are equal, False otherwise.

        Notes
        -----
        Equality requires matching:

        1. Type (both must be the same DNA_part subclass)
        2. Name (identical names)
        3. Assembly/parent (same parent construct or both have None)
        4. Direction (both forward or both reverse)
        5. Position (same position in parent construct)

        Parts are considered equal even if their parent constructs are
        different objects, as long as the string representations of the
        parents match.

        """
        if type(other) is type(self):
            if self.name == other.name:
                if self.assembly is not None and other.assembly is not None:
                    if str(self.assembly) == str(other.assembly):
                        return True
                elif self.assembly is not None or other.assembly is not None:
                    # if one has an assembly and the other doesn't, then
                    # they aren't the same!!
                    return False
                elif (self.parent is None and other.parent is None) or (
                    self.parent is not None
                    and other.parent is not None
                    and str(self.parent) == str(other.parent)
                ):
                    # this is for when we are using the OrderedMonomer for
                    # its intended function
                    if (
                        self.direction == other.direction
                        and self.position == other.position
                    ):
                        return True
        return False

    def clone(self, position, direction, parent_dna):
        """Attach this part to a specific position in a DNA construct.

        Parameters
        ----------
        position : int
            Position in the parent DNA where this part should be placed.
        direction : str
            Orientation of the part: 'forward' or 'reverse'.
        parent_dna : DNA_construct or OrderedPolymer
            The DNA construct that will contain this part.

        Returns
        -------
        DNA_part
            Returns self after setting position and parent.

        Notes
        -----
        This method establishes the relationship between a part and its
        containing construct, setting the part's position and orientation.

        """
        # Define where the part is in what piece of DNA.
        # TODO add warning if DNA_part is not cloned
        self.insert(parent_dna, position, direction)
        return self

    def unclone(self):
        """Remove this part from its parent construct.

        Detaches the part from any parent construct or assembly, resetting
        its position and parent references.

        Returns
        -------
        DNA_part
            Returns self after removal from parent.

        Notes
        -----
        This method calls the `remove` method from the `OrderedMonomer` base
        class to detach the part from its parent polymer structure.

        After calling this method, the part becomes "orphaned" and can be
        attached to a different construct using `clone`.

        See Also
        --------
        clone : Attach the part to a construct at a specific position.

        """
        self.remove()
        return self

    def reverse(self):
        """Reverse the orientation of this DNA part.

        Flips the direction of the part between 'forward' and 'reverse'.

        Returns
        -------
        DNA_part
            Returns self after reversing direction.

        Notes
        -----
        This method is typically called when a containing construct is
        reversed, ensuring all parts maintain proper relative orientation.

        """
        OrderedMonomer.reverse(self)
        return self
