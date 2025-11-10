"""Polymer support module.

The classes `OrderedPolymer` and `OrderedMonomer` are data structures used
to represent Polymers and their associated components.

These classes are used by `ChemicalReactionNetwork` species as well as certain
components such as DNA_construct.

"""

import copy
from warnings import warn


class MonomerCollection:
    """Collection of ordered monomers without any particular structure.

    A base container class that holds a collection of `OrderedMonomer`
    objects without imposing any ordering or structural constraints. This
    class serves as a parent class for more structured polymer
    representations like `OrderedPolymer`.

    Parameters
    ----------
    monomers : list of OrderedMonomer
        List of `OrderedMonomer` objects to include in the collection. Each
        monomer is copied and linked to this collection as its parent.

    Attributes
    ----------
    monomers : tuple of OrderedMonomer
        Tuple containing copies of the input monomers, with each monomer's
        parent set to this collection.

    See Also
    --------
    OrderedPolymer : A polymer with ordered monomers and directionality.
    OrderedMonomer : A unit that can belong to a MonomerCollection.

    Notes
    -----
    Monomers are stored as a tuple (immutable) to prevent direct
    modification. Each monomer is copied during initialization to ensure
    the collection maintains its own references.

    Examples
    --------
    Create a collection of monomers:

    >>> mon1 = bcp.OrderedMonomer()
    >>> mon2 = bcp.OrderedMonomer()
    >>> collection = bcp.MonomerCollection([mon1, mon2])
    >>> len(collection.monomers)
    2

    """

    def __init__(self, monomers):
        self.monomers = monomers

    @property
    def monomers(self):
        return self._monomers

    @monomers.setter
    def monomers(self, monomers):
        """Set the monomers in the collection.

        Parameters
        ----------
        monomers : list of OrderedMonomer
            List of `OrderedMonomer` objects to store in the collection.
            Each monomer is copied and has its parent set to this
            collection.

        Raises
        ------
        AssertionError
            If any element in `monomers` is not an `OrderedMonomer`.

        """
        mon_list = []
        for monomer in monomers:
            assert isinstance(monomer, OrderedMonomer)
            mon_copy = copy.copy(monomer)
            mon_copy.parent = self
            mon_list.append(mon_copy)
        self._monomers = tuple(mon_list)


class OrderedPolymer(MonomerCollection):
    """A polymer made up of OrderedMonomers with a specific order.

    Represents a linear sequence of monomers where each monomer has a
    defined position and direction. This class extends `MonomerCollection`
    to provide ordered, directional polymer structures commonly used to
    represent DNA constructs, RNA sequences, and protein chains.

    Parameters
    ----------
    parts : list or tuple
        Sequence of parts to add to the polymer. Each element can be:

        - An `OrderedMonomer` object (uses existing direction)
        - A list/tuple `[OrderedMonomer, direction]` specifying monomer
          and its direction explicitly

    default_direction : str or int, optional
        Default direction for monomers when not explicitly specified.
        Common values include 'forward', 'reverse', 0, 1, or None.

    Attributes
    ----------
    polymer : tuple of OrderedMonomer
        Ordered tuple of monomers in this polymer. Alias for `monomers`
        property inherited from `MonomerCollection`.
    default_direction : str, int, or None
        Default direction assigned to monomers lacking explicit direction.

    See Also
    --------
    NamedPolymer : An OrderedPolymer with an associated name.
    OrderedMonomer : A unit that belongs to an OrderedPolymer.
    MonomerCollection : Base class for monomer collections.

    Notes
    -----
    Directions indicate the orientation of monomers in the polymer:

    - 'forward' or 0: Standard/positive orientation
    - 'reverse' or 1: Inverted/negative orientation
    - None: No specified orientation

    The polymer tuple is immutable, but monomers can be added via
    `insert`, `append`, or `replace` methods. Direct assignment to
    positions uses `__setitem__` which calls `replace`.

    All monomers are deep-copied when added to ensure the polymer
    maintains independent references. This prevents external modifications
    from affecting the polymer structure.

    Examples
    --------
    Create a polymer from monomers:

    >>> mon1 = bcp.OrderedMonomer()
    >>> mon2 = bcp.OrderedMonomer()
    >>> polymer = bcp.OrderedPolymer(
    ...     parts=[mon1, mon2],
    ...     default_direction='forward'
    ... )
    >>> len(polymer)
    2

    Create a polymer with explicit directions:

    >>> polymer = bcp.OrderedPolymer(
    ...     parts=[[mon1, 'forward'], [mon2, 'reverse']]
    ... )
    >>> polymer[0].direction
    'forward'
    >>> polymer[1].direction
    'reverse'

    """

    def __init__(self, parts, default_direction=None):
        self.default_direction = default_direction
        self.polymer = parts

    @property
    def polymer(self):
        return self._monomers

    @polymer.setter
    def polymer(self, parts):
        """Set the polymer sequence from a list of parts.

        Parameters
        ----------
        parts : list or tuple
            Sequence of parts to add to the polymer. Each element can be:

            - An `OrderedMonomer` object
            - A list/tuple `[OrderedMonomer, direction]`

        Raises
        ------
        AssertionError
            If `parts` is not a list or tuple.
        ValueError
            If any element is not an `OrderedMonomer` or valid part
            specification.

        """
        polymer = []
        assert isinstance(
            parts, (list, tuple)
        ), 'OrderedPolymer must be instantiated with a list'
        for item in parts:
            if isinstance(item, (list, tuple)):
                part = item[0]
                if len(item) > 1:
                    partdir = item[1]
                else:
                    partdir = None
            elif isinstance(item, OrderedMonomer):
                part = item
                partdir = item.direction
            else:
                raise ValueError(
                    f"{str(item)} is not an OrderedMonomer or "
                    "a list of the form [OrderedMonomer,direction]"
                )

            # OrderedMonomers are always copied when inserted into an
            # OrderedPolymer
            part_copy = copy.copy(part)
            polymer += [part_copy]
            position = len(polymer) - 1
            if partdir is None:
                partdir = self.default_direction
            part_copy.monomer_insert(self, position, partdir)

        self._monomers = tuple(polymer)

    def __hash__(self):
        hval = 0
        if not hasattr(self, '_polymer') or len(self._polymer) == 0:
            hval = 0
        else:
            hval = sum([a.subhash() for a in self._polymer])
        if hasattr(self, 'name'):
            hval += hash(self.name)

        return hval

    def changed(self):
        """Callback method invoked whenever the polymer structure changes.

        This method is called after operations that modify the polymer,
        such as `insert`, `replace`, `delpart`, or `reverse`.
        Subclasses can override this to implement custom behavior when
        the polymer is modified.

        Notes
        -----
        The base implementation does nothing. Override in subclasses to add
        functionality like name regeneration, validation, or notifications.

        """
        # runs whenever anything changed
        pass

    def insert(self, position, part, direction=None):
        """Insert a monomer at a specific position in the polymer.

        Inserts a copy of the given monomer at the specified position,
        shifting all subsequent monomers to higher positions. Calls the
        `changed` callback after insertion.

        Parameters
        ----------
        position : int
            Index at which to insert the monomer. Must be between 0 and
            `len(polymer)` (inclusive).
        part : OrderedMonomer
            The monomer to insert. A copy of this monomer will be added.
        direction : str, int, or None, optional
            Direction for the inserted monomer. If None, uses the monomer's
            existing direction.

        See Also
        --------
        append : Add a monomer to the end of the polymer.
        replace : Replace a monomer at a specific position.

        Notes
        -----
        The monomer is deep-copied before insertion to maintain
        independence from the original object.

        """
        # OrderedMonomers are always copied when inserted into an
        # OrderedPolymer
        part_copy = copy.copy(part)

        if direction is None:
            direction = part.direction

        part_copy.monomer_insert(self, position, direction)
        for subsequent_part in self.polymer[position:]:
            subsequent_part.position += 1
        self.polymer = (
            self.polymer[:position] + (part_copy,) + self.polymer[position:]
        )
        self.changed()

    def replace(self, position, part, direction=None):
        """Replace a monomer at a specific position in the polymer.

        Removes the monomer at the given position and inserts a copy of
        the new monomer in its place. Calls the `changed` callback after
        replacement.

        Parameters
        ----------
        position : int
            Index of the monomer to replace. Must be a valid position in
            the polymer.
        part : OrderedMonomer
            The monomer to insert. A copy of this monomer will be added.
        direction : str, int, or None, optional
            Direction for the new monomer. If None, uses the monomer's
            existing direction.

        See Also
        --------
        insert : Insert a monomer at a specific position.
        delpart : Remove a monomer from the polymer.

        Notes
        -----
        The removed monomer's `remove` method is called to clear its
        parent, position, and direction attributes.

        """
        # OrderedMonomers are always copied when inserted into an
        # OrderedPolymer
        part_copy = copy.copy(part)

        if direction is None:
            direction = part.direction

        self.polymer[position].remove()
        part_copy.monomer_insert(self, position, direction)
        self.polymer = (
            self.polymer[:position]
            + (part_copy,)
            + self.polymer[position + 1 :]
        )
        self.changed()

    def append(self, part, direction=None):
        """Add a monomer to the end of the polymer.

        Appends a copy of the given monomer to the end of the polymer
        sequence by calling `insert` at the final position.

        Parameters
        ----------
        part : OrderedMonomer
            The monomer to append. A copy of this monomer will be added.
        direction : str, int, or None, optional
            Direction for the appended monomer. If None, uses the monomer's
            existing direction if available.

        See Also
        --------
        insert : Insert a monomer at a specific position.

        Examples
        --------
        >>> polymer = bcp.OrderedPolymer(parts=[])
        >>> mon = bcp.OrderedMonomer()
        >>> polymer.append(mon, direction='forward')
        >>> len(polymer)
        1

        """
        # OrderedMonomers are always copied when inserted into an
        # OrderedPolymer
        part_copy = copy.copy(part)

        if direction is None:
            if hasattr(part, 'direction'):
                direction = part_copy.direction
            else:
                direction = None
        pos = len(self.polymer)
        self.insert(pos, part_copy, direction)

    def __repr__(self):
        outstr = 'polymer('
        for part in self.polymer:
            outstr += str(part) + ', direction = ' + str(part.direction) + ','
        if outstr[:-1] == ',':
            outstr = outstr[:-1]
        outstr += ')'
        return outstr

    def direction_invert(self, dirname):
        """Invert a direction value.

        Converts a direction to its opposite orientation. Used during
        polymer reversal operations.

        Parameters
        ----------
        dirname : str, int, or None
            The direction to invert. Supported values:

            - 'forward' <--> 'reverse'
            - 0 <--> 1
            - None -> None

        Returns
        -------
        str, int, or None
            The inverted direction. Returns the input unchanged if it
            cannot be inverted.

        Warns
        -----
        UserWarning
            If the direction value is not recognized.

        Examples
        --------
        >>> polymer = bcp.OrderedPolymer(parts=[])
        >>> polymer.direction_invert('forward')
        'reverse'
        >>> polymer.direction_invert(0)
        1

        """
        if dirname == 'forward':
            return 'reverse'
        elif dirname == 'reverse':
            return 'forward'
        elif dirname == 0:
            return 1
        elif dirname == 1:
            return 0
        elif dirname is None:
            return None
        else:
            warn("didn't know how to invert {}".format(str(dirname)))
            return dirname

    def __len__(self):
        """Return the number of monomers in the polymer.

        Returns
        -------
        int
            The number of monomers in the polymer sequence.

        """
        return len(self.polymer)

    def __getitem__(self, ii):
        """Get a monomer or slice of monomers from the polymer.

        Parameters
        ----------
        ii : int or slice
            Index or slice to retrieve from the polymer.

        Returns
        -------
        OrderedMonomer or tuple
            The monomer at the given index, or a tuple of monomers for a
            slice.

        """
        return self.polymer[ii]

    def __setitem__(self, ii, val):
        """Replace a monomer at a specific position.

        Parameters
        ----------
        ii : int
            Index at which to replace the monomer.
        val : OrderedMonomer
            The new monomer to insert at the position.

        Notes
        -----
        Internally calls `replace` with the monomer's existing direction.

        """
        self.replace(ii, val, val.direction)

    def __eq__(self, other):
        """Check equality with another OrderedPolymer.

        Two polymers are equal if they have the same length and each
        corresponding pair of monomers has the same direction, position,
        and type.

        Parameters
        ----------
        other : OrderedPolymer
            The polymer to compare with.

        Returns
        -------
        bool
            True if polymers are equal, False otherwise.

        """
        if isinstance(other, OrderedPolymer):
            for item1, item2 in zip(self.polymer, other.polymer):
                if (
                    item1.direction == item2.direction
                    and item1.position == item2.position
                    and type(item1) is type(item2)
                ):
                    pass
                else:
                    return False
            if len(self.polymer) == len(other.polymer):
                return True
        return False

    def __contains__(self, item):
        """Check if a monomer is in the polymer.

        Parameters
        ----------
        item : OrderedMonomer
            The monomer to search for.

        Returns
        -------
        bool
            True if the monomer is in the polymer, False otherwise.

        """
        if item in self.polymer:
            return True
        else:
            return False

    def delpart(self, position):
        """Remove a monomer from the polymer at a specific position.

        Removes the monomer at the given position, shifts all subsequent
        monomers to lower positions, and calls the `changed` callback.

        Parameters
        ----------
        position : int
            Index of the monomer to remove. Must be a valid position in
            the polymer.

        See Also
        --------
        replace : Replace a monomer at a specific position.
        insert : Insert a monomer at a specific position.

        Notes
        -----
        The removed monomer's `remove` method is called to clear its
        parent, position, and direction. If the polymer has a `name`
        attribute and a `make_name` method, the name is regenerated
        after deletion.

        """
        part = self.polymer[position]
        part.remove()
        for subsequent_part in self.polymer[position + 1 :]:
            subsequent_part.position -= 1
        self.polymer = self.polymer[:position] + self.polymer[position + 1 :]
        self.changed()
        if hasattr(self, 'name') and hasattr(self, 'make_name'):
            self.name = self.make_name()

    def reverse(self):
        """Reverse the order and directions of all monomers in the polymer.

        Reverses the polymer sequence and inverts the direction of each
        monomer. Updates all monomer positions to reflect their new
        locations. Calls the `changed` callback after reversal.

        Notes
        -----
        This operation modifies the polymer in place. All monomers have
        their directions inverted using `direction_invert` and their
        positions updated to match the reversed sequence.

        Examples
        --------
        >>> mon1 = bcp.OrderedMonomer()
        >>> mon2 = bcp.OrderedMonomer()
        >>> polymer = bcp.OrderedPolymer(
        ...     parts=[[mon1, 'forward'], [mon2, 'reverse']]
        ... )
        >>> polymer.reverse()
        >>> polymer[0].direction
        'forward'
        >>> polymer[1].direction
        'reverse'

        """
        self.polymer = self.polymer[::-1]
        for ind, part in enumerate(self.polymer):
            part.position = ind
            part.direction = self.direction_invert(part.direction)
        self.changed()


class NamedPolymer(OrderedPolymer):
    """An OrderedPolymer with an associated name and circularity flag.

    Extends `OrderedPolymer` to include a name identifier and optional
    circular topology flag. Commonly used to represent named biological
    constructs like plasmids, chromosomes, or specific DNA/RNA sequences.

    Parameters
    ----------
    parts : list or tuple
        Sequence of parts to add to the polymer. See `OrderedPolymer` for
        format details.
    name : str
        Name identifier for the polymer.
    default_direction : str, int, or None, optional
        Default direction for monomers when not explicitly specified.
    circular : bool, default=False
        If True, indicates the polymer has circular topology (e.g., a
        plasmid). If False, the polymer is linear.

    Attributes
    ----------
    name : str
        Name identifier for the polymer.
    circular : bool
        Flag indicating circular (True) or linear (False) topology.
    polymer : tuple of OrderedMonomer
        Ordered tuple of monomers in this polymer (inherited).
    default_direction : str, int, or None
        Default direction for monomers (inherited).

    See Also
    --------
    OrderedPolymer : Base class for ordered polymer structures.
    OrderedMonomer : A unit that belongs to an OrderedPolymer.

    Notes
    -----
    The `circular` attribute is primarily informational and does not
    automatically enforce circular topology constraints in polymer
    operations. Subclasses or external code must handle circular semantics
    as needed.

    Examples
    --------
    Create a linear named polymer:

    >>> mon1 = bcp.OrderedMonomer()
    >>> mon2 = bcp.OrderedMonomer()
    >>> polymer = bcp.NamedPolymer(
    ...     parts=[mon1, mon2],
    ...     name='my_construct',
    ...     default_direction='forward'
    ... )
    >>> polymer.name
    'my_construct'
    >>> polymer.circular
    False

    Create a circular polymer (plasmid):

    >>> plasmid = bcp.NamedPolymer(
    ...     parts=[mon1, mon2],
    ...     name='pUC19',
    ...     circular=True
    ... )
    >>> plasmid.circular
    True

    """

    def __init__(self, parts, name, default_direction=None, circular=False):
        self.name = name
        self.circular = circular
        OrderedPolymer.__init__(
            self=self, parts=parts, default_direction=default_direction
        )


class OrderedMonomer:
    """A unit that belongs to an OrderedPolymer.

    Represents a single monomer unit within a polymer structure. Each
    monomer tracks its position in the polymer, its directional
    orientation, and maintains a reference to its parent polymer. This
    class is used as a base for representing DNA parts, RNA components,
    amino acids, and other polymer building blocks.

    Parameters
    ----------
    direction : str, int, or None, optional
        Directional orientation of the monomer in the polymer. Common
        values include 'forward', 'reverse', 0, 1, or None. Default is
        None.
    position : int or None, optional
        Index position of the monomer within its parent polymer. Must be
        non-None if the monomer belongs to a polymer. Default is None.
    parent : MonomerCollection or None, optional
        Reference to the parent `MonomerCollection` or `OrderedPolymer`
        containing this monomer. Default is None.

    Attributes
    ----------
    direction : str, int, or None
        Directional orientation of the monomer.
    position : int or None
        Position index within the parent polymer.
    parent : MonomerCollection or None
        Reference to the parent collection or polymer.
    is_polymer_component : bool
        Flag indicating whether this monomer is part of a polymer
        structure. Set to True when inserted into a polymer.

    See Also
    --------
    OrderedPolymer : A polymer made up of OrderedMonomers.
    MonomerCollection : Base class for monomer collections.

    Notes
    -----
    OrderedMonomers are deep-copied when inserted into polymers to ensure
    independence. Use `get_orphan` to create a copy without a parent
    reference, or `get_removed` to create a fully detached copy.

    The `is_polymer_component` flag and `find_polymer_component` method
    support scenarios where monomers may be nested within complex species.

    Examples
    --------
    Create a standalone monomer:

    >>> monomer = bcp.OrderedMonomer(direction='forward')
    >>> monomer.direction
    'forward'
    >>> monomer.parent is None
    True

    Add a monomer to a polymer:

    >>> polymer = bcp.OrderedPolymer(parts=[])
    >>> monomer = bcp.OrderedMonomer()
    >>> polymer.append(monomer, direction='forward')
    >>> polymer[0].position
    0
    >>> polymer[0].parent is polymer
    True

    """

    def __init__(self, direction=None, position=None, parent=None):
        """Initialize an OrderedMonomer.

        Parameters
        ----------
        direction : str, int, or None, optional
            Directional orientation of the monomer. Default is None.
        position : int or None, optional
            Position index within the parent polymer. Default is None.
        parent : MonomerCollection or None, optional
            Reference to the parent collection. Default is None.

        """
        # The default is that the monomer is not part of a polymer.
        self.parent = None
        self.direction = None
        # Set position to prevent weird testing errors of not having
        # attributes
        self.position = None
        # by default, we assume that an orderedmonomer is not part of a
        # polymer
        self.is_polymer_component = False
        # Set properties correctly
        self.parent = parent
        self.direction = direction
        self.position = position

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        """Set the parent collection for this monomer.

        Parameters
        ----------
        parent : MonomerCollection or None
            The parent collection or polymer to assign to this monomer.

        Raises
        ------
        ValueError
            If `parent` is not None and not a `MonomerCollection` instance.

        """
        if parent is None or isinstance(parent, MonomerCollection):
            self._parent = parent
        else:
            raise ValueError(
                f"parent must be an MonomerCollection. Recieved {parent}"
            )

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, direction):
        """Set the directional orientation of the monomer.

        Parameters
        ----------
        direction : str, int, or None
            The direction to assign. Common values include 'forward',
            'reverse', 0, 1, or None.

        """
        self._direction = direction

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        """Set the position index of the monomer.

        Parameters
        ----------
        position : int or None
            The position index to assign. Must be non-None if the monomer
            has a parent polymer.

        Raises
        ------
        ValueError
            If the monomer has a parent but position is None.

        """
        if self.parent is not None and position is None:
            raise ValueError(f"{self} is part of a polymer with no position!")
        else:
            self._position = position

    def find_polymer_component(self):
        """Find the polymer component within this monomer or its species.

        Searches this monomer and, if it is a `ComplexSpecies`, its
        constituent species to find which one is marked as a polymer
        component.

        Returns
        -------
        OrderedMonomer or None
            The monomer that is part of a polymer structure, or None if no
            polymer component is found.

        Raises
        ------
        ValueError
            If multiple species are marked as polymer components in the
            same location.

        Notes
        -----
        This method is primarily used internally to handle complex species
        that may contain monomers as part of larger structures.

        """
        from .species import ComplexSpecies

        outpolymer = None
        if isinstance(self, ComplexSpecies):
            for specie in self.species:
                if specie.is_polymer_component:
                    if outpolymer is not None:
                        raise ValueError(
                            "multiple species are part of the polymer "
                            "in the same place!!"
                        )
                    else:
                        outpolymer = specie
        if self.is_polymer_component:
            if outpolymer is not None:
                raise ValueError(
                    "multiple species are part of the polymer "
                    "in the same place!!"
                )
            else:
                outpolymer = self
        return outpolymer

    def monomer_insert(
        self, parent: OrderedPolymer, position: int, direction=None
    ):
        """Insert this monomer into a polymer at a specific position.

        Sets the monomer's parent, position, and direction attributes to
        reflect its insertion into the polymer. Marks the monomer (or its
        polymer component if it is a complex species) as a polymer
        component.

        Parameters
        ----------
        parent : OrderedPolymer
            The polymer to insert this monomer into.
        position : int
            The position index where this monomer is being inserted.
        direction : str, int, or None, optional
            The direction for this monomer. If None, uses the monomer's
            existing direction.

        Raises
        ------
        ValueError
            If position is None, or if parent is None.

        """
        if position is None:
            raise ValueError(f"{self} has no position to be inserted at!")
        if direction is None:
            if self.direction is not None:
                direction = self.direction
        if parent is None:
            raise ValueError(f"{self} is trying to be inserted into nothing!")
        if self.is_polymer_component is False:
            if self.find_polymer_component() is None:
                self.is_polymer_component = True
        self.parent = parent
        self.position = position
        self.direction = direction

    def set_dir(self, direction):
        """Set the direction of the monomer and return self.

        Convenience method for setting direction in a fluent interface
        style.

        Parameters
        ----------
        direction : str, int, or None
            The direction to assign to this monomer.

        Returns
        -------
        OrderedMonomer
            Returns self for method chaining.

        Examples
        --------
        >>> monomer = bcp.OrderedMonomer().set_dir('forward')
        >>> monomer.direction
        'forward'

        """
        self.direction = direction
        return self

    def remove(self):
        """Remove this monomer from its parent polymer.

        Clears the monomer's parent, position, and direction attributes,
        effectively detaching it from any polymer structure.

        Returns
        -------
        OrderedMonomer
            Returns self for method chaining.

        See Also
        --------
        get_removed : Create a fully detached copy of the monomer.
        get_orphan : Create a copy with parent removed but position and
                     direction preserved.

        """
        self.parent = None
        self.position = None
        self.direction = None

        return self

    def get_orphan(self):
        """Create a copy of this monomer without a parent reference.

        Returns a copy that retains position and direction but has no
        parent polymer. Useful for temporarily working with monomers
        outside their polymer context.

        Returns
        -------
        OrderedMonomer
            A copy of this monomer with parent set to None but position
            and direction preserved.

        See Also
        --------
        get_removed : Create a fully detached copy.
        remove : Remove this monomer from its parent in place.

        Notes
        -----
        This is a shallow copy of the monomer object itself, though the
        parent reference is explicitly cleared.

        """
        # Returns a copy of this monomer, except with no parent.
        # But it still has a position and direction.
        copied_monomer = copy.copy(self)
        copied_monomer.parent = None
        return copied_monomer

    def get_removed(self):
        """Create a fully detached copy of this monomer.

        Returns a copy with all polymer-related attributes (parent,
        position, direction) cleared. Also removes 'forward' and 'reverse'
        attributes if present.

        Returns
        -------
        OrderedMonomer
            A copy of this monomer with no parent, position, or direction,
            and with directional attributes removed.

        See Also
        --------
        get_orphan : Create a copy without parent but with position and
                     direction.
        remove : Remove this monomer from its parent in place.

        Notes
        -----
        This method is useful for creating completely independent copies of
        monomers that can be reused in different contexts without any
        polymer associations.

        """
        copied_part = copy.copy(self)
        copied_part.parent = None
        copied_part.direction = None
        copied_part.position = None
        if hasattr(copied_part, '_attributes'):
            copied_part.remove_attribute('forward')
            copied_part.remove_attribute('reverse')
        return copied_part

    def __repr__(self):
        txt = (
            'OrderedMonomer(direction='
            + str(self.direction)
            + ',position='
            + str(self.position)
            + ')'
        )
        return txt

    def __eq__(self, other):
        """Check equality with another OrderedMonomer.

        Two monomers are equal if they have the same direction, position,
        and parent.

        Parameters
        ----------
        other : OrderedMonomer
            The monomer to compare with.

        Returns
        -------
        bool
            True if monomers are equal, False otherwise.

        """
        if isinstance(other, OrderedMonomer):
            if (
                self.direction == other.direction
                and self.position == other.position
                and self.parent == other.parent
            ):
                return True
        return False

    def __hash__(self):
        """Compute hash value for this monomer.

        Returns
        -------
        int
            Hash value based on position, direction, name (if present),
            and parent.

        """
        hval = 0
        hval += self.subhash()

        if self.parent is not None:
            hval += hash(self.parent)

        return hval

    def subhash(self):
        """Compute hash contribution from monomer properties.

        Computes a hash value based on the monomer's position, direction,
        and name (if present), excluding the parent reference.

        Returns
        -------
        int
            Hash value based on monomer-specific properties.

        Notes
        -----
        This method is used by `__hash__` to compute the monomer's hash
        contribution. It excludes the parent to avoid circular dependencies
        in hash computation.

        """
        hval = 0
        hval += hash(self.position)
        hval += hash(self.direction)
        if hasattr(self, 'name'):
            hval += hash(self.name)
        return hval
