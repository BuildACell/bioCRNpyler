"""The classes OrderedPolymer and OrderedMonomer are datastructures used to represent Polymers and their associatd components.

These classes are used by Chemical Reaction Network Species as well as certain Components such as DNA_construct.
"""
import copy
from warnings import warn


class MonomerCollection:
    """
    A class used to represent a collection of OrderedMonomers without any particular structure
    """

    def __init__(self, monomers):
        self.monomers = monomers

    @property
    def monomers(self):
        return self._monomers

    @monomers.setter
    def monomers(self, monomers):
        mon_list = []
        for monomer in monomers:
            assert isinstance(monomer, OrderedMonomer)
            mon_copy = copy.copy(monomer)
            mon_copy.parent = self
            mon_list.append(mon_copy)
        self._monomers = tuple(mon_list)


class OrderedPolymer(MonomerCollection):

    """a polymer made up of OrderedMonomers that has a specific order"""

    def __init__(self, parts, default_direction=None):
        """parts can be a list of lists containing 
        [[OrderedMonomer,direction],[OrderedMonomer,direction],...]
        alternatively, you can have a regular list, and the direcitons
        will end up being None"""
        self.default_direction = default_direction
        self.polymer = parts

    @property
    def polymer(self):
        return self._monomers

    @polymer.setter
    def polymer(self, parts):
        polymer = []
        assert (type(parts) == list or type(parts) ==
                tuple), "OrderedPolymer must be instantiated with a list"
        for item in parts:
            if (isinstance(item, list) or isinstance(item, tuple)):
                part = item[0]
                if (len(item) > 1):
                    partdir = item[1]
                else:
                    partdir = None
            elif (isinstance(item, OrderedMonomer)):
                part = item
                partdir = item.direction
            else:
                raise ValueError(
                    "{} is not an OrderedMonomer or a list of the form [OrderedMonomer,direction]".format(str(item)))
            # OrderedMonomers are always copied when inserted into an OrderedPolymer
            part_copy = copy.copy(part)
            polymer += [part_copy]
            position = len(polymer)-1
            if (partdir == None):
                partdir = self.default_direction
            part_copy.monomer_insert(self, position, partdir)

        self._monomers = tuple(polymer)

    def __hash__(self):
        hval = 0
        if (not hasattr(self, "_polymer") or len(self._polymer) == 0):
            hval = 0
        else:
            hval = sum([a.subhash() for a in self._polymer])
        if (hasattr(self, "name")):
            hval += hash(self.name)

        return hval

    def changed(self):
        # runs whenever anything changed
        pass

    def insert(self, position, part, direction=None):
        # OrderedMonomers are always copied when inserted into an OrderedPolymer
        part_copy = copy.copy(part)

        if (direction is None):
            direction = part.direction

        part_copy.monomer_insert(self, position, direction)
        for subsequent_part in self.polymer[position:]:
            subsequent_part.position += 1
        self.polymer = self.polymer[:position] + \
            (part_copy,)+self.polymer[position:]
        self.changed()

    def replace(self, position, part, direction=None):
        # OrderedMonomers are always copied when inserted into an OrderedPolymer
        part_copy = copy.copy(part)

        if (direction is None):
            direction = part.direction

        self.polymer[position].remove()
        part_copy.monomer_insert(self, position, direction)
        self.polymer = self.polymer[:position] + \
            (part_copy,)+self.polymer[position+1:]
        self.changed()

    def append(self, part, direction=None):
        # OrderedMonomers are always copied when inserted into an OrderedPolymer
        part_copy = copy.copy(part)

        if (direction is None):
            if (hasattr(part, "direction")):
                direction = part_copy.direction
            else:
                direction = None
        pos = len(self.polymer)
        self.insert(pos, part_copy, direction)

    def __repr__(self):
        outstr = "polymer("
        for part in self.polymer:
            outstr += str(part)+", direction = "+str(part.direction)+","
        if (outstr[:-1] == ","):
            outstr = outstr[:-1]
        outstr += ")"
        return outstr

    def direction_invert(self, dirname):
        if (dirname == "forward"):
            return "reverse"
        elif (dirname == "reverse"):
            return "forward"
        elif (dirname == 0):
            return 1
        elif (dirname == 1):
            return 0
        elif (dirname is None):
            return None
        else:
            warn("didn't know how to invert {}".format(str(dirname)))
            return dirname

    def __len__(self):
        return len(self.polymer)

    def __getitem__(self, ii):
        return self.polymer[ii]

    def __setitem__(self, ii, val):
        self.replace(ii, val, val.direction)

    def __eq__(self, other):
        if (isinstance(other, OrderedPolymer)):
            for item1, item2 in zip(self.polymer, other.polymer):
                if (item1.direction == item2.direction and item1.position == item2.position and type(item1) == type(item2)):
                    pass
                else:
                    return False
            if (len(self.polymer) == len(other.polymer)):
                return True
        return False

    def __contains__(self, item):
        if (item in self.polymer):
            return True
        else:
            return False

    def delpart(self, position):
        part = self.polymer[position]
        part.remove()
        for subsequent_part in self.polymer[position+1:]:
            subsequent_part.position -= 1
        self.polymer = self.polymer[:position] + self.polymer[position+1:]
        self.changed()
        if (hasattr(self, "name") and hasattr(self, "make_name")):
            self.name = self.make_name()

    def reverse(self):
        self.polymer = self.polymer[::-1]
        for ind, part in enumerate(self.polymer):
            part.position = ind
            part.direction = self.direction_invert(part.direction)
        self.changed()


class NamedPolymer(OrderedPolymer):
    """The same as an OrderedPolymer but it has a name"""

    def __init__(self, parts, name, default_direction=None, circular=False):
        self.name = name
        self.circular = circular
        OrderedPolymer.__init__(self=self, parts=parts,
                                default_direction=default_direction)


class OrderedMonomer:
    """a unit that belongs to an OrderedPolymer. Each unit has a direction, a location, and a link back to its parent"""

    def __init__(self, direction=None, position=None, parent=None):
        """the default is that the monomer is not part of a polymer"""

        self.parent = None
        self.direction = None
        self.position = None  # Prevents weird testing errors of not having attributes
        # by default, we assume that an orderedmonomer is not part of a polymer
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
        if parent is None or isinstance(parent, MonomerCollection):
            self._parent = parent
        else:
            raise ValueError(
                f"parent must be an MonomerCollection. Recieved {parent}")

    @property
    def direction(self):
        return self._direction

    @direction.setter
    def direction(self, direction):
        self._direction = direction

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        if self.parent is not None and position is None:
            raise ValueError(
                "{} is part of a polymer with no position!".format(self))
        else:
            self._position = position

    def find_polymer_component(self):
        from .species import ComplexSpecies
        outpolymer = None
        if (isinstance(self, ComplexSpecies)):
            for specie in self.species:
                if (specie.is_polymer_component):
                    if (outpolymer is not None):
                        raise ValueError(
                            "multiple species are part of the polymer in the same place!!")
                    else:
                        outpolymer = specie
        if (self.is_polymer_component):
            if (outpolymer is not None):
                raise ValueError(
                    "multiple species are part of the polymer in the same place!!")
            else:
                outpolymer = self
        return outpolymer

    def monomer_insert(self, parent: OrderedPolymer, position: int, direction=None):
        if (position is None):
            raise ValueError(
                "{} has no position to be inserted at!".format(self))
        if (direction is None):
            if (self.direction is not None):
                direction = self.direction
        if (parent is None):
            raise ValueError(
                "{} is trying to be inserted into nothing!".format(self))
        if (self.is_polymer_component is False):
            if (self.find_polymer_component() is None):
                self.is_polymer_component = True
        self.parent = parent
        self.position = position
        self.direction = direction

    def set_dir(self, direction):
        self.direction = direction
        return (self)

    def remove(self):
        self.parent = None
        self.position = None
        self.direction = None

        return (self)

    def get_orphan(self):
        """returns a copy of this monomer, except with no parent. But it still has a position and direction"""
        copied_monomer = copy.copy(self)
        copied_monomer.parent = None
        return (copied_monomer)

    def get_removed(self):
        copied_part = copy.copy(self)
        copied_part.parent = None
        copied_part.direction = None
        copied_part.position = None
        if (hasattr(copied_part, "_attributes")):
            copied_part.remove_attribute("forward")
            copied_part.remove_attribute("reverse")
        return copied_part

    def __repr__(self):
        txt = "OrderedMonomer(direction="+str(self.direction)+",position=" +\
            str(self.position)+")"
        return txt

    def __eq__(self, other):
        if (isinstance(other, OrderedMonomer)):
            if (self.direction == other.direction and self.position == other.position and self.parent == other.parent):
                return True
        return False

    def __hash__(self):
        hval = 0
        hval += self.subhash()

        if (self.parent is not None):
            hval += hash(self.parent)

        return hval

    def subhash(self):
        hval = 0
        hval += hash(self.position)
        hval += hash(self.direction)
        if (hasattr(self, "name")):
            hval += hash(self.name)
        return hval
