"""
The classes OrderedPolymer and OrderedMonomer are datastructures used to represent Polymers and their associatd components.
These classes are used by Chemical Reaction Network Species as well as certain Components such as DNA_construct.
"""
class OrderedPolymer:

    """a polymer made up of OrderedMonomers that has a specific order"""
    def __init__(self,parts):
        """parts can be a list of lists containing 
        [[OrderedMonomer,direction],[OrderedMonomer,direction],...]
        alternatively, you can have a regular list, and the direcitons
        will end up being None"""
        polymer = []
        assert(type(parts)==list or type(parts)==tuple), "OrderedPolymer must be instantiated with a list"
        for item in parts:
            if(isinstance(item,list)):
                part = item[0]
                if(len(item)>1):
                    partdir = item[1]
                else:
                    partdir = None
            elif(isinstance(item,OrderedMonomer)):
                part = item
                partdir = item.direction
            else:
                raise ValueError("{} is not an OrderedMonomer or a list of the form [OrderedMonomer,direction]".format(str(part)))
            
            polymer += [part]
            position = len(polymer)-1
            direction = partdir
            part.monomer_insert(self,position,direction)

        self._polymer = tuple(polymer)

    def __hash__(self):
        return hash(self._polymer)

    def insert(self,position,part,direction=None):
        part.monomer_insert(self,position,direction)
        for subsequent_part in self._polymer[position:]:
            subsequent_part.position += 1
        self._polymer = self._polymer[:position]+(part,)+self._polymer[position:]
        if(hasattr(self,"name") and hasattr(self,"make_name")):
            self.name = self.make_name()

    def replace(self,position,part,direction=None):
        if(direction is None):
            direction = part.direction
        self._polymer[position].remove()
        part.monomer_insert(self,position,direction)
        self._polymer = self._polymer[:position]+(part,)+self._polymer[position+1:]
        if(hasattr(self,"name") and hasattr(self,"make_name")):
            self.name = self.make_name()

    def append(self,part,direction=None):
        if(direction is None):
            if(hasattr(part,"direction")):
                direction = part.direction
            else:
                direction = None
        pos = len(self._polymer)
        self.insert(pos,part,direction)

    def __repr__(self):
        outstr = "polymer("
        for part in self._polymer:
            outstr += str(part)+", direction = "+str(part.direction)+","
        if(outstr[:-1]==","):
            outstr = outstr[:-1]
        outstr += ")"
        return outstr

    def direction_invert(self,dirname):
        if(dirname == "forward"):
            return "reverse"
        elif(dirname == "reverse"):
            return "forward"
        elif(dirname == 0):
            return 1
        elif(dirname == 1):
            return 0
        elif(dirname is None):
            return None
        else:
            warn("didn't know how to invert {}".format(str(dirname)))
            return dirname

    def __len__(self):
        return len(self._polymer)

    def __getitem__(self,ii):
        return self._polymer[ii]

    def __setitem__(self,ii,val):
        self.replace(ii,val,val.direction)

    def __eq__(self,other):
        if(isinstance(other,OrderedPolymer)):
            for item1,item2 in zip(self._polymer,other._polymer):
                if(item1.direction==item2.direction and item1.position==item2.position and type(item1)==type(item2)):
                    pass
                else:
                    return False
            if(len(self._polymer)==len(other._polymer)):
                return True
        return False

    def __contains__(self,item):
        if(item in self._polymer):
            return True
        else:
            return False

    def delpart(self,position):
        part = self._polymer[position]
        part.remove()
        for subsequent_part in self._polymer[position+1:]:
            subsequent_part.position -= 1
        self._polymer = self._polymer[:position] + self._polymer[position+1:]
        if(hasattr(self,"name") and hasattr(self,"make_name")):
            self.name = self.make_name()

    def reverse(self):
        self._polymer = self._polymer[::-1]
        for ind,part in enumerate(self._polymer):
            part.position = ind
            part.direction = self.direction_invert(part.direction)


class OrderedMonomer:
    """a unit that belongs to an OrderedPolymer. Each unit has a direction, a location, and a link back to its parent"""
    def __init__(self,direction=None,position=None,parent=None):
        """the default is that the monomer is not part of a polymer"""
        if(parent is not None):
            assert(isinstance(parent,OrderedPolymer)),"parent of OrderedMonomer must be an OrderedPolymer"
        self.direction = direction
        self.parent = parent
        self.position = position
        if self.parent is None:
            if self.position is not None:
                raise ValueError("{} cannot have a position if it has no parent".format(self))
        else:
            if(self.position is None):
                raise ValueError("{} is part of a polymer with no position!".format(self))
                
    def monomer_insert(self,parent:OrderedPolymer,position:int,direction=None):
        if(position is None):
            raise ValueError("{} has no position to be inserted at!".format(self))
        if(direction is None):
            if(self.direction is not None):
                direction  = self.direction
        if(parent is None):
            raise ValueError("{} is trying to be inserted into nothing!".format(self))
        self.parent = parent
        self.position = position
        self.direction = direction

    def set_dir(self,direction):
        self.direction = direction
        return(self)

    def remove(self):
        self.parent = None
        self.direction = None
        self.position = None

    def __repr__(self):
        txt = "OrderedMonomer(direction="+str(self.direction)+",position="+\
                                str(self.position)+")"
        return txt
    def __eq__(self,other):
        if(isinstance(other,OrderedMonomer)):
            if(self.direction == other.direction and self.position == other.position and self.parent == other.parent):
                return True
        return False