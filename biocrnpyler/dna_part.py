################################################################
#       DNA_part: a component-like intermediate class necessary for DNA_construct
#       Author: Andrey Shur
#       Latest update: 6/4/2020
#
#
################################################################


# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from .component import Component
from .polymer import OrderedMonomer, OrderedPolymer
from .species import Species
import copy


class DNA_part(Component, OrderedMonomer):
    def __init__(self,name, **keywords):
        """this represents a modular component sequence. These get compiled into working components"""

        if "initial_concentration" in keywords:
            raise AttributeError("DNA_part should not recieve initial_concentration keyword. Pass this into the DNAassembly or DNA_construct instead.")

        Component.__init__(self=self, name = name, **keywords)
        
        self.name = name
        assembly = None #this is already covered in self.assembly. Is this needed?
        direction = None #orientation
        pos = None #position in the dna_construct
        self.sequence = None #nucleotide sequence
        #most parts have stop codons
        #if you want your part to not have stop codons, put "forward" and/or "reverse"
        self.no_stop_codons = [] #some parts will set this, others won't. Default is that it has stop codons
        self.material_type = "part"
        for key, value in keywords.items():
            #goes through any extra parameters and sets them!
            if(key=="assembly"):
                assembly = value
            elif(key=="direction"):
                direction = value
            elif(key=="pos"):
                pos = value
            elif(key=="sequence"):
                self.sequence = value
            elif(key=="no_stop_codons"):
                if(value is not None):
                    self.no_stop_codons = value
            elif(key=="material_type"):
                self.material_type = value
        if(isinstance(assembly,OrderedPolymer)):
            OrderedMonomer.__init__(self,position=pos,parent=assembly,direction=direction)
        else:
            self.assembly = assembly
            OrderedMonomer.__init__(self,position=pos,direction=direction)

    @property
    def dna_species(self):
        return Species(self.name, material_type="part")

    def __repr__(self):
        myname = self.name
        if(self.position is not None):
            myname+="_"+str(self.position)
        if(self.direction =="reverse"):
            myname += "_r"
        return myname
    def __hash__(self):
        return OrderedMonomer.__hash__(self)+hash(self.name)
    def __eq__(self,other):
        if(type(other)==type(self)):
            if(self.name==other.name):
                if(self.assembly is not None and other.assembly is not None):
                    if(str(self.assembly)==str(other.assembly)):
                        return True
                elif((self.assembly is not None) or (other.assembly is not None)):
                    #if one has an assembly and the other doesn't, then they aren't the same!!
                    return False
                elif(((self.parent is None) and (other.parent is None)) or \
                    ((self.parent is not None) and (other.parent is not None) and (str(self.parent)==str(other.parent)))):
                    #this is for when we are using the OrderedMonomer for its intended function
                    if(self.direction==other.direction and self.position==other.position):
                        return True
        return False

    def clone(self,position,direction,parent_dna):
        """this defines where the part is in what piece of DNA"""
        #TODO add warning if DNA_part is not cloned
        self.insert(parent_dna,position,direction)
        return self

    def unclone(self):
        """removes the current part from anything"""
        self.remove()
        return self
        
    def reverse(self):
        OrderedMonomer.reverse(self)
        return self
