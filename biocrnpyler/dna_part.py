################################################################
#       DNA_part: a component-like intermediate class necessary for DNA_construct
#       Author: Andrey Shur
#       Latest update: 6/4/2020
#
#
################################################################

from .component import Component
from warnings import warn

class DNA_part(Component):
    def __init__(self,name,mechanisms={},parameters={},**keywords):
        """this represents a modular component sequence. These get compiled into working components"""
        Component.__init__(self=self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)
        self.name = name
        self.assembly = None #this is already covered in self.assembly. Is this needed?
        self.direction=None #orientation
        self.color = None #this will be taken from higher up
        self.color2 = None #this will be taken from higher up and only defined for attL/R sites
        self.pos = None #position in the dna_construct
        self.sequence = None #nucleotide sequence
        self.no_stop_codons = [] #some parts will set this, others won't. Default is that it has stop codons
        for key, value in keywords.items():
            #goes through any extra parameters and sets them!
            if(key=="assembly"):
                self.assembly = value
            elif(key=="direction"):
                self.direction = value
            elif(key=="color"):
                self.color = value
            elif(key=="color2"):
                self.color2 = value
            elif(key=="pos"):
                self.pos = value
            elif(key=="sequence"):
                self.sequence = value
            elif(key=="no_stop_codons"):
                self.no_stop_codons = value
    def __repr__(self):
        myname = self.name
        if(self.pos is not None):
            myname+="_"+str(self.pos)
        if(self.direction =="reverse"):
            myname += "_r"
        return myname
    def clone(self,position,direction,parent_dna):
        """this defines where the part is in what piece of DNA"""
        self.pos = position
        self.direction = direction
        self.assembly = parent_dna
        #TODO make this a copy function, so we don't have circular references
        return self
    def unclone(self):
        """removes the current part from anything"""
        self.pos = None
        self.direction = None
        self.assembly = None
        return self
    def reverse(self):
        if(self.direction=="forward"):
            self.direction = "reverse"
        elif(self.direction=="reverse"):
            self.direction = "forward"
        elif(self.direction==None):
            warn(self.name+" has no direction. Perhaps it wasn't cloned?")
        else:
            raise ValueError("direction is not forward or reverse! It's "+self.direction)
        return self