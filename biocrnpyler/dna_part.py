################################################################
#       DNA_part: a component-like intermediate class necessary for DNA_construct
#       Author: Andrey Shur
#       Latest update: 6/4/2020
#
#
################################################################

from .component import Component
from .chemical_reaction_network import OrderedMonomer, Species
from warnings import warn

class DNA_part(Component,OrderedMonomer):
    def __init__(self,name,mechanisms={},parameters={},**keywords):
        """this represents a modular component sequence. These get compiled into working components"""
        Component.__init__(self=self, name = name, mechanisms = mechanisms,
                           parameters = parameters, **keywords)
        
        self.name = name
        #self.assembly = property(self._get_assembly,self._set_assembly)
        assembly = None #this is already covered in self.assembly. Is this needed?
        direction = None #orientation
        self.color = None #this will be taken from higher up
        self.color2 = None #this will be taken from higher up and only defined for attL/R sites
        pos = None #position in the dna_construct
        self.sequence = None #nucleotide sequence
        self.no_stop_codons = [] #some parts will set this, others won't. Default is that it has stop codons
        
        for key, value in keywords.items():
            #goes through any extra parameters and sets them!
            if(key=="assembly"):
                assembly = value
            elif(key=="direction"):
                direction = value
            elif(key=="color"):
                self.color = value
            elif(key=="color2"):
                self.color2 = value
            elif(key=="pos"):
                pos = value
            elif(key=="sequence"):
                self.sequence = value
            elif(key=="no_stop_codons"):
                self.no_stop_codons = value
        OrderedMonomer.__init__(self,data=None,position=pos,parent=assembly,direction=direction)
    @property
    def dna_species(self):
        return Species(self.name, material_type="dna")
    def __repr__(self):
        myname = self.name
        if(self.position is not None):
            myname+="_"+str(self.position)
        if(self.direction =="reverse"):
            myname += "_r"
        return myname
    def clone(self,position,direction,parent_dna):
        """this defines where the part is in what piece of DNA"""
        #TODO add warning if DNA_part is not cloned
        self.insert(parent_dna,position,direction)
        #if(self.assembly is not None):
        #    warn(str(self) + " already belongs to "+str(self.assembly.name)+"! It will now be part of the new assembly")
        #    self.unclone()

        #self.pos = position
        #self.direction = direction
        #self.assembly = parent_dna
        return self
    def _set_assembly(self,assembly):
        """set the "assembly" variable"""
        self.parent = assembly
        if(assembly is None):
            self.dna_to_bind = None
        else:
            self.dna_to_bind = assembly.dna
    def _get_assembly(self):
        return self.parent
    def unclone(self):
        """removes the current part from anything"""
        self.remove()
        #if(self.assembly is not None):
        #    rightparts = self.assembly.parts_list[self.pos+1:]
        #    for part in rightparts:
        #        part.pos -= 1
        #    self.assembly.parts_list = self.assembly.parts_list[:self.pos]+rightparts
        #self.pos = None
        #self.direction = None
        #self.assembly = None
        return self
    def reverse(self):
        self.reverse()
        #if(self.direction=="forward"):
        #    self.direction = "reverse"
        #elif(self.direction=="reverse"):
        #    self.direction = "forward"
        #elif(self.direction==None):
        #    warn(str(self)+" has no direction. Perhaps it wasn't cloned?")
        #else:
        #    raise ValueError("direction is not forward or reverse! It's "+self.direction)
        return self