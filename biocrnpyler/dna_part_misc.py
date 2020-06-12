from .chemical_reaction_network import Species
from .mechanisms_binding import One_Step_Cooperative_Binding
from .dna_part import DNA_part
from warnings import warn

integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]

class AttachmentSite(DNA_part):
    def __init__(self,name,site_type,integrase="int1",integrase_cooperativity=2,dinucleotide=1,no_stop_codons=["forward","reverse"],**keywords):
        """an integrase attachment site binds to integrase"""
        self.mechanisms = {"binding":One_Step_Cooperative_Binding()}
        DNA_part.__init__(self,name,no_stop_codons=no_stop_codons,mechanisms = self.mechanisms,**keywords)
        self.name = name
        self.assembly = None
        self.integrase_cooperativity=integrase_cooperativity
        if(isinstance(integrase,Species)):
            self.integrase_species = integrase
            self.integrase = integrase.name
        elif(isinstance(integrase,str)):
            self.integrase_species = Species(integrase,material_type="protein")
            self.integrase = integrase
        #self.integrase = integrase
        self.dinucleotide = dinucleotide
        self.site_type = site_type
    def __repr__(self):
        myname = self.name
        if(self.site_type in integrase_sites):
            myname += "_" + self.integrase
            if(self.dinucleotide != 1):
                myname += "_"+str(self.dinucleotide) 
        else:
            warn("warning! site {} has site_type {} which is not recognized".format(self.name,self.site_type))
        return myname
    def update_species(self):
        if(self.assembly is None):
            return [self.integrase_species]
        else:
            mech_b = self.mechanisms["binding"]
            spec = mech_b.update_species(self.integrase_species,self.assembly.dna,component=self,part_id = self.site_type+"-"+self.integrase)
            spec += [self.integrase_species]
        return spec
    def update_reactions(self):
        rxns = []
        if(self.assembly is None):
            return []
        else:
            mech_b = self.mechanisms["binding"]
            rxns = mech_b.update_reactions(self.integrase_species,self.assembly.dna,component=self,part_id = self.site_type+"-"+self.integrase)
        return rxns