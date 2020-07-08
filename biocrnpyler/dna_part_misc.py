from .chemical_reaction_network import Species, make_species
from .mechanisms_binding import One_Step_Cooperative_Binding
from .dna_part import DNA_part
from warnings import warn

integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]
class DNABindingSite(DNA_part):
    def __init__(self,name,binders,no_stop_codons=[],assembly = None,**keywords):
        """an integrase attachment site binds to integrase"""
        self.binders = make_species(binders)
        if(not isinstance(self.binders,list)):
            self.binders = [self.binders]
        self.mechanisms = {"binding":One_Step_Cooperative_Binding()}
        DNA_part.__init__(self,name,no_stop_codons=no_stop_codons,mechanisms = self.mechanisms,assembly = assembly,**keywords)
        self.name = name
        #self.assembly = None
        
    def __repr__(self):
        myname = self.name
        return myname
    def update_species(self):
        spec = []
        spec += self.binders
        if(self.dna_to_bind is not None):
            mech_b = self.mechanisms["binding"]
            for binder in self.binders:
                spec += mech_b.update_species(binder,self.dna_to_bind,component=self,part_id = self.name)
                #TODO: different proteins probably have different affinity to bind this sequence
        return spec
    def update_reactions(self):
        rxns = []
        if(self.dna_to_bind is not None):
            mech_b = self.mechanisms["binding"]
            rxns = mech_b.update_reactions(self.binders,self.dna_to_bind,component=self,part_id = self.name)
        return rxns
class AttachmentSite(DNABindingSite):
    #TODO generalize to "DNA binding site"
    def __init__(self,name, site_type = "attB",integrase = "int1", dinucleotide = 1,no_stop_codons=[],**keywords):
        if(isinstance(integrase,Species)):
            self.integrase_species = integrase
            self.integrase = integrase.name
        elif(isinstance(integrase,str)):
            self.integrase_species = Species(integrase,material_type="protein")
            self.integrase = integrase
        #self.integrase = integrase
        self.dinucleotide = dinucleotide
        self.site_type = site_type
        DNABindingSite.__init__(self,name,self.integrase_species,no_stop_codons=no_stop_codons,**keywords)
    def __repr__(self):
        myname = self.name
        if(self.site_type in integrase_sites):
            myname += "_" + self.integrase
            if(self.dinucleotide != 1):
                myname += "_"+str(self.dinucleotide) 
        else:
            warn("warning! site {} has site_type {} which is not recognized".format(self.name,self.site_type))
        return myname