from .species import Species
from .mechanisms_binding import One_Step_Cooperative_Binding
from .dna_part import DNA_part
from warnings import warn
import copy

integrase_sites = ["attB","attP","attL","attR","FLP","CRE"]
class DNABindingSite(DNA_part):
    def __init__(self,name,binders,no_stop_codons=None,assembly = None,**keywords):
        """an integrase attachment site binds to integrase"""

        if(isinstance(binders,list)):
            self.binders = [self.set_species(a) for a in binders]
        else:
            self.binders = [binders]
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
            for binder in self.binders:
                rxns += mech_b.update_reactions(binder,self.dna_to_bind,component=self,part_id = self.name)
        return rxns
    def update_component(self,dna,rnas=None,proteins=None,mypos = None):
        """returns a copy of this component, except with the proper fields updated"""
        out_component = copy.deepcopy(self)
        if(mypos is not None):
            out_component.dna_to_bind = dna[mypos]
        else:
            out_component.dna_to_bind = dna
        if(dna.material_type == "rna"):
            #DNA binding sites only work with DNA
            return None
        elif(dna.material_type == "dna"):
            return out_component
class AttachmentSite(DNABindingSite):
    #TODO generalize to "DNA binding site"
    def __init__(self,name, site_type = "attB",integrase = "int1", dinucleotide = 1,no_stop_codons=None,**keywords):
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
