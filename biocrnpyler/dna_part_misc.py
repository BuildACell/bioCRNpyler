
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from .dna_part import DNA_part
from .mechanisms_binding import One_Step_Cooperative_Binding
from .species import Species, ComplexSpecies

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
        self.dna_to_bind = None
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
    def update_component(self,internal_species=None,**keywords):
        """returns a copy of this component, except with the proper fields updated"""
        from .dna_construct import RNA_construct,DNA_construct
        if(isinstance(self.parent,DNA_construct)):
            out_component = copy.copy(self)
            out_component.dna_to_bind = internal_species
            return out_component
        else:
            return None
class AttachmentSite(DNABindingSite):
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
        self.other_dna = None
        self.linked_sites = {}
        self.complexed_version = None
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
    def update_species(self):
        print(f"update species of {self}")
        species = DNABindingSite.update_species(self)
        self.complexed_version = None
        for spec in species:
            if(isinstance(spec,ComplexSpecies) and self.complexed_version is None):
                self.complexed_version = spec
            elif(isinstance(spec,ComplexSpecies) and spec != self.complexed_version):
                raise ValueError(f"in {self}, {spec} is a complex but I already found a complex called {self.complexed_version}")
        #go through all the sites i am linked to
        complex_parent = self.complexed_version.parent
        for site in self.linked_sites:
            print(f"other dnas are {self.linked_sites[site][1]}")
            integrated_dnas = []
            for integrase_function in self.linked_sites[site][0]:
            #if they already got processed then they will have deposited their info here
                if(integrase_function.number_of_inputs == 1):
                    #this is an intramolecular reaction so we don't care about the other guy's DNAs
                    if(isinstance(complex_parent[self.position],ComplexSpecies) and \
                                isinstance(complex_parent[site.position],ComplexSpecies)):
                        print(f"recombining {complex_parent}")
                        integr = integrase_function.create_polymer(complex_parent)
                        print(f"obtained {integr}")
                        integrated_dnas += [integr]
                else:
                    #if we do care, then calculate for each RNA
                    for other_dna in self.linked_sites[site][1]:
                        #for every result that we could get, generate it
                        print(f"recombining {complex_parent} with {other_dna.parent}")
                        integrated_dnas += [integrase_function.create_polymer(complex_parent,other_dna.parent)]
            species += integrated_dnas
            #populate their sites
            #this should probably happen in update_reactions
            if(complex_parent not in site.linked_sites[self][1]):
                print(f"seeding {site} with {complex_parent}")
                site.linked_sites[self][1] += [complex_parent]
        return species