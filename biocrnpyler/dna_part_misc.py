
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from .dna_part import DNA_part
from .mechanisms_binding import One_Step_Cooperative_Binding
from .species import Species, ComplexSpecies, Complex
from .mechanisms_integrase import BasicIntegration
from .component import Component

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
        self.integrase = Component.set_species(integrase,material_type='protein')
        #self.integrase = integrase
        self.dinucleotide = dinucleotide
        self.site_type = site_type
        self.other_dna = None
        self.linked_sites = {}
        self.complexed_version = None
        DNABindingSite.__init__(self,name,self.integrase,no_stop_codons=no_stop_codons,**keywords)
        self.add_mechanism(BasicIntegration(self.integrase.name))
    def __repr__(self):
        myname = self.name
        if(self.site_type in integrase_sites):
            myname += "_" + self.integrase.name
            if(self.dinucleotide != 1):
                myname += "_"+str(self.dinucleotide) 
        else:
            warn("warning! site {} has site_type {} which is not recognized".format(self.name,self.site_type))
        return myname
    def get_complexed_species(self,dna):
        recomp = Complex([dna,self.integrase,self.integrase],material_type=None)
        return recomp
    def update_reactions(self):
        #TODO implement "integrase" mechanism
        reactions = DNABindingSite.update_reactions(self)
        complex_parent = self.get_complexed_species(self.dna_to_bind).parent
        int_mech = self.get_mechanism('integration')
        #this next part generates integrase reactions
        for site in self.linked_sites:
            integrated_dnas = []
            #each integrase reaction will have an entry here
            
            #now we must know if the site pair got processed. If it did,
            #then generate the reaction. If it didn't then update the pair

            #however, if we are dealing with an intramolecular reaction, then
            #only return the reaction if the pair HASN'T been processed
            dict_inputs = []
            for transformation in self.linked_sites[site][0]:
                for pspec in transformation.parentsdict:
                    dict_inputs += [transformation.parentsdict[pspec]]

            if(len(set(dict_inputs))== 1):
                #this is an intramolecular reaction so we don't care about the other site's DNAs
                integrated_dnas = []
                if(isinstance(complex_parent[self.position],ComplexSpecies) and \
                                isinstance(complex_parent[site.position],ComplexSpecies) and \
                                self.linked_sites[site][1]==[]):
                    #make sure that both this site and the other site are bound
                    for integrase_function in self.linked_sites[site][0]:
                        integr = integrase_function.create_polymer(complex_parent)
                        integrated_dnas += [integr]
                
                reactions += int_mech.update_reactions([complex_parent],integrated_dnas,component=self,part_id = self.integrase.name)
            else:
                
                #go through all possible "other" dnas then calculate for each RNA
                for other_dna in self.linked_sites[site][1]:
                    integrated_dnas = []
                    for integrase_function in self.linked_sites[site][0]:
                        #for every result that we could get, generate it
                        integrated_dnas += [integrase_function.create_polymer(complex_parent,other_dna)]
                    reactions += int_mech.update_reactions([complex_parent,other_dna],integrated_dnas,component=self,part_id = self.integrase.name)
            #next part updates the linked site
            if(complex_parent is not None and complex_parent not in site.linked_sites[self][1]):
                site.linked_sites[self][1] += [complex_parent]
        return reactions
            
