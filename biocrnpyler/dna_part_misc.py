
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import copy
from warnings import warn

from .dna_part import DNA_part
from .mechanisms_binding import One_Step_Cooperative_Binding
from .species import Species, ComplexSpecies, Complex
from .mechanisms_integrase import BasicIntegration
from .component import Component
from .components_basic import DNA

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
        if(isinstance(self.parent,DNA)):
            out_component = copy.copy(self)
            out_component.dna_to_bind = internal_species
            return out_component
        else:
            return None
class IntegraseSite(DNABindingSite):
    def __init__(self,name, site_type = "attB",integrase = "int1", dinucleotide = 1,no_stop_codons=None,**keywords):
        self.update_integrase(integrase)
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
    def update_integrase(self,int_name):
        self.integrase = Component.set_species(int_name,material_type='protein')
    def __hash__(self):
        sumhash = DNABindingSite.__hash__(self) + self.dna_to_bind.__hash__()
        return sumhash
    def get_complexed_species(self,dna):
        recomp = Complex([dna,self.integrase,self.integrase])
        return recomp
    def update_component(self,internal_species=None,**keywords):
        """returns a copy of this component, except with the proper fields updated"""
        newcomp = DNABindingSite.update_component(self,internal_species=internal_species)
        if(newcomp is None):
            return None
        elif("practice_run" in keywords and keywords["practice_run"]):
            for othersite in self.linked_sites:
                mystuff = copy.copy(othersite.linked_sites[self])
                othersite.linked_sites[newcomp] = mystuff
                del othersite.linked_sites[self]
            return newcomp
        else:
            for othersite in self.linked_sites:
                populate = True #by default, add the appropriate data members to the other site
                if(self.count_inputs(othersite)==1):
                    #the reaction with the site in question is intramolecular
                    #that means we should only populate the other site if our internal_species
                    #has the proper location bound by integrase
                    otherisbound = othersite.integrase in internal_species.parent[othersite.position]
                    if(not otherisbound):
                        populate = False
                if(populate):
                    mystuff = copy.copy(othersite.linked_sites[self])
                    del othersite.linked_sites[self]
                    othersite.linked_sites[newcomp]=mystuff

                    assert(othersite in newcomp.linked_sites)
            return newcomp
    def count_inputs(self,site):
        assert(site in self.linked_sites)
        dict_inputs = []
        for transformation in self.linked_sites[site][0]:
            for pspec in transformation.parentsdict:
                dict_inputs += [transformation.parentsdict[pspec]]
        return len(set(dict_inputs))
    def update_reactions(self):
        reactions = DNABindingSite.update_reactions(self)
        if(self.linked_sites == {}):
            return reactions
        complex_parent = self.get_complexed_species(self.dna_to_bind).parent
        int_mech = self.get_mechanism('integration')

        #this next part generates integrase reactions
        for site in self.linked_sites:
            if(site.dna_to_bind is None or site.dna_to_bind.parent is None):
                #skip sites which do not know who they bind to
                continue
            integrated_dnas = []
            #each integrase reaction will have an entry here
            
            #now we must know if the site pair got processed. If it did,
            #then generate the reaction. If it didn't then update the pair

            #however, if we are dealing with an intramolecular reaction, then
            #only return the reaction if the pair HASN'T been processed
            populate = True
            if(self.count_inputs(site)== 1):
                #this is an intramolecular reaction so we don't care about the other site's DNAs
                integrated_dnas = []
                if(site.integrase in complex_parent[site.position] and \
                                complex_parent in self.linked_sites[site][1]):
                    #make sure that both this site and the other site are bound
                    for integrase_function in self.linked_sites[site][0]:
                        integr = integrase_function.create_polymer([complex_parent])
                        integrated_dnas += [integr]
                
                    reactions += int_mech.update_reactions([complex_parent],integrated_dnas,component=self,part_id = self.integrase.name)
                    populate = False
            else:
                #this is for intermolecular reactions. now "other_dna" is possible
                #go through all possible "other" dnas then calculate for each RNA
                for other_dna in self.linked_sites[site][1]:
                    #I am not convinced this loop runs more than once EVER
                    integrated_dnas = []
                    for integrase_function in self.linked_sites[site][0]:
                        #for every result that we could get, generate it
                        #the next line generates the OrderedPolymerSpecies which results from recombination
                        integrated_dnas += [integrase_function.create_polymer([complex_parent,other_dna])]
                    reactions += int_mech.update_reactions([complex_parent,other_dna],integrated_dnas,component=self,part_id = self.integrase.name)
                    populate = False #since the other site already populated us, no reason to populate it
            #next part updates the linked site
            if(populate and (complex_parent is not None) and self in site.linked_sites and (complex_parent not in site.linked_sites[self][1])):
                site.linked_sites[self][1] += [complex_parent]
        return reactions
            
