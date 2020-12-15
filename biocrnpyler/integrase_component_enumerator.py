# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from .component_enumerator import ComponentEnumerator
from .dna_part import DNA_part
from .dna_construct import DNA_construct
from .dna_part_misc import AttachmentSite
from .species import (OrderedComplexSpecies, Species, ComplexSpecies)
from .components_basic import Protein
import itertools as it

import copy
def combine_dictionaries(dict1,dict2):
    """append lists that share the same key, and add new keys"""
    outdict = copy.deepcopy(dict1)
    for key in dict2:
        if key in outdict:
            assert(isinstance(dict2[key],list))
            assert(isinstance(outdict[key],list))
            outdict[key] += dict2[key]
        else:
            outdict[key] = dict2[key]
    return outdict

def all_paths(prototype_list):
    """recursively enumerate all paths through a list"""
    if(len(prototype_list)==1):
        #base case
        return [[a] for a in prototype_list[0]]
    elif(len(prototype_list)>1):
        rest = all_paths(prototype_list[1:])
        #this is the recursive part
        retlist = []
        for a in prototype_list[0]:
            #for every member of the first list
            for rlist in rest:
                #put in every possible rest of the list
                retlist+=[[a]+rlist]
    return retlist

def list_integrase(construct):
    """lists all the parts that can be acted on by integrases"""
    int_dict = {}
    for part in construct.parts_list:
        if(isinstance(part,AttachmentSite) and part.integrase is not None):
            if(part.integrase in int_dict):
                int_dict.update({part.integrase:int_dict[part.integrase]+[part]})
            else:
                int_dict[part.integrase]=[part]
    return int_dict

class Integrase(ComponentEnumerator):
    def __init__(self,name:str,enumerator_type="integrase"):
        ComponentEnumerator.__init__(self,name=name,enumerator_type=enumerator_type)
    def explore_integrases(self,construct_list,int_mechanisms):
        """this explores all the possible integrase-motivated DNA configurations. If some
        integrases aren't present, then define intnames to be a list of names of the
        integrases which are present.
        An integrase can act in different ways. 
        * serine integrases recombine B and P sites that turn into L and R sites, 
                    and only sites with the same dinucleotide can be recombined.
        * serine integrases with directionality factors recombine L and R sites
                    with the same dinucleotide
        * Invertases only do flipping reactions
        * resolvases only do deletion reactions
        * FLP or CRE react with homotypic sites, so site1+site1 = site1+site1. But
                    there are still different types of sites which are orthogonal. For
                    example, a CRE type 1 or a CRE type 2 site. The sites can also be palindromic,
                    which means that they can react in either direction.
        """
        int_dict = {}
        for construct in construct_list:
            #list each integrase that exists and which sites they react with
            con_dict = list_integrase(construct)
            
            int_dict = combine_dictionaries(int_dict,con_dict)
        constructlist = []
        rxn_list = []
        integ_component_list = [] #not sure what this is
        #print(self.integraselist)
        for integrase in int_dict:
            if(integrase in int_mechanisms):
                int_mech = int_mechanisms[integrase]
                #now, going through each one, generate the reactions and species that arise
                attsites = int_dict[integrase]
                #but now we need to know what kind of integrase reactions are possible
                attcombos = [a for a in it.combinations(attsites,2)]
                #print(attcombos)
                for combo in attcombos:
                    #first question: is this combo legal?
                    if(tuple([a.site_type for a in combo]) in int_mech.reactions):
                        #this means the reaction can exist
                        results,int_rxn = int_mech.integrate(combo[0],combo[1])
                        #self.add_integration(combo[0],combo[1],s1,s2)
                        #integ_component_list = 
                        #second question: what are the resulting dna_constructs and reaction-primed parts
                        constructlist += copy.deepcopy(results)
                        rxn_list += [int_rxn]
        global_int = GlobalIntegraseMechanism(rxn_list)
        return constructlist,global_int




class DNA_function:
    def __init__(self,outputlists,productsites,intbound = None,integrase_species=None):
        """this creates one sequence out of one or two input sequences
        reactants: Species of the reactants
        outputlists: dictionary of indexlists for the products. key is the output Species, value is the [[index,direction, identity],[index1,direction1,identity],...]
        productsites: {output_name1:[site1,site2,site3,...],output_name2:[site1,site2,...]} each site is a list which has [species,position]
        intbound: 
        """
        self.outputlists = outputlists
        self.integrase_species = integrase_species
        self.intbound = intbound
        self.update_reactants_products()
        self.productsites = productsites
    def replace_dna_name(self,namedict):
        """replace internal representation of dna names with whatever the user wants"""
        new_outputlists = {}
        new_productsites = {}
        for outputname in self.outputlists:
            if(outputname in namedict):
                #this means we are replacing the name of an output.
                #it occurs in two places: outputlists and also productsites
                new_outputlists[namedict[outputname]]=self.outputlists[outputname]
                new_productsites[namedict[outputname]]=self.productsites[outputname]
            else:
                #if we aren't replacing it, still remember that it exists
                new_outputlists[outputname] = self.outputlists[outputname]
                new_productsites[outputname]= self.productsites[outputname]
        for new_name in new_outputlists:
            #maybe we need to replace an input for some reason?
            dna_prototype = []
            for element in new_outputlists[new_name]:
                if(element[2] in namedict):
                    dna_prototype +=[[element[0],element[1],namedict[element[2]]]]
                else:
                    dna_prototype += [element]
            new_outputlists[new_name]=dna_prototype
        self.outputlists = new_outputlists
        self.productsites = new_productsites
        self.update_reactants_products()
    def update_reactants_products(self):
        """update reactants and products to coincide with outputlists"""
        identities = []
        for outputSpecies in self.outputlists:
            identities+=[a[2] for a in self.outputlists[outputSpecies]]
        self.number_of_inputs = len(list(set(identities)))
        self.reactants =  list(set(identities))
        self.products = [a for a in self.outputlists.keys()]
    def transform_ordered_complex(self,input_complexes):
        """take an input complex and convert it into a product
        inputComplexes is a list of OrderedComplexes"""
        assert(len(input_complexes)==self.number_of_inputs,"wrong number of complexes provided")
        #make sure we have the right number of inputs
        #TODO this does not work with two copies of the same thing that react together, does it!?
        for a in input_complexes:
            assert(isinstance(a,OrderedComplexSpecies), "{} is not an OrderedComplexSpecies".format(a))
        new_input_complexes = copy.deepcopy(input_complexes)
        input_complex_dict = {a.species[-1].name:a.species[:-1] for a in new_input_complexes}
        out_dict = self.transform_list(input_complex_dict)
        outproducts = []
        for out_spec_name in out_dict:
            out_spec = Species(out_spec_name,material_type="dna")
            #print("out_spec is "+str(type(out_spec)))
            unordered = [a[0] for a in out_dict[out_spec_name]]
            modified_parts_list = []
            for spec in unordered:
                if(isinstance(spec,DNA_part)):
                    # you can't have DNA_parts in an OrderedComplexSpecies,
                    # but in this case the only thing that can happen 
                    # is to make a species with integrase bound.
                    #TODO put in a make_complex() here?
                    partspec = Species(spec.name,material_type="dna")
                    coop = spec.integrase_cooperativity
                    part_complex = ComplexSpecies([partspec]+[spec.integrase_species]*coop)
                    modified_parts_list += [part_complex]
                elif(isinstance(spec,str) or isinstance(spec,Species)):
                    modified_parts_list += [spec]
            outproducts +=[OrderedComplexSpecies(modified_parts_list+[out_spec],\
                                material_type="ordered_complex")] #TODO just assume we want the material_type of the first complex..?
        return outproducts
    def transform_DNA_parts_list(self,input_complexes):
        """take an input complex and convert it into a product
        input_complexes is a dictionary of lists of DNA_part objects"""
        
        assert(len(input_complexes)==self.number_of_inputs,"wrong number of complexes provided")
        #make sure we have the right number of inputs
        #TODO this does not work with two copies of the same thing that react together, does it!?
        for input_sublist in input_complexes:
            #everything should be a DNA part
            for a in input_complexes[input_sublist]:
                assert(isinstance(a,DNA_part),"{} in {} is not a DNA_part".format(a,input_complexes[input_sublist]))
            #assert(sum([isinstance(a,DNA_part) for a in input_complexes[input_sublist]])==len(input_complexes[input_sublist]), "{} is not made of DNA_parts")
        new_input_complexes = copy.deepcopy(input_complexes)
        out_dict = self.transform_list(new_input_complexes)
        outproducts = {}
        for out_spec in out_dict:
            speclist = []
            for a in out_dict[out_spec]:
                if(a[1]=="reverse"):
                    a[0].reverse()
                a[1] = a[0].direction
                speclist += [a]
            outproducts[out_spec] =speclist
        return outproducts
    def is_valid_reactant(self,ordered_species):
        if(self.intbound is None):
            raise ValueError("intbound is not defined")
        if(isinstance(ordered_species,OrderedComplexSpecies)):
            pass #in this case keep going
        else:
            #only OrderedComplexSpecies can be valid reactants anyway
            return False
        dna_name = ordered_species.species[-1].name
        if(dna_name in self.reactants):
            #is the name of this dna something we recognize?
            bound_amt = len(self.intbound[dna_name])
            #this is an accumulator for counting how many valid integrase bound locations we have
            for bound_loc in self.intbound[dna_name]:
                #go through each location
                bound_spec = self.intbound[dna_name][bound_loc]
                #this is the species that should be bound at this location
                if(isinstance(ordered_species.species[bound_loc],ComplexSpecies) and bound_spec in ordered_species.species[bound_loc]):
                    #if integrase is bound, count success!
                    #TODO I don't care about how many copies of integrase are bound but that might be important
                    bound_amt-=1 #count it
            if(bound_amt < 0):
                #this should not happen?
                raise ValueError("there are too many integrases bound to {}! That seems wrong".format(ordered_species))
            elif(bound_amt == 0):
                #there are exactly the right amount of integrases bound
                return True
            else:
                return False

    def __repr__(self):
        outstr = "DNA_func "
        dup_reactants = copy.deepcopy(self.reactants)
        dup_prod = copy.deepcopy(self.products)
        #print(self.intbound)
        #print(self.productsites)
        for input_spec in self.reactants:
            intbounds = ', '.join([str(a) for a in self.intbound[input_spec]])
            outstr+=str(input_spec)+" + "
        outstr = outstr[:-3]+" "+intbounds+" ===> "
        for output_spec in self.products:
            #print(self.productsites)
            intbounds = ', '.join([str(a[1]) for a in self.productsites[output_spec]])
            outstr+=str(output_spec)+" "+intbounds+" + "
        outstr = outstr[:-3]
        return outstr
    def transform_list(self,input_lists):
        """input_lists is a dictionary of lists. {dna1:[part1,part2,..],dna2:[part1,part2,...]]}"""
        #in this case we have a list of parts as input. the keys are DNA species that match the ones
        #we gave when we created this instance

        #assert(sum([isinstance(a,Species) for a in input_lists])==len([a for a in input_lists]))
        outlst = {a:[] for a in self.outputlists}
        #outproducts_dict = {}
        #this is what the outputs look like
        for out_spec in self.outputlists:
            output_list = self.outputlists[out_spec]
            #print(list(output_list))
            for index,direction,identity in output_list:
                #print([input_lists[identity][index],direction])
                #print(outlst)
                outlst[out_spec] += [[input_lists[identity][index],direction]]
            for site in self.productsites[out_spec]:
                #now we swap out the sites on the output for the ones that should be there
                outlst[out_spec][site[1]]=[site[0],"forward"]
            #outproducts_dict[out_spec] = outlst
        return outlst



class IntegraseMechanism:
    def __init__(self,name=None,reactions={("attB","attP"):"attL",("attP","attB"):"attR"}):
        """The integrase mechanism is a mechanism at the level of DNA. It creates DNA species which
        the integrase manipulations would lead to. This mechanism does not create any reaction rates.
        We need to figure out how integrase binding will work before being able to create
        reactions and their corresponding rates
        This is intended to be a GLOBAL MECHANISM"""
        if(name is None):
            self.name = "int1"
        else:
            self.name = name
        self.integrase_species = Protein(self.name).get_species()
        self.reactions = reactions
        self.attsites = []
        for reactants in reactions:
            self.attsites+=list(reactants)+list(reactions[reactants])
        self.integrations_to_do = [] #these are the reactions that will be performed at compile time
    def binds_to(self):
        return self.attsites
    def generate_products(self,site1,site2):
        """generates DNA_part objects corresponding to the products of recombination"""
        #the sites should have the same integrase and dinucleotide, otherwise it won't work
        assert(site1.integrase == site2.integrase)
        assert(site1.dinucleotide == site2.dinucleotide)
        integrase = site1.integrase
        dinucleotide = site1.dinucleotide
        rxn1 = (site1.site_type,site2.site_type)
        rxn2 = (site2.site_type,site1.site_type)
        try:
            prod1 = self.reactions[rxn1]
        except KeyError:
            #this means the reaction is not possible!!
            raise KeyError("{} not found to react with {} in {}".format(site1,site2,self.reactions))
        try:
            prod2 = self.reactions[rxn2]
        except KeyError:
            raise KeyError("{} not found to react with {} in {}".format(site2,site1,self.reactions))
        
        part_prod1 = AttachmentSite(prod1,prod1,dinucleotide=dinucleotide,
                                    integrase=integrase,direction=site1.direction,
                                            color=site1.color,color2=site2.color)
        part_prod2 = AttachmentSite(prod2,prod2,dinucleotide=dinucleotide,
                                    integrase=integrase,direction=site2.direction,
                                            color=site2.color,color2=site1.color)
        if(site1.direction=="forward"):
            return (part_prod1,part_prod2)
        else:
            part_prod2.direction = site1.direction
            part_prod1.direction = site2.direction
            return(part_prod2,part_prod1)
    
    def integrate(self,site1,site2):
        """perform an integration reaction between the chosen sites and make new DNA_constructs"""
        #if one of the sites is not part of a construct then raise an error!
        if(not(isinstance(site1.assembly,DNA_construct))):
            raise ValueError("{} not part of a construct".format(site1))
        elif(not(isinstance(site2.assembly,DNA_construct))):
            raise ValueError("{} not part of a construct".format(site2))
        site1_initial_facing = copy.deepcopy(site1.direction) #keep track of what we reverse so we can flip it back!
        site2_initial_facing = copy.deepcopy(site2.direction)
        
        cutpos1 = site1.pos
        cutpos2 = site2.pos
        #below are the references to the sites in the products
        finalprod1 = None
        finalprod2 = None
        if(site1.assembly==site2.assembly):
            #these sites are part of the same piece of DNA, so they are going to do an intramolecular reaction
            if(site1.pos > site2.pos):
                #we reverse which position is where
                cutpos2 = site1.pos
                cutpos1 = site2.pos
                prod1,prod2 = self.generate_products(site2,site1)
                #integrase sites are converted into different sites according to this function
            else:
                prod1,prod2 = self.generate_products(site1,site2)
                #integrase sites are converted into different sites according to this function
            
            dna = copy.deepcopy(site1.assembly)
            dna_list = list(range(len(dna.parts_list)))
            circularity = dna.circular
            
            if(site1.direction == site2.direction):
                #if the sites point in the same direction, then we are doing a deletion reaction
                cutdna_list_parts = dna_list[:cutpos1]+[dna_list[cutpos1]]+dna_list[cutpos2+1:] #delete
                newdna_list_parts = [dna_list[cutpos1]]+dna_list[1+cutpos1:cutpos2]
                cutdna_list = [[a,"forward",dna.name] for a in cutdna_list_parts]
                newdna_list = [[a,"forward",dna.name] for a in newdna_list_parts]
                integ_func = DNA_function({"cutdna":cutdna_list,"newdna":newdna_list},\
                                                            {"cutdna":[[prod1,cutpos1]],"newdna":[[prod2,0]]})
                #we want to feed out this integ_func
                outdna = integ_func.transform_DNA_parts_list({dna.name:dna.parts_list})
                cutdna_transformed = outdna["cutdna"]
                newdna_transformed = outdna["newdna"]

                cutdna = DNA_construct(cutdna_transformed,circular=circularity) 
                #the part we cut out of stays the same topology
                finalprod1 = cutdna.parts_list[cutpos1]
                if(site1_initial_facing != prod1.direction):
                    cutdna.reverse()
                
                newdna2 = DNA_construct(newdna_transformed,circular=True) 
                integ_func.replace_dna_name({"cutdna":cutdna.name,"newdna":newdna2.name})
                #update integ_func with the actual names
                
                #the part that gets cut out is always a circle
                finalprod2 = newdna2.parts_list[0]
                if(site2_initial_facing != prod2.direction):
                    newdna2.reverse()
                newdna = [newdna2,cutdna]
            else:
                #this means we are dealing with an inversion
                inv_segment = [[a,"reverse"] for a in dna_list[cutpos1+1:cutpos2][::-1]]
                #the inverted segment is reverse relative to everything else
                
                invertdna_list = [[a,"forward"] for a in dna_list[:cutpos1+1]] + \
                                        inv_segment+ \
                                        [[a,"forward"] for a in dna_list[cutpos2:]]
                invertdna_list = [[a[0],a[1],dna.name] for a in invertdna_list]
                #above makes a list of indexes that tells the DNA_function below
                #how to assemble to new list of DNA_part objects
                
                integ_func = DNA_function({"invdna":invertdna_list},\
                                                {"invdna":[[prod1,cutpos1],[prod2,cutpos2]]})
                #once the new list has been assembled, we extract the correct DNA_part list
                outdna = integ_func.transform_DNA_parts_list({dna.name:dna.parts_list})
                inv_dna = outdna["invdna"]
                #now the DNA_part list is converted into a DNA_construct
                inverteddna = DNA_construct(inv_dna,circular=circularity)
                integ_func.replace_dna_name({"invdna":inverteddna.name})
                #update integ_func with the actual names
                #this just pulls out the attL and attR sites TODO unnecessary, now that we are using DNA_function?
                finalprod1 = inverteddna.parts_list[cutpos1]
                finalprod2 = inverteddna.parts_list[cutpos2]
                if(site1_initial_facing != prod1.direction):
                    #this is making sure the final DNA_construct is inverted into the correct orientation
                    inverteddna.reverse()
                newdna = [inverteddna]
        else:
            #TODO didn't modify intermolecular reactions to use the DNA_function yet
            #otherwise these sites are on different pieces of DNA, so they are going to combine
            #don't mangle the inputs
            dna1 = copy.deepcopy(site1.assembly)
            dna2 = copy.deepcopy(site2.assembly)
            #fix references to sites
            site1 = dna1.parts_list[site1.pos]
            site2 = dna2.parts_list[site2.pos]
            #make sure everyone is forwards
            if(site1.direction=="reverse"):
                dna1.reverse()
            if(site2.direction=="reverse"):
                dna2.reverse()
            
            
            dna1_list = [[a,dna1.name] for a in list(range(len(dna1.parts_list)))] #list of indexes, for the DNA_function
            dna2_list = [[a,dna2.name] for a in list(range(len(dna2.parts_list)))]
            circ1 = dna1.circular
            circ2 = dna2.circular
            prod1,prod2 = self.generate_products(site1,site2)

            #generate attL and attR sites
            prod1.direction = "forward"
            prod2.direction = "forward"
            #direction of everything should be forward
            #now, cut the DNAs
            #flip the DNAs around so that they are pointing forwards
            
            
            if(circ2==True):
                #in this case we are combining a circular plasmid with a circular or linear plasmid
                #either way the result is basically the same, except the result is either linear or circular
                new_construct_list_ind = dna1_list[:site1.pos+1]+dna2_list[site2.pos+1:]+dna2_list[:site2.pos]+dna1_list[site1.pos:]
                new_construct_list_dir = [[a[0],"forward",a[1]] for a in new_construct_list_ind]
                new_construct_list = {"newdna":new_construct_list_dir}
                #computing the integrase site positions
                prod1_pos = site1.pos
                prod2_pos = site1.pos+len(dna2.parts_list)
                #constructing DNA transformation function
                integ_func = DNA_function(new_construct_list,{"newdna":[[prod1,prod1_pos],[prod2,prod2_pos]]})
                #applying DNA transformation function
                new_parts_lists = integ_func.transform_DNA_parts_list({dna1.name:dna1.parts_list,dna2.name:dna2.parts_list})
                #make a DNA_construct out of the transformed list
                if(circ1==False):
                    #if only the first one was a circle then the result is linear
                    new_construct = DNA_construct(new_parts_lists["newdna"],circular=False)
                elif(circ1==True):
                    #if both were circles then make the result a circle
                    new_construct = DNA_construct(new_parts_lists["newdna"],circular=True)
                integ_func.replace_dna_name({"newdna":new_construct.name})
                #store results for returning later
                newdna = [new_construct]
                finalprod1 = new_construct.parts_list[prod1_pos]
                finalprod2 = new_construct.parts_list[prod2_pos]
            elif(circ2 ==False and circ1 == True):
                #if the sites are backwards just reverse everything
                return self.integrate(site2,site1)
            elif(circ1==False and circ1==circ2):
                #here we are recombining two linear plasmids, so two linear plasmids are produced
                newdna_1_list = dna1_list[:site1.pos]+dna2_list[site2.pos:]
                newdna_2_list = dna2_list[:site2.pos]+dna1_list[site1.pos:]
                newdna1_list_dir = [[a[0],"forward",a[1]] for a in newdna_1_list]
                newdna2_list_dir = [[a[0],"forward",a[1]] for a in newdna_2_list]
                prod1_pos = site1.pos
                prod2_pos = site2.pos
                integ_func = DNA_function({"newdna1":newdna1_list_dir,"newdna2":newdna2_list_dir},
                                            {"newdna1":[[prod1,prod1_pos]],"newdna2":[[prod2,prod2_pos]]})
                reactant_dicts = {dna1.name:dna1.parts_list,dna2.name:dna2.parts_list}
                new_parts = integ_func.transform_DNA_parts_list(reactant_dicts)
                newdna_1 = DNA_construct(new_parts["newdna1"],circular=False)
                newdna_2 = DNA_construct(new_parts["newdna2"],circular=False)
                #create the new constructs
                integ_func.replace_dna_name({"newdna1":newdna_1.name,"newdna2":newdna_2.name})
                #now that we know what the new dnas will be called, put their names in there!
                finalprod1 = newdna_1.parts_list[site1.pos]
                
                finalprod2 = newdna_2.parts_list[site2.pos]
                #flip DNA so it's oriented correctly!
                if(site1_initial_facing != prod1.direction):
                    newdna_1.reverse()
                if(site2_initial_facing != prod2.direction):
                    newdna_2.reverse()
                newdna = [newdna_1,newdna_2]
        intbound = {}
        #now we are remembering where the integrases should be bound for the reactants
        for site in [site1,site2]:
            dna_name = site.assembly.name
            boundpos = site.pos
            integrase = Protein(site.integrase).get_species()
            if(dna_name in intbound):
                intbound[dna_name][boundpos] = integrase
            else:
                intbound[dna_name] = {boundpos:integrase}
        integ_func.intbound = intbound
        integ_func.integrase_species = self.integrase_species
        return newdna,integ_func






#integrase_construct doesn't work because it has no way of knowing about 
#other dna_constructs
class GlobalIntegraseMechanism(GlobalMechanism):
    def __init__(self,int_rxns,name= "integrases",mechanism_type = "integration", 
                filter_dict = {},
                 default_on = True):
        self.int_rxns = int_rxns
        self.filter_dict = filter_dict
        self.default_on = default_on
        GlobalMechanism.__init__(self, name=name, 
                                mechanism_type=mechanism_type,
                                default_on = default_on,
                                filter_dict=filter_dict)
    def append_rxns(self,rxn_list):
        """adds integrase reactions to this object"""
        if(isinstance(rxn_list,list)):
            for a in rxn_list:
                if(not isinstance(a,DNA_function)):
                    raise ValueError(str(a) + " is not a DNA_function")
            self.int_rxns += rxn_list
        elif(isinstance(rxn_list,DNA_function)):
            self.int_rxns += [rxn_list]
        else:
            raise ValueError(str(rxn_list) + " is not a DNA_function or a list of such")
        
    
    def update_reactions(self,rxnset,rxn,parameters):
        """actually creates the reactions"""
        rxnlist = []
        #for rxnset,rxn in rxn_dict_list:
        rxn_prototype = [rxnset[a] for a in rxnset]
        #print("rxn prototype")
        #print()
        #print(rxn_prototype)
        #print()
        rxn_allcomb = all_paths(rxn_prototype)
        rxn_k = self.get_parameter(rxn.integrase_species,parameters,"kint")
        #get the forward integration rate
        try:
            #if the reverse rate is not specified, assume it's zero
            rxn_k_rev = self.get_parameter(rxn.integrase_species,parameters,"kintrev")
        except ValueError:
            rxn_k_rev = 0
        for combo in rxn_allcomb:
            #
            outputs = rxn.transform_ordered_complex(combo)

            myrxn = Reaction(inputs=combo,outputs=outputs,k=rxn_k,k_rev=rxn_k_rev,details="integrase")
            rxnlist+=[myrxn]
        return rxnlist


    def update_reactions_global(self,species_list,parameters):
        """this function runs once, and generates all the appropriate integrase
        reactions between properly bound DNA + integrase molecules. This reaction
        does NOT account for tetramer formation"""

        rxnlist = []
        #print()
        #print("integrase reactions!")
        #print()
        #print(self.int_rxns)
        #print()
        unique_species = list(set(species_list))
        for rxn in self.int_rxns:
            rxnset = {}
            #this is going to run once for each reaction
            for specie in unique_species:
                if(rxn.is_valid_reactant(specie)):
                    #TODO can we do this with dictionaries? Is that faster at all?
                    specname = specie.species[-1].name
                    #here what we are doing is finding all the complexes
                    #that involve the same specie and putting them in a dictionary
                    #filed under that specie.
                    if(specname in rxnset):
                        rxnset[specname] += [specie]
                    else:
                        rxnset[specname] = [specie]
            #there should be one "rxn" for each unique integration event
            rxnlist += self.update_reactions(rxnset,rxn,parameters)

            #rxn_dict_list += [[rxnset,rxn]]
        return rxnlist