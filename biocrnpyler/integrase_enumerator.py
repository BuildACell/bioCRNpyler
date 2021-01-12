#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.
from .species import Species
from .polymer import OrderedPolymer, OrderedMonomer, NamedPolymer
from .dna_construct import Construct, DNA_construct
from .dna_part_misc import AttachmentSite
from .component_enumerator import GlobalComponentEnumerator
import itertools as it
import copy

class Polymer_transformation:
    def __init__(self,partslist,circular=False,parentsdict = None):
        """A Polymer transformation is like a generic transformation of a polymer sequence.
        You specify a parts list that would make up the output polymer. This list can contain:
        parts from ordered polymers
        parts that aren't in any polymers (have parent = None)
        parts from ordered polymers are considered as "placeholders"
        parts with parent = None are inserted into the new polymer.

        Also you can specify if the output should be circular or not.
        """
        if(parentsdict is None):
            parentsdict = {} #this is a conversion of the parents of the inputted parts into "blank" dummy polymers
                            #each unique parent becomes a different dummy polymer
        inputcount = 1
        actual_partslist = []
        partdir = -1 #-1 = not set. Valid values are "forward" "reverse" and None
        #go through the parts
        inputname = None
        for part in partslist:
            if(type(part)==list or type(part) ==tuple):
                #if the part is a tuple that means it looks like [OrderedMonomer,"direction"]
                partdir = part[1]
                part_ref = part[0]
            else:
                part_ref = part
                partdir = part_ref.direction
            #if the parent is populated, it means this part should be a placeholder
            if(part_ref.parent is not None):
                #if we haven't tracked it already
                if(part_ref.parent not in parentsdict):
                    #create a 'blank' polymer based on the input, and give it a generic name
                    inputname = "input"+str(inputcount)
                    dummyPolymer = self.dummify(part_ref.parent,inputname)
                    #set the value in the dictionary
                    parentsdict[part_ref.parent]=dummyPolymer.name
                    #increment for the next time this happens
                    
                    inputcount+=1
                elif(inputname is None or inputname != parentsdict[part_ref.parent]):
                    inputname = parentsdict[part_ref.parent]
                    dummyPolymer = self.dummify(part_ref.parent,inputname)
                #this variable is the actual parts list that will be stored. It has mostly dummy parts
                actual_partslist += [[dummyPolymer[part_ref.position],partdir]]
            else:
                #if the part has no parent, copy it and put it in
                actual_partslist += [[copy.deepcopy(part_ref),partdir]]
        self.number_of_inputs = len(parentsdict)
        self.parentsdict = parentsdict
        self.partslist = actual_partslist
        self.circular = circular
    def inverted(self):
        """return an "inverted" version of yourself where input1 = input2
        WARNING: this will not work with more than 2 sites!"""
        new_parentsdict = {}
        if(len(self.parentsdict)==1):
            #in this case we have only "input1"
            return self
        else:
            for part in self.partslist:
                if(part[0].parent is not None and part[0].parent not in new_parentsdict):
                    #this next line outputs input2 when given input1
                    newname = ["input1","input2"][part[0].parent.name=="input1"]
                    new_parentsdict[part[0].parent]=newname
            return Polymer_transformation(self.partslist,self.circular,parentsdict=new_parentsdict)
    def create_polymer(self,*args,bound=True,**keywords):
        inputcount = 1
        for arg in args:
            inputname = "input"+str(inputcount)
            assert(inputname not in keywords)
            keywords[inputname]=arg
            inputcount += 1
        assert(sum(["input" in a for a in keywords])>=self.number_of_inputs)
        outlst = []
        for part_list in self.partslist:
            part = part_list[0]
            partdir = part_list[1]
            outpart = None
            
            if(part.parent is not None):
                outpart = keywords[part.parent.name][part.position] #grab the part from the proper input
            else:
                #parts that don't come from input1 or input2.
                #they need to be either DNA_part or species objects, depending on
                #what kind of object we are making in the end.
                if(isinstance(keywords["input1"],Construct)):
                    outpart = part
                else:
                    #this is a new part which wasn't part of a polymer (like an attL site)
                    if(bound):
                        outpart = part.get_complexed_species(part.dna_species)
                    else:
                        outpart = part.dna_species
            #assuming the stored parts have a valid direction
            if(isinstance(outpart,AttachmentSite)):
                outpart = copy.copy(outpart)
                outpart.linked_sites = {} #make sure that any integrase sites we copy this way have no
                                            #linked sites, as those would not be links created by the integrate() function
            outlst += [[outpart,partdir]]
        outpolymer = keywords["input1"].__class__(outlst,circular = self.circular)
        return outpolymer
    def dummify(self,in_polymer,name):
        out_list = []    
        for element in in_polymer:
            out_list += [OrderedMonomer(direction=element.direction)]
        return NamedPolymer(out_list,name)
    def __repr__(self):
        part_texts = []
        for plist in self.partslist:
            part = plist[0]
            part_dir = plist[1]
            if(part.parent is not None):
                part_texts += [[part.parent.name,part.position,part_dir]]
            else:
                part_texts += [[part.name,part_dir]]
        out_txt = "Polymer transformation = "
        for part_text in part_texts:
            out_txt += "("+str(part_text[0])
            if(len(part_text)==3):
                out_txt +="-"+str(part_text[1])
            if(part_text[-1]!="forward"):
                out_txt += "-r"
            out_txt += ")"
        return out_txt

class IntegraseMechanism:
    def __init__(self,name=None,reactions={("attB","attP"):"attL",("attP","attB"):"attR"}):
        """The integrase mechanism is a mechanism at the level of DNA. It creates DNA species which
        the integrase manipulations would lead to. This mechanism does not create any reaction rates.
        We need to figure out how integrase binding will work before being able to create
        reactions and their corresponding rates"""
        if(name is None):
            self.name = "int1"
        else:
            self.name = name
        self.integrase_species = Species(self.name,material_type="protein")
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
                                            color=site2.color,color2=site1.color)
        part_prod2 = AttachmentSite(prod2,prod2,dinucleotide=dinucleotide,
                                    integrase=integrase,direction=site2.direction,
                                            color=site1.color,color2=site2.color)
        if(site1.direction=="forward"):
            return (part_prod1,part_prod2)
        else:
            part_prod2.direction = site1.direction
            part_prod1.direction = site2.direction
            return(part_prod2,part_prod1)
    
    def integrate(self,site1,site2):
        """perform an integration reaction between the chosen sites and make new DNA_constructs"""
        #if one of the sites is not part of a construct then raise an error!
        if(not(isinstance(site1.parent,Construct))):
            raise ValueError("{} not part of a construct".format(site1))
        elif(not(isinstance(site2.parent,Construct))):
            raise ValueError("{} not part of a construct".format(site2))
        cutpos1 = site1.position
        cutpos2 = site2.position
        #below are the references to the sites in the products
        dna_inputs = []
        if(site1.parent==site2.parent):
            #these sites are part of the same piece of DNA, so they are going to do an intramolecular reaction
            if(site1.position > site2.position):
                #we reverse which position is where
                cutpos2 = site1.position
                cutpos1 = site2.position
                prod1,prod2 = self.generate_products(site2,site1)
                #integrase sites are converted into different sites according to this function
            else:
                prod1,prod2 = self.generate_products(site1,site2)
                #integrase sites are converted into different sites according to this function
            
            dna = site1.parent
            dna_inputs = [dna]
            #dna = copy.deepcopy(site1.assembly)
            #dna_list = list(range(len(dna.parts_list)))
            circularity = dna.circular
            
            if(site1.direction == site2.direction):
                #if the sites point in the same direction, then we are doing a deletion reaction
                #direction doesn't matter so we don't need to flip anything
                cutdna_list_parts = list(dna[:cutpos1])+[[prod1,site1.direction]]+list(dna[cutpos2+1:]) #delete
                newdna_list_parts = [[prod2,site2.direction]]+list(dna[1+cutpos1:cutpos2])

                integ_funcs = [Polymer_transformation(cutdna_list_parts,circular = circularity),\
                                        Polymer_transformation(newdna_list_parts,circular=True)]
            else:
                #this means we are dealing with an inversion
                inv_segment = []
                
                [[a,] for a in dna[cutpos1+1:cutpos2][::-1]]
                for a in dna[cutpos1+1:cutpos2][::-1]:
                    inv_segment += [[a,["forward","reverse"][a.direction=="forward"]]]
                #the inverted segment is reversed
                
                invertdna_list = list(dna[:cutpos1])+\
                                [[prod1,site1.direction]]+\
                                inv_segment+ \
                                [[prod2,site2.direction]]+\
                                list(dna[cutpos2+1:])

                integ_funcs = [Polymer_transformation(invertdna_list,circular=circularity)]
        else:
            #otherwise these sites are on different pieces of DNA, so they are going to combine
            dna1 = site1.parent
            dna2 = site2.parent
            dna_inputs = [dna1,dna2]
            pdict = {a[1]:"input"+str(a[0]+1) for a in enumerate(dna_inputs)}
            #make sure everyone is forwards
            sites = [site1,site2]
            site_halves = []
            for dna_num,site_num in zip(dna_inputs,sites):
                if(site_num.direction=="reverse"):
                    dnanum_beginning = [[a,["forward","reverse"][a.direction=="forward"]] for a in dna_num[site_num.position+1:][::-1]]
                    dnanum_end = [[a,["forward","reverse"][a.direction=="forward"]] for a in dna_num[:site_num.position][::-1]]
                else:
                    dnanum_beginning = dna_num[:site_num.position]
                    dnanum_end = dna_num[site_num.position+1:]
                site_halves += [[list(dnanum_beginning),list(dnanum_end)]]
            dna1_halves = site_halves[0]
            dna2_halves = site_halves[1]

            circ1 = dna1.circular
            circ2 = dna2.circular
            prod1,prod2 = self.generate_products(site1,site2)

            #direction of everything should be forward

            if(circ2==True):
                #in this case we are combining a circular plasmid with a circular or linear plasmid
                #either way the result is basically the same, except the result is either linear or circular
                #result is ONE PIECE OF DNA
                result = dna1_halves[0]+[[prod1,"forward"]]+dna2_halves[1]+dna2_halves[0]+[[prod2,"forward"]]+dna1_halves[1]
                integ_funcs = [Polymer_transformation(result,circ1,parentsdict=pdict)]
            elif(circ2 ==False and circ1 == True):
                #if the sites are backwards just reverse everything
                return self.integrate(site2,site1)
            elif(circ1==False and circ1==circ2):
                #here we are recombining two linear plasmids, so two linear plasmids are produced

                result1 = dna1_halves[0]+[[prod1,"forward"]]+dna2_halves[1]
                result2 = dna2_halves[0]+[[prod2,"forward"]]+dna1_halves[1]
                integ_funcs = [Polymer_transformation(result1,parentsdict=pdict),Polymer_transformation(result2,parentsdict=pdict)]
        
        #newdna = [a.create_polymer(*dna_inputs) for a in integ_funcs]

        site1.linked_sites[site2] = [integ_funcs,[]]
        site2.linked_sites[site1] = [[a.inverted() for a in integ_funcs],[]]
        return integ_funcs





class Integrase_Enumerator(GlobalComponentEnumerator):
    def __init__(self,name:str,int_mechanisms = None):
        if(int_mechanisms is None):
            int_mechanisms={"int1":IntegraseMechanism()}
        self.int_mechanisms = int_mechanisms
        GlobalComponentEnumerator.__init__(self,name=name)
    def combine_dictionaries(self,dict1,dict2):
        """append lists that share the same key, and add new keys"""
        outdict = dict1
        for key in dict2:
            if key in outdict:
                assert(isinstance(dict2[key],list))
                assert(isinstance(outdict[key],list))
                outdict[key] += dict2[key]
            else:
                outdict[key] = dict2[key]
        return outdict
    def list_integrase(self,construct):
        """lists all the parts that can be acted on by integrases"""
        int_dict = {}
        for part in construct.parts_list:
            if(isinstance(part,AttachmentSite) and part.integrase is not None):
                part_integrase = part.integrase.name
                if(part_integrase in int_dict):
                    int_dict.update({part_integrase:int_dict[part_integrase]+[part]})
                else:
                    int_dict[part_integrase]=[part]
        return int_dict
    def enumerate_components(self,components = None, **keywords):
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
        
        construct_list = []
        for component in components:
            if(isinstance(component,DNA_construct)):
                construct_list += [component]
       
        int_dict = {}
        for construct in construct_list:
            #list each integrase that exists and which sites they react with
            con_dict = self.list_integrase(construct)
            
            int_dict = self.combine_dictionaries(int_dict,con_dict)
        constructlist = []
        for integrase in int_dict:
            if(integrase in self.int_mechanisms):
                int_mech = self.int_mechanisms[integrase]
                #now, going through each one, generate the reactions and species that arise
                attsites = int_dict[integrase]
                #but now we need to know what kind of integrase reactions are possible
                attcombos = [a for a in it.combinations(attsites,2)]
                #print(attcombos)
                for combo in attcombos:
                    #first question: is this combo legal?
                    if(tuple([a.site_type for a in combo]) in int_mech.reactions):
                        #this means the reaction can exist
                        int_functions = int_mech.integrate(combo[0],combo[1])
                        new_dnas = []
                        for a in int_functions:
                            new_dna = a.create_polymer(combo[0].parent,combo[1].parent)
                            new_dnas += [new_dna]
                        constructlist += new_dnas
        return constructlist

