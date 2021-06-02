#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.
from typing import Dict, List
from .species import Species, ComplexSpecies
from .polymer import OrderedPolymer, OrderedMonomer, NamedPolymer
from .dna_construct import Construct, DNA_construct
from .dna_part_misc import IntegraseSite
from .component_enumerator import GlobalComponentEnumerator
import itertools as it
from .utils import combine_dictionaries
import copy
import time
class Polymer_transformation:
    def __init__(self,partslist,circular=False,parentsdict = None,material_type="dna"):
        """A Polymer transformation is like a generic transformation of a polymer sequence.
        You specify a parts list that would make up the output polymer. This list can contain:
        parts from ordered polymers
        parts that aren't in any polymers (have parent = None)
        parts from ordered polymers are considered as "placeholders"
        parts with parent = None are inserted into the new polymer.

        Also you can specify if the output should be circular or not.

        example:

        valid partslist:
        partslist = [Monomer(forward,3,"input1"),Monomer(reverse,1,"input2"),Promoter(forward,None,None)]
        
        then, self.create_polymer([polymer1,polymer2])
        takes element 3 from polymer1 and puts it forward, element 1 from polymer 2, and a Promoter object and creates
        a new polymer by feeding these three monomers into a polymer constructor.

        new_polymer.parts_list = [polymer1[3].setdir("forward",polymer2[1].setdir("reverse"),[Promoter,"forward"]]
        """
        if(parentsdict is None):
            parentsdict = {} #the input to this function is a list of monomers that belong to various parents.
                             #each different parent that is represented in these monomers is converted into "input#" arbitrarily
                             #to keep it consistant, a link of that parent to the appropriate "input#" is kept in this dictionary.
                             #you can also pass a pre-populated dictionary into this function in order to control which parent gets
                             #which "input#". This is essential if we want to properly "reverse" a transformation (input1 becomes input2 and input2 becomes input1)
        inputcount = 1
        actual_partslist = []
        partdir = -1 #-1 = not set. Valid values are "forward" "reverse" and None
        #go through the parts
        inputname = None
        new_parentsdict = {}
        for part in partslist:
            if(type(part)==list or type(part) ==tuple):
                #if the part is a list that means it looks like [OrderedMonomer,"direction"]
                part_ref = part[0]
                partdir = part[1]
            else:
                part_ref = part
                partdir = part_ref.direction
            #if the parent is populated, it means this part should be a placeholder
            if(part_ref.parent is not None):
                #if we haven't tracked it already
                if((part_ref.parent not in new_parentsdict) and (part_ref.parent not in parentsdict)):
                    #create a 'blank' polymer based on the input, and give it a generic name
                    inputname = "input"+str(inputcount)
                    while(inputname in list(new_parentsdict.values())):
                        inputcount+=1
                        inputname = "input"+str(inputcount)
                    dummyPolymer = self.dummify(part_ref.parent,inputname)
                    #set the value in the dictionary
                    new_parentsdict[part_ref.parent]=dummyPolymer.name
                    #increment for the next time this happens
                    inputcount+=1
                    
                elif((part_ref.parent not in new_parentsdict) and (part_ref.parent in parentsdict)):
                    inputname = parentsdict[part_ref.parent]
                    dummyPolymer = self.dummify(part_ref.parent,inputname)
                    new_parentsdict[part_ref.parent]=dummyPolymer.name
                elif(inputname is None or inputname != new_parentsdict[part_ref.parent]):
                    inputname = new_parentsdict[part_ref.parent]
                    dummyPolymer = self.dummify(part_ref.parent,inputname)
                #this variable is the actual parts list that will be stored. It has mostly dummy parts
                if(part_ref.parent[part_ref.position]!=part_ref):
                    #this will be true only for situations where a component needs to be transformed into
                    #another component but still needs to keep a reference to where it was in the parent
                    copied_part = part_ref.get_orphan()
                    copied_part.parent = dummyPolymer
                    actual_partslist += [[copied_part,partdir]]
                else:
                    actual_partslist += [[dummyPolymer[part_ref.position],partdir]]
            else:
                #if the part has no parent, copy it and put it in
                actual_partslist += [[copy.deepcopy(part_ref),partdir]]
        self.number_of_inputs = len(new_parentsdict)
        self.parentsdict = new_parentsdict
        self.partslist = actual_partslist
        self.circular = circular
        self.material_type = material_type
    def renumber_output(self,output_renumbering_function):
        """change the ordering of the output list, using the output_renumbering_function which takes
        in an int and returns an int which is the new index of the part"""
        new_partslist = []
        for i in range(len(self.partslist)):
            new_i,direc = output_renumbering_function(i)
            part = self.partslist[new_i]
            if(direc == "r"):
                new_partslist+= [[part[0],["forward","reverse"][part[1]=="forward"]]]
            else:
                new_partslist+= [part]
        self.partslist = new_partslist

    def get_renumbered(self,output_renumbering_function):
        """return a copy of this transformation, but the output indexes are renumbered"""
        rxn_copied = copy.copy(self)
        rxn_copied.renumber_output(output_renumbering_function)
        return rxn_copied
    def reversed(self):
        """return a circularly permuted version of self. That means the inputs are shuffled around
        For example, we had input1, input2, input3. Now we will have input1=input2, input2=input3, input3=input1."""
        new_parentsdict = {}

        new_name_dict = {}
        for inputnum in range(self.number_of_inputs):
            if(inputnum == self.number_of_inputs-1):
                new_name_dict["input"+str(inputnum+1)]="input1"
            else:
                new_name_dict["input"+str(inputnum+1)]="input"+str(inputnum+2)
            
        if(len(self.parentsdict)==1):
            #in this case we have only "input1"
            return self
        else:
            for part in self.partslist:
                if(part[0].parent is not None and part[0].parent not in new_parentsdict):
                    #this next line outputs input2 when given input1
                    
                    newname = new_name_dict[part[0].parent.name]
                    new_parentsdict[part[0].parent]=newname
            return Polymer_transformation(self.partslist,self.circular,parentsdict=new_parentsdict)
    def create_polymer(self,polymer_list,**keywords):
        """this function creates a new polymer from the template saved inside this class.
        A polymer_list is a list of polymers from which the resulting polymer is made. Some of
        the parts which compose the output polymer don't have a parent, and therefore are new parts.
        Usually, these parts will have a "bound" form, which is basically the version of them which
        has proteins bound. if bound=False, the unbound form of these parts will be used."""
        polymer_dict = {"input"+str(a+1):b for a,b in enumerate(polymer_list)}
        assert(len(polymer_list)>=self.number_of_inputs)
        outlst = []
        for part_list in self.partslist:
            part = part_list[0]
            partdir = part_list[1]
            outpart = None
            
            if((part.parent is not None) and not hasattr(part,"name")):
                outpart = polymer_dict[part.parent.name][part.position] #grab the part from the proper input
            elif((part.parent is not None) and hasattr(part,"name")):
                #in this case we are transforming a part into a different part, taking the complexes from
                #the previous position
                
                if(isinstance(polymer_dict["input1"],Construct)):
                    outpart = part.get_orphan()
                else:
                    template_part = polymer_dict[part.parent.name][part.position]
                    if(isinstance(template_part,ComplexSpecies)):
                        #in this case we have to replace the thing inside the complex with the new part, but leave
                        #everything else the same.
                        old_species = template_part.get_species(recursive=True)
                        core_parts = []
                        for spec in old_species:
                            if(spec.material_type=="part" and not spec == template_part):#if you have a material type of part then you are
                                                           #part of a polymer and thus will be replaced
                                core_parts += [spec]
                        #now we should have found only one core part. If there are multiple core parts then
                        #we might have to get crazy.
                        new_part = part.dna_species
                        if(len(core_parts)==1):
                            new_part = template_part.replace_species(core_parts[0],part.dna_species)
                        else:
                            raise KeyError(f"{template_part} contained more than one species with material_type=\"part\", they were {core_parts}")
                        outpart = new_part
                    else:
                        outpart = part.dna_species
            else:
                #parts that don't come from input1 or input2.
                #they need to be either DNA_part or species objects, depending on
                #what kind of object we are making in the end.
                if(isinstance(polymer_dict["input1"],Construct)):
                    outpart = part
                else:
                    outpart = part.dna_species
            #assuming the stored parts have a valid direction
            if(hasattr(outpart,"linked_sites")):
                outpart = copy.copy(outpart)
                outpart.linked_sites = {} #make sure that any integrase sites we copy this way have no
                                            #linked sites, as those would not be links created by the integrate() function
            outlst += [[outpart,partdir]]
        if(hasattr(outlst[0][0],"material_type") and any(["complex" in a[0].material_type for a in outlst])):
            outpolymer = polymer_dict["input1"].__class__(outlst,circular = self.circular,material_type = "ordered_polymer")
        else:
            outpolymer = polymer_dict["input1"].__class__(outlst,circular = self.circular,material_type=self.material_type)
        return outpolymer
    @classmethod
    def dummify(cls,in_polymer,name):
        """creates a simplified polymer that has the same number of monomers, direction of monomers,
        and name as the input polymer, but is otherwise disconnected."""
        #this is used specifically with polymerTransformation. Dummified version of polymers are stored in
        #polymerTransformation as generic "slots" for monomers to be properly placed into
        out_list = []    
        for element in in_polymer:
            out_list += [OrderedMonomer(direction=element.direction)]
        circular = False
        if(hasattr(in_polymer,"circular")):
            circular = in_polymer.circular
        return NamedPolymer(out_list,name,circular=circular)
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
        if(self.circular):
            out_txt+="(circular)"
        return out_txt

class IntegraseRule:
    def __init__(self,name=None,reactions=None,allow_deletion=True,allow_integration=True,allow_inversion=True):
        """The integrase mechanism is a mechanism at the level of DNA. It creates DNA species which
        the integrase manipulations would lead to. This mechanism does not create any reaction rates.
        We need to figure out how integrase binding will work before being able to create
        reactions and their corresponding rates"""
        if(reactions is None):
            reactions = {("attB","attP"):"attL",("attP","attB"):"attR"}
        if(name is None):
            self.name = "int1"
        else:
            self.name = name
        self.allow_deletion = allow_deletion
        self.allow_integration = allow_integration
        self.allow_inversion = allow_inversion
        self.integrase_species = Species(self.name,material_type="protein")
        self.reactions = reactions
        self.attsites = []
        for reaction in reactions:
            self.attsites+=list(reaction)+[reactions[reaction]]
        self.attsites = list(set(self.attsites))
        self.integrations_to_do = [] #these are the reactions that will be performed at compile time
    def binds_to(self):
        return self.attsites
    def reaction_allowed(self,site1,site2):
        assert(isinstance(site1,IntegraseSite))
        assert(isinstance(site2,IntegraseSite))
        assert site1.integrase==site2.integrase
        assert site1.integrase==self.integrase_species
        if(tuple([site1.site_type,site2.site_type]) in self.reactions):
            return True
        return False
    def reactive_sites(self):
        """returns a list of attachment sites (strings) that participate in integrase reactions"""
        attsites = []
        for reaction in self.reactions:
            attsites+=list(reaction)
        attsites = list(set(attsites))
        return attsites
    def generate_products(self,site1,site2):
        """generates DNA_part objects corresponding to the products of recombination"""
        #the sites should have the same integrase and dinucleotide, otherwise it won't work
        assert (site1.integrase == site2.integrase)
        assert(site1.integrase==self.integrase_species), f"sites have integrase {site1.integrase} but should be {self.integrase_species}"
        assert(site1.dinucleotide == site2.dinucleotide)
        integrase = site1.integrase
        dinucleotide = site1.dinucleotide
        rxn1 = (site1.site_type,site2.site_type)
        rxn2 = (site2.site_type,site1.site_type)
        #this next part checks if these parts can even react
        if(rxn1 in self.reactions):
            prod1 = self.reactions[rxn1]
        else:
            raise KeyError("{} not found to react with {} in {}".format(site1,site2,self.reactions))
        
        if(rxn2 in self.reactions):
            prod2 = self.reactions[rxn2]
        else:
            raise KeyError("{} not found to react with {} in {}".format(site2,site1,self.reactions))
            
        
        part_prod1 = IntegraseSite(prod1,prod1,dinucleotide=dinucleotide,integrase=integrase,
                            direction=site1.direction,integrase_binding=site1.integrase_binding,
                            material_type = site1.material_type)
        
        part_prod2 = IntegraseSite(prod2,prod2,dinucleotide=dinucleotide,integrase=integrase,
                            direction=site2.direction,integrase_binding=site1.integrase_binding,
                            material_type = site2.material_type)
        if(site1.direction=="forward"):
            part_prod1.position = site1.position
            part_prod1.parent = site1.parent
            part_prod2.position = site2.position
            part_prod2.parent = site2.parent
            return (part_prod1,part_prod2)
        else:
            part_prod2.direction = site1.direction
            part_prod2.position = site1.position
            part_prod2.parent = site1.parent
            part_prod1.direction = site2.direction
            part_prod1.position = site2.position
            part_prod1.parent = site2.parent
            
            return(part_prod2,part_prod1)
    
    def integrate(self,site1,site2,also_inter=True,force_inter=False, existing_dna_constructs = None):
        """perform an integration reaction between the chosen sites and make new DNA_constructs
        site1 and site2 are integrase site dna_parts which have parents that are DNA_constructs.
        
        There are four possible reactions:
        1) inversion
           two sites are part of the same dna construct
           the result is another dna construct with the same circularity and the region in between the sites flipped
        2) deletion
           two sites are part of the same dna construct
           the result is two dna constructs: one with the same circularity but the region between the sites deleted, and another
           circular dna construct that contains the deleted portion
        3) integration
           the sites are on two different dna constructs
           the result is a single dna construct
        4) recombination
           the sites are on two different dna constructs
           the results are two different dna constructs with the proper portions swapped
        
        after the correct dna constructs are generated, the reactions which were done to produce them
        are encoded into polymer_transformations and "baked into" the integrase sites themselves. So,
        each integrase site knows which specific integrase reactions it should produce when it comes
        time to update_reactions.
        
        also_inter controls whether intramolecular reactions should also generate intermolecular reactions
        that occur between two copies of the same plasmid.

        force_inter forces a reaction to be intermolecular even if the two sites are on the same plasmid
        """
        intermolecular = True #by default, the reaction is intermolecular
        #if one of the sites is not part of a construct then raise an error!
        integ_funcs = []
        new_dna_constructs = [] #new dna constructs made by this function!
        if(existing_dna_constructs is None):
            existing_dna_constructs = []
        if(not(isinstance(site1.parent,Construct))):
            raise ValueError("{} not part of a construct".format(site1))
        elif(not(isinstance(site2.parent,Construct))):
            raise ValueError("{} not part of a construct".format(site2))
        cutpos1 = site1.position
        cutpos2 = site2.position
        #below are the references to the sites in the products
        dna_inputs = []
        if(site1.parent==site2.parent and not force_inter):
            intermolecular = False
            #these sites are part of the same piece of DNA, so they are going to do an intramolecular reaction
            contains_no_inter = any(['no_inter' in a.attributes for a in site1.parent])
            if(also_inter and not contains_no_inter):
                #we should generate the intermolecular reaction also!
                #in this case we are generating multiple integration reactions at the same time
                #i think the right thing to do is NOT return it?
                new_dna_constructs += self.integrate(site1,site2,force_inter=True,existing_dna_constructs = existing_dna_constructs)
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
            circularity = dna.circular
            
            if(site1.direction == site2.direction):
                if(self.allow_deletion):
                    #case 2: deletion
                    #if the sites point in the same direction, then we are doing a deletion reaction
                    #direction doesn't matter so we don't need to flip anything
                    cutdna_list_parts = list(dna[:cutpos1])+[[prod1,dna[cutpos1].direction]]+list(dna[cutpos2+1:]) #delete
                    newdna_list_parts = [[prod2,dna[cutpos2].direction]]+list(dna[1+cutpos1:cutpos2])

                    integ_funcs += [Polymer_transformation(cutdna_list_parts,circular = circularity),\
                                            Polymer_transformation(newdna_list_parts,circular=True)]
            else:
                if(self.allow_inversion):
                    #case 1: inversion
                    #this means we are dealing with an inversion
                    inv_segment = []
                    
                    [[a,] for a in dna[cutpos1+1:cutpos2][::-1]]
                    for a in dna[cutpos1+1:cutpos2][::-1]:
                        inv_segment += [[a,["forward","reverse"][a.direction=="forward"]]]
                    #the inverted segment is reversed
                    
                    invertdna_list = list(dna[:cutpos1])+\
                                    [[prod1,dna[cutpos1].direction]]+\
                                    inv_segment+ \
                                    [[prod2,dna[cutpos2].direction]]+\
                                    list(dna[cutpos2+1:])

                    integ_funcs += [Polymer_transformation(invertdna_list,circular=circularity)]
        else:
            if(self.allow_integration):
                #otherwise these sites are on different pieces of DNA, so they are going to combine
                dna1 = site1.parent
                dna2 = site2.parent
                circ1 = dna1.circular
                circ2 = dna2.circular
                if(dna1 == dna2):
                    #this will happen if we trying to do an intermolecular reaction between two copies of the same thing
                    dna2 = copy.copy(dna1)
                    dna2.name = dna2.name+"_duplicate"
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

                
                
                if(site1.direction=="reverse" and site2.direction=="reverse"):
                    prod1,prod2 = self.generate_products(site2,site1)
                else:
                    prod1,prod2 = self.generate_products(site1,site2)
                #direction of everything should be forward

                if(circ2==True):
                    #case 3: integration
                    #in this case we are combining a circular plasmid with a circular or linear plasmid
                    #either way the result is basically the same, except the result is either linear or circular
                    #result is ONE PIECE OF DNA
                    result = dna1_halves[0]+[[prod1,"forward"]]+dna2_halves[1]+dna2_halves[0]+[[prod2,"forward"]]+dna1_halves[1]
                    integ_funcs += [Polymer_transformation(result,circ1,parentsdict=pdict)]
                elif(circ2 ==False and circ1 == True):
                    #if the sites are backwards just reverse everything
                    new_dna_constructs += self.integrate(site2,site1,force_inter=force_inter,existing_dna_constructs=existing_dna_constructs)
                    #the above already populates the sites, so then we don't need to
                elif(circ1==False and circ1==circ2):
                    #case 4: recombination
                    #here we are recombining two linear dnas, so two linear dnas are produced
                    result1 = dna1_halves[0]+[[prod1,"forward"]]+dna2_halves[1]
                    result2 = dna2_halves[0]+[[prod2,"forward"]]+dna1_halves[1]
                    integ_funcs += [Polymer_transformation(result1,parentsdict=pdict),Polymer_transformation(result2,parentsdict=pdict)]
        if(len(integ_funcs)>0):
            for integ_func in integ_funcs:
                #generate new dna constructs and check if we already made them before
                created_dna = integ_func.create_polymer([site1.parent,site2.parent])
                output = Integrase_Enumerator.find_dna_construct(created_dna,existing_dna_constructs+new_dna_constructs)
                if(output is not None):
                    #if the construct has already been made, then fix the integ_func so it outputs the right thing
                    integ_func.renumber_output(output[1])
                else:
                    #otherwise, add it to the list
                    new_dna_constructs += [created_dna]
            #link the two sites and give them the adjusted integ_funcs
            site1.linked_sites[(site2,intermolecular)] = [integ_funcs,[]]
            site2.linked_sites[(site1,intermolecular)] = [[a.reversed() for a in integ_funcs],[]]
        #the return value of this function is used mostly only for generating constructs
        return new_dna_constructs





class Integrase_Enumerator(GlobalComponentEnumerator):
    def __init__(self,name:str,int_mechanisms = None):
        if(int_mechanisms is None):
            int_mechanisms={"int1":IntegraseRule()}
        self.int_mechanisms = int_mechanisms
        GlobalComponentEnumerator.__init__(self,name=name)
    
    def list_integrase(self,construct):
        """lists all the parts that can be acted on by integrases"""
        int_dict = {}
        for part in construct.parts_list:
            if(isinstance(part,IntegraseSite) and part.integrase is not None):
                part_integrase = part.integrase.name
                if(part_integrase in int_dict):
                    int_dict.update({part_integrase:int_dict[part_integrase]+[part]})
                else:
                    int_dict[part_integrase]=[part]
        return int_dict
    def reset(self,components=None, **keywords):
        """this resets the linked_sites member in any attachment sites"""
        for component in components:
            if(hasattr(component,"parts_list")):
                for part in component:
                    if(hasattr(part,"linked_sites")):
                        part.linked_sites = {}
    @classmethod
    def find_dna_construct(cls,construct:Construct,conlist:List[Construct]):
        """find a construct that matches the input "construct", but can be reverse or circularly permuted
        returns: found_construct, index_function(index) => (new_index,"f" or "r" if it must be reversed)
        or,
        None, if no matching construct is found"""
        matched_construct = None
        for other_construct in conlist:
            if(not isinstance(other_construct,Construct)):
                continue
            if(construct.directionless_hash==other_construct.directionless_hash):
                if(matched_construct is not None):
                    #a construct must not match two constructs, since they are generated and checked in order
                    raise KeyError(f"{construct} matches with {matched_construct} but also {other_construct}")
                matched_construct = other_construct
        if(matched_construct is None):
            return None
        other_construct = matched_construct
        other_indexes = list(range(len(other_construct)))
        if(construct.circular):
            for pos,_ in enumerate(construct):
                cp_construct = construct.get_circularly_permuted(pos)
                
                if(cp_construct.get_species() == other_construct.get_species()):
                    return other_construct, (lambda a:((other_indexes[pos:]+other_indexes[:pos])[a],"f"))
                elif(cp_construct.get_reversed().get_species()==other_construct.get_species()):
                    return other_construct, (lambda a:((other_indexes[pos:]+other_indexes[:pos])[::-1][a],"r"))
        else:
            if(construct.get_species() == other_construct.get_species()):
                return other_construct, (lambda a: (a,"f"))
            elif(construct.get_reversed().get_species() == other_construct.get_species()):
                return other_construct, (lambda a: (other_indexes[::-1][a],"r"))

    def enumerate_components(self,components = None,previously_enumerated = None, **keywords):
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
        if(previously_enumerated is None):
            previously_enumerated = []
        construct_list = []
        if(previously_enumerated is None):
            previously_enumerated = []
        for component in components:
            if(isinstance(component,DNA_construct)):
                construct_list += [component]
        int_dict = {}
        for construct in construct_list:
            #list each integrase that exists and which sites they react with
            
            con_dict = self.list_integrase(construct)
            
            int_dict = combine_dictionaries(int_dict,con_dict)
        constructlist = []
        for integrase in int_dict:
            if(integrase in self.int_mechanisms):
                int_mech = self.int_mechanisms[integrase]
                #now, going through each one, generate the reactions and species that arise
                attsites = int_dict[integrase]
                #but now we need to know what kind of integrase reactions are possible
                reactive_sites = int_mech.reactive_sites()
                attcombos = [a for a in it.combinations(attsites,2) if ((a[0].site_type in reactive_sites) and (a[1].site_type in reactive_sites))]
                for combo in attcombos:
                    #first question: is this combo legal?
                    if(int_mech.reaction_allowed(combo[0],combo[1])):
                        #this means the reaction can exist
                        #integrate now
                        new_dnas = int_mech.integrate(combo[0],combo[1],existing_dna_constructs = (previously_enumerated+constructlist))
                        constructlist += new_dnas
                
        return constructlist

