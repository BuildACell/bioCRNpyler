from .mechanism import *
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer

class Reversible_Bimolecular_Binding(Mechanism):
    def __init__(self, name="reversible_bimolecular_binding",
                 mechanism_type="bimolecular_binding"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, s1, s2, **keywords):
        complexS = ComplexSpecies([s1, s2])
        return [s1, s2, complexS]

    def update_reactions(self, s1, s2, component = None, kb = None, ku = None, \
                                              part_id = None,complex=None, **keywords):

        #Get Parameters
        if part_id is None:
            repr(s1)+"_"+repr(s2)
        if kb is None and component is not None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku is None and component is not None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if component is None and (kb is None or ku is None):
            raise ValueError("Must pass in a Component or values for kb, ku.")
        if(complex==None):
            complexS = ComplexSpecies([s1, s2])
        else:
            complexS = complex
        rxns = [Reaction([s1, s2], [complexS], k=kb, k_rev=ku)]
        return rxns


class One_Step_Cooperative_Binding(Mechanism):
    """A reaction where n binders (A) bind to 1 bindee (B) in one step
       n A + B <--> nA:B
    """
    def __init__(self, name="one_step_cooperative_binding",
                 mechanism_type="cooperative_binding"):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, binder, bindee, complex_species = None, cooperativity=None, component = None, part_id = None, **kwords):

        if part_id is None:
            part_id = repr(binder)+"-"+repr(bindee)

        if cooperativity is None and component != None:
            cooperativity = component.get_parameter("cooperativity", part_id = part_id, mechanism = self)
        elif component is None and cooperativity is None:
            raise ValueError("Must pass in a Component or values for cooperativity")

        complexS = None
        if complex_species is None:
            complex_name = None
            material_type = None
        elif isinstance(complex_species, str):
            complex_name = complex_species
            material_type = None
        elif isinstance(complex_species, Species):
            complexS = complex_species
            material_type = complex_species.material_type
        else:
            raise TypeError("complex_species keyword must be a str, Species, or None.")

        if complexS is None:
            complexS = ComplexSpecies([binder]*int(cooperativity)+[bindee], name = complex_name, material_type = material_type)

        
        return [binder, bindee, complexS]

    def update_reactions(self, binder, bindee, complex_species = None, component = None, kb = None, ku = None, part_id = None, cooperativity=None, **kwords):

        complexS = self.update_species(binder, bindee, cooperativity = cooperativity, part_id = part_id, component = component, **kwords)[-1]

        #Get Parameters
        if part_id is None:
            part_id = repr(binder)+"-"+repr(bindee)
        if kb is None and component is not None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        if ku is None and component is not None:
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        if cooperativity is None and component is not None:
            cooperativity = component.get_parameter("cooperativity", part_id = part_id, mechanism = self)
        if component is None and (kb is None or ku is None or cooperativity is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")


        rxns = []
        rxns += [
            Reaction(inputs=[binder, bindee], outputs=[complexS],
                     input_coefs=[cooperativity, 1], output_coefs=[1], k=kb,
                     k_rev=ku)]
        return rxns


class Two_Step_Cooperative_Binding(Mechanism):
    """A reaction where n binders (s1) bind to 1 bindee (s2) in two steps
       n A <--> nx_A
       nx_A <--> nx_A:B
    """
    def __init__(self, name="two_step_cooperative_binding",
                 mechanism_type="cooperative_binding"):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, binder, bindee, component = None, complex_species = None, n_mer_species = None, cooperativity=None, part_id = None, **keywords):

        if part_id is None:
            part_id = repr(binder)+"-"+repr(bindee)

        if cooperativity is None and component != None:
            cooperativity = component.get_parameter("cooperativity", part_id = part_id, mechanism = self)
        elif component is None and cooperativity is None:
            raise ValueError("Must pass in a Component or values for cooperativity")

        n_mer = None
        if n_mer_species is None:
            n_mer_name = binder.name
            n_mer_material = binder.material_type
        elif isinstance(n_mer_species, str):
            n_mer_name = n_mer_species
            n_mer_material = "complex"
        elif isinstance(n_mer_species, Species):
            n_mer = n_mer_species
        else:
            raise TypeError("n_mer_species keyword nust be a str, Species, or None. Not "+str(n_mer_species))

        if n_mer is None:
            n_mer = Multimer(binder, cooperativity, name = n_mer_name, material_type = n_mer_material)

        complexS = None
        if complex_species is None:
            complex_name = None
            material_type = "complex"
        elif isinstance(complex_species, str):
            complex_name = complex_species
            material_type = "complex"
        elif isinstance(complex_species, Species):
            complexS = complex_species
        else:
            raise TypeError("complex_species keyword must be a str, Species, or None. Not "+str(complex_species))

        if complexS is None:
            complexS = ComplexSpecies([n_mer, bindee], name = complex_name)
        return [binder, bindee, complexS, n_mer]

    def update_reactions(self, binder, bindee, kb = None, ku = None, component = None, part_id = None, cooperativity=None, complex_species = None, n_mer_species = None, **keywords):
        """
        Returns reactions:
        cooperativity binder <--> n_mer, kf = kb1, kr = ku1
        n_mer + bindee <--> complex, kf = kb2, kr = ku2
        :param s1:
        :param s2:
        :param kb:
        :param ku:
        :param cooperativity:
        :param keywords:
        :return:
        """

        if part_id is None:
            repr(binder)+"-"+repr(bindee)
        if (kb is None or ku is None or cooperativity is None) and Component != None:
            kb1 = component.get_parameter("kb1", part_id = part_id, mechanism = self)
            kb2 = component.get_parameter("kb2", part_id = part_id, mechanism = self)
            ku1 = component.get_parameter("ku1", part_id = part_id, mechanism = self)
            ku2 = component.get_parameter("ku2", part_id = part_id, mechanism = self)
            cooperativity = component.get_parameter("cooperativity", part_id = part_id, mechanism = self)
        elif component is None and (kb is None or ku is None or cooperativity is None):
            raise ValueError("Must pass in a Component or values for kb, ku, and cooperativity")
        elif len(kb) != len(ku) != 2:
            raise ValueError("kb and ku must contain 2 values each for "
                             "two-step binding")
        else:
            kb1, kb2 = kb
            ku1, ku2 = ku
        n_mer_name = f"{cooperativity}x_{binder.material_type}_{binder.name}"
        n_mer = ComplexSpecies([binder], name = n_mer_name)

        binder, bindee, complexS, n_mer = self.update_species(binder, bindee, component = component, complex_species = complex_species, n_mer_species = n_mer_species, cooperativity=cooperativity, part_id = part_id, **keywords)

        rxns = [
            Reaction(inputs=[binder], outputs=[n_mer],
                     input_coefs=[cooperativity], output_coefs=[1], k=kb1,
                     k_rev=ku1),
            Reaction(inputs=[n_mer, bindee], outputs=[complexS], k=kb2,
                     k_rev=ku2)]

        return rxns

class Combinatorial_Cooperative_Binding(Mechanism):
    """a reaction where some number of binders bind combinatorially to a bindee"""
    def __init__(self,name="Combinatorial_Cooperative_binding",
                            mechanism_type="cooperative_binding"):
        Mechanism.__init__(self,name,mechanism_type)
    def make_cooperative_complex(self,combo,bindee,cooperativity):
        """given a list of binders and their cooperativities, make a complex
        that contains the binders present N number of times where N is
        each one's cooperativity"""
        complexed_species_list = []
        for binder in combo:
            binder_cooperativity = int(cooperativity[binder.name])
            #I hope that cooperativity is an int! what if it isn't
            if(binder_cooperativity > 1):
                complexed_species_list += [binder]*binder_cooperativity
            else:
                complexed_species_list += [binder]
        complexed_species_list += [bindee]
        if(len(complexed_species_list)==1):
            myspecies = complexed_species_list[0]
        else:
            myspecies = ComplexSpecies(complexed_species_list)
        return myspecies

    def update_species(self,binders,bindee,cooperativity=None,\
                            component = None, part_id = None, **kwords):
        cooperativity_dict = {}
        for binder in binders:
            binder_partid = part_id+"_"+binder.name
            if ((cooperativity is None) or (type(cooperativity)==dict and binder_partid not in cooperativity) \
                                                                                        and (component is not None)):
                #here we are extracting the relevant cooperativity value from the dictionary which should be passed
                #in as the cooperativity argument
                coop_val = component.get_parameter("cooperativity", part_id = binder_partid, mechanism = self)
            elif type(cooperativity)==dict and binder_partid in cooperativity:
                coop_val = cooperativity[binder_partid]
            if component is None and ( cooperativity is None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            cooperativity_dict[binder.name]=coop_val
        out_species = []
        for i in range(1, len(binders)+1):
            for combo in it.combinations(binders,i):
                #go through every possible combination of reactants and dna and make
                #all the complexes
                out_species += [self.make_cooperative_complex(combo,bindee,cooperativity_dict)]
        return out_species
    
    def update_reactions(self,binders,bindee,component=None,kbs=None,kus=None,\
                                                            part_id = None,cooperativity=None,**kwords):
        binder_params = {}
        for binder in binders:
            binder_partid = part_id+"_"+binder.name
            if ((type(kbs)==dict and binder not in kbs) or (type(kbs)!=dict and component is not None)):
                kb = component.get_parameter("kb", part_id = binder_partid, mechanism = self)
            elif(type(kbs)==dict and binder in kbs):
                kb = kbs[binder.name]
            elif(type(kbs)!=dict and component is None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            if ((type(kus)==dict and binder not in kus) or (kus is None and component is not None)):
                ku = component.get_parameter("ku", part_id = binder_partid, mechanism = self)
            elif(type(kus)==dict and binder in kus):
                ku = kus[binder.name]
            elif(type(kus)!=dict and component is None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            if ((cooperativity is None) or (type(cooperativity)==dict and binder.name not in cooperativity)  \
                                                                                        and component is not None):
                coop_val = component.get_parameter("cooperativity", part_id = binder_partid, mechanism = self)
            elif type(cooperativity)==dict and binder.name in cooperativity:
                coop_val = cooperativity[binder.name]
            if component is None and (kb is None or ku is None or cooperativity is None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            binder_params[binder] = {"kb":kb,"ku":ku,"cooperativity":coop_val}
        #out_rxns = []
        rxndict = {}
        coop_dict = {a.name:binder_params[a]["cooperativity"] for a in binder_params}
        for i in range(1, len(binders)+1):
            for combo in it.combinations(binders,i):
                #come up with all combinations of binders
                
                product = self.make_cooperative_complex(combo,bindee,coop_dict)
                #this is the complex which becomes the product
                for binder in combo:
                    #then in each case, one binder can be added to make
                    #this combination
                    reactant = tuple(set(combo)-set([binder]))
                    rxn_prototype = (binder,reactant)
                    
                    #this part makes a describer of the reaction; which reactants are combining?
                    #print(rxn_prototype)
                    #print(rxndict)
                    if(rxn_prototype in rxndict):
                        #if we already did this reaction then forget about it
                        continue
                    else:
                        reactant_complex = self.make_cooperative_complex(reactant,bindee,coop_dict)
                        reaction = Reaction(inputs=[binder, reactant_complex], outputs=[product],
                        input_coefs=[binder_params[binder]["cooperativity"], 1], output_coefs=[1], \
                                         k=binder_params[binder]["kb"],k_rev=binder_params[binder]["ku"])
                        rxndict[rxn_prototype]=reaction
        return [rxndict[a] for a in rxndict]

        
class One_Step_Binding(Mechanism):
    def __init__(self, name="one_step_binding",
                 mechanism_type="binding"):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, species, component = None, complex_species = None, part_id = None, **keywords):
        if part_id is None:
            part_id = ""
            for s in species:
                part_id += s.name+"_"
            part_id = part_id[:-1]

        if complex_species is None:
            complex_species = ComplexSpecies(species)

        return species + [complex_species]


    def update_reactions(self, species, component = None, complex_species = None, part_id = None, kb = None, ku = None, **keywords):
        if part_id is None:
            part_id = ""
            for s in species:
                part_id += s.name+"_"
            part_id = part_id[:-1]

        if (kb is None or ku is None) and Component != None:
            kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
            ku = component.get_parameter("ku", part_id = part_id, mechanism = self)
        elif component is None and (kb is None or ku is None):
            raise ValueError("Must pass in a Component or values for kb and ku")

        if complex_species is None:
            complex_species = ComplexSpecies(species)

        return [Reaction(inputs = species, outputs = [complex_species], k = kb, k_rev = ku)]


    


