# mechanism.py - mechanism class for implementing TX-TL mechanisms
# RMM, 16 Aug 2018
#
# Mechanisms the means by which all reactions in a TX-TL reaction are
# established.  Mechanisms can be overridden to allow specialized
# processing of core reactions (eg, adding additional detail, using
# simplified models, etc.
#
# Mechanisms are established in the following order (lower levels
# override higher levels):
#
# Default extract mechanisms
#   Default mechanisms
#       Mechanisms passed to Component() [eg DNA Assembly]
#         Mechanisms based to Sub) [eg, DNA elements]
#
# This hierarchy allows reactions to be created without the user
# having to specify any alternative mechanisms (defaults will be
# used), but also allows the user to override all mechanisms used for
# every (e.g, by giving an alternative transcription
# mechanisms when setting up an extract) or individual mechanisms for
# a given (by passing an alternative mechanism just to that
# .
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from warnings import warn
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer
from .component import Component
import itertools as it

class Mechanism(object):
    """Mechanism class for core mechanisms

    Core mechanisms within a mixture (transcription, translation, etc)

    The Mechanism class is used to implement different core
    mechanisms in TX-TL.  All specific core mechanisms should be
    derived from this class.

    """
    def __init__(self, name, mechanism_type=""):
        self.name = name
        self.mechanism_type = mechanism_type
        if mechanism_type == "" or mechanism_type is None:
            warn(f"Mechanism {name} instantiated without a type. This could "
                 "prevent the mechanism from being inherited properly.")

    def update_species(self):
        """
        the child class should implement this method
        :return: empty list
        """
        warn(f"Default Update Species Called for Mechanism = {self.name}.")
        return []

    def update_reactions(self, component = None, part_id = None):
        """
        the child class should implement this method
        :return: empty list
        """
        warn(f"Default Update Reactions Called for Mechanism = {self.name}.")
        return []

    def __repr__(self):
        return self.name


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
            if ((cooperativity == None) or (type(cooperativity)==dict and binder_partid not in cooperativity) \
                                                                                        and (component is not None)):
                #here we are extracting the relevant cooperativity value from the dictionary which should be passed
                #in as the cooperativity argument
                coop_val = component.get_parameter("cooperativity", part_id = binder_partid, mechanism = self)
            elif type(cooperativity)==dict and binder_partid in cooperativity:
                coop_val = cooperativity[binder_partid]
            if component is None and ( cooperativity == None):
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
            elif(type(kbs)!=dict and component == None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            if ((type(kus)==dict and binder not in kus) or (kus == None and component is not None)):
                ku = component.get_parameter("ku", part_id = binder_partid, mechanism = self)
            elif(type(kus)==dict and binder in kus):
                ku = kus[binder.name]
            elif(type(kus)!=dict and component == None):
                raise ValueError("Must pass in a Component or values for kb, ku, and coopertivity.")
            if ((cooperativity == None) or (type(cooperativity)==dict and binder.name not in cooperativity)  \
                                                                                        and component is not None):
                coop_val = component.get_parameter("cooperativity", part_id = binder_partid, mechanism = self)
            elif type(cooperativity)==dict and binder.name in cooperativity:
                coop_val = cooperativity[binder.name]
            if component is None and (kb == None or ku == None or cooperativity == None):
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
class EmptyMechanism(Mechanism):
    def __init__(self, name, mechanism_type):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, **keywords):
        return []

    def update_reactions(self, **keywords):
        return []