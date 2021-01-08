
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import itertools as it
from warnings import warn

from .dna_part import DNA_part
from .mechanisms_binding import (Combinatorial_Cooperative_Binding,
                                 One_Step_Cooperative_Binding)
from .mechanisms_txtl import (NegativeHillTranscription,
                              PositiveHillTranscription)
from .species import Species


class Promoter(DNA_part):
    """A basic Promoter class with no regulation. Needs to be included in a DNAassembly or DNAconstruct to function.
    """
    def __init__(self, name, assembly=None,
                 transcript=None, length=0,
                 mechanisms=None, parameters=None, protein=None,dna_to_bind=None, **keywords):
        self._dna_bind = dna_to_bind
        self.length = length

        if transcript is None and assembly is None:
            self.transcript = None
        elif transcript is None:
            self.transcript = Species(assembly.name, material_type="rna")
        elif type(transcript) == list:
            self.transcript = transcript
        else:
            self.transcript = self.set_species(transcript, material_type='rna')
        
        if(isinstance(protein,str)):
            self.protein = [self.set_species(protein,material_type="protein")]
        elif(isinstance(protein,Species)):
            self.protein = [protein]
        elif(isinstance(protein,list)):
            self.protein = protein
        else:
            self.protein = None
        #Promoter should not have initial conditions. These need to be in DNAAssembly or DNAConstruct
        if "initial_concentration" in keywords.values() and keywords["initial_concentration"] is not None:
            raise AttributeError("Cannot set initial_concentration of a Promoter. Must set initial_concentration for the DNAassembly or DNAConstruct.")
        if "initial_condition_dictionary" in keywords.values() and keywords["initial_condition_dictionary"] is not None:
            raise AttributeError("Cannot set initial_condition_dictionary of a Promoter. Must set initial_condition_dictionary for the DNAassembly or DNAconstruct.")

        DNA_part.__init__(self, name = name, mechanisms = mechanisms,
                           parameters = parameters, assembly=assembly, **keywords)
    def update_species(self):
        mech_tx = self.get_mechanism("transcription")
        species = []

        species += mech_tx.update_species(dna = self.dna_to_bind, \
            transcript = self.transcript, protein = self.get_protein_for_expression(),
            component = self, part_id = self.name)

        return species
    @property
    def dna_to_bind(self):
        if(self._dna_bind is None):
            if(self.assembly is None):
                return None
            else:
                self._dna_bind = self.assembly.dna
        return self._dna_bind
    @dna_to_bind.setter
    def dna_to_bind(self,value):
        self._dna_bind = value
    def get_species(self):
        return None
    def update_reactions(self):
        mech_tx = self.get_mechanism("transcription")
        reactions = []

        reactions += mech_tx.update_reactions(dna = self.dna_to_bind, \
                        component = self, part_id = self.name, complex = None,
                        transcript = self.transcript, protein = self.get_protein_for_expression())
        return reactions
    def update_component(self,dna,rnas,proteins,mypos = None):
        """returns a copy of this component, except with the proper fields updated"""
        if(dna.material_type == "rna"):
            #Promoters only work with DNA
            return None
        out_component = copy.deepcopy(self)
        if(mypos is not None):
            out_component.dna_to_bind = dna[mypos]
        else:
            out_component.dna_to_bind = dna
        myrna = rnas[self]
        out_component.transcript = myrna.get_species()
        rbslist = proteins[myrna]
        proteinlist = []
        for a in rbslist:
            proteinlist += rbslist[a]
        out_component.protein = proteinlist
        return out_component


    #Used for expression mixtures where transcripts are replaced by proteins
    def get_protein_for_expression(self):
        if self.transcript is None:
            return self.protein
        else:
            return None

    @classmethod
    def from_promoter(cls, name, assembly, transcript, protein):
        """Helper function to initialize a promoter instance from another promoter or str.

        :param name: either string or an other promoter instance
        :param assembly:
        :param transcript:
        :param protein:
        :return: Promoter instance
        """
        if isinstance(name, Promoter):
            promoter_instance = copy.deepcopy(name)
            promoter_instance.assembly = assembly
            promoter_instance.transcript = transcript
            promoter_instance.protein = protein
        elif isinstance(name, str):
            promoter_instance = cls(name=name, assembly=assembly,
                                    transcript=transcript, protein=protein)
        else:
            raise TypeError(f'Promoter can be initialized from string or another promoter! We got {type(name)}')
        return promoter_instance


class RegulatedPromoter(Promoter):
    """
    A Promoter class with simple regulation.
    regulators = [list of species]
    Each regulator binds independently to the Promoter to regulate it.
    """
    def __init__(self, name: str, regulators, leak=True, assembly=None,
                 transcript=None, length=0, mechanisms=None,
                 parameters=None , **keywords):

        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters,
                          **keywords)

        if not isinstance(regulators, list):
            regulators = [regulators]

        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]

        self.leak = leak

        self.add_mechanism(One_Step_Cooperative_Binding(), "binding")
        self.complexes = []
        
    def update_species(self):
        mech_tx = self.get_mechanism("transcription")
        mech_b = self.get_mechanism('binding')
        species = []

        
        self.complexes = []
        if self.leak is not False:
            species += mech_tx.update_species(dna = self.dna_to_bind, transcript = self.transcript,\
                             component = self, part_id = self.name+"_leak", protein = self.get_protein_for_expression())

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]

            species_b = mech_b.update_species(regulator, self.dna_to_bind, part_id = self.name+"_"+regulator.name,\
                             component = self, protein = self.get_protein_for_expression())
            species += species_b

            #Find complexes containing DNA and the regulator
            #DNA should be *not* a part of an OrderedPolymer for this to work
            dna_simple = copy.deepcopy(self.dna_to_bind)
            dna_simple.remove()
            for s in species_b:
                if dna_simple in s and regulator in s:
                    self.complexes += [s]

                    species += mech_tx.update_species(dna = s, transcript = self.transcript, \
                            protein = self.get_protein_for_expression(), part_id = self.name+"_"+regulator.name, component = self)
        return species

    def update_reactions(self):
        reactions = []
        mech_tx = self.get_mechanism("transcription")
        mech_b = self.get_mechanism('binding')

        if self.leak is not False:
            reactions += mech_tx.update_reactions(dna = self.dna_to_bind, component = self, part_id = self.name+"_leak", \
                                                            transcript = self.transcript, protein = self.get_protein_for_expression())

        for i in range(len(self.regulators)):
            regulator = self.regulators[i]
            complex_ = self.complexes[i]

            reactions += mech_b.update_reactions(regulator, self.dna_to_bind, component = self, \
                                                                    part_id = self.name+"_"+regulator.name,protein=self.protein)
            
            reactions += mech_tx.update_reactions(dna = complex_, component = self, part_id = self.name+"_"+regulator.name, \
                                            transcript = self.transcript, protein = self.get_protein_for_expression())

        return reactions


class ActivatablePromoter(Promoter):
    """
    A class for a promoter which can be activated by a single species, modelled as a positive hill function
    """
    def __init__(self, name, activator, transcript = None, leak = False, **keywords):
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, **keywords)

        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        self.activator = self.set_species(activator) 

        self.leak = leak #toggles whether or not there is a leak reaction
        
        #Non-default Mechanisms are added to the Component with .add_mechanism
        self.add_mechanism(PositiveHillTranscription(), "transcription", overwrite = True)
        

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.get_mechanism("transcription")
        
        species = [] #A list of species must be returned
        species += mech_tx.update_species(dna = self.dna_to_bind, transcript = self.transcript, regulator = self.activator, 
            part_id = self.name+"_"+self.activator.name, leak = self.leak, component = self, protein = self.get_protein_for_expression())
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.get_mechanism("transcription")
        
        reactions = [] #a list of reactions must be returned


        reactions += mech_tx.update_reactions(dna = self.dna_to_bind, transcript = self.transcript, regulator = self.activator, 
                                             component = self, part_id = self.name+"_"+self.activator.name, leak = self.leak, protein = self.get_protein_for_expression())

        
        return reactions

class RepressiblePromoter(Promoter):
    """
    A class for a promoter which can be repressed by a single species, modelled as a negative hill function
    """
    def __init__(self, name, repressor, transcript = None, leak = False, **keywords):
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, **keywords)

        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        self.repressor = self.set_species(repressor) 

        self.leak = leak #toggles whether or not there is a leak reaction
        
        #Mechanisms are inherited from the Mixture unless set specifically in self.default_mechanisms.
        self.add_mechanism(NegativeHillTranscription(), "transcription", overwrite = True)
        
        

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.get_mechanism("transcription")
        
        species = [] #A list of species must be returned
        species += mech_tx.update_species(dna = self.dna_to_bind, transcript = self.transcript, regulator = self.repressor, component = self, 
            part_id = self.name+"_"+self.repressor.name, leak = self.leak, protein = self.get_protein_for_expression(), **keywords)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.get_mechanism("transcription")
        
        reactions = [] #a list of reactions must be returned

        reactions += mech_tx.update_reactions(dna = self.dna_to_bind, transcript = self.transcript, regulator = self.repressor, 
                                             component = self, part_id = self.name+"_"+self.repressor.name, leak = self.leak, protein = self.get_protein_for_expression(), **keywords)
        return reactions


class CombinatorialPromoter(Promoter):
    def __init__(self, name, regulators, leak = False, assembly = None,
                 transcript = None, length = 0, mechanisms = None,
                 parameters = None,protein=None,tx_capable_list = None,cooperativity = None, **keywords):
        """
        A combinatorial promoter is something where binding multiple regulators result in
        qualitatively different transcription behaviour. For example, maybe it's an AND
        gate promoter where it only transcribes if two regulators are bound, but not if 
        either one is bound.
        
        =============
        inputs
        =============
        name: the name of the promoter
        regulators: a list of strings or species indicating all the possible regualtors that can bind
        
        leak: if true, then a promoter with nothing bound will transcribe
        
        assembly: a DNA_assembly object that contains this promoter
        
        transcript: the transcript that this promoter makes
        
        length: the length in nt? I don't think this is used for anything at the moment
        
        mechanisms: additional mechanisms. formatted with {"mechanism_type":mechanismObject(),...}
        
        parameters: promoter-specific parameters. Formatted as {("identifier1","identifier2"):value,...}
        
        tx_capable_list: list of which combination of regulators bound will lead to transcription. 
                        formatted as [["regulator1","regulator2"],["regulator1"],...] regulators
                        can be strings or Species
        cooperativity: a dictionary of cooperativity values. For example, {"regulator":2,"regulator2":1,....}
        """

        Promoter.__init__(self, name = name, assembly = assembly,
                          transcript = transcript, length = length,
                          mechanisms = mechanisms, parameters = parameters, protein=protein,
                          **keywords)

        if not isinstance(regulators, list):
            #you could give one string as a regulator
            regulators = [regulators]
        self.cooperativity = cooperativity
        self.regulators = []
        for regulator in regulators:
            self.regulators += [self.set_species(regulator, material_type = "protein")]
            
        #after we've sanitized the inputs, then sort
        self.regulators = sorted(self.regulators)
        #now let's work out the tx_capable_list
        if tx_capable_list is None:
            #if nothing is passed, that means everything transcribes
            allcomb = []
            for r in range(1,len(self.regulators)+1):
                #make all combinations of regulators
                allcomb += [set(a) for a in it.combinations([a.name for a in self.regulators],r)]
            self.tx_capable_list = allcomb
        elif isinstance(tx_capable_list, list):
            #if the user passed a list then the user knows what they want
            newlist = []
            #this part converts any species in the tx_capable_list into a string
            for element in tx_capable_list:
                sublist = []
                for specie in element:
                    if(isinstance(specie,Species)):
                        sublist += [specie.name]
                    else:
                        sublist += [specie]
                newlist+=[sublist]
            self.tx_capable_list = [set(a) for a in newlist]

        self.leak = leak
        self.complex_combinations = {}
        self.tx_capable_complexes = []
        self.leak_complexes = []
        self.add_mechanism(Combinatorial_Cooperative_Binding(), "binding", overwrite = True)
        
        

    def update_species(self):

        mech_tx = self.get_mechanism("transcription")
        mech_b = self.get_mechanism('binding')
        #set the tx_capable_complexes to nothing because we havent updated species yet!
        self.tx_capable_complexes = []
        self.leak_complexes = []
        species = [self.dna_to_bind]+self.regulators
        self.complexes = []
        bound_species = mech_b.update_species(self.regulators,self.dna_to_bind,\
                        component = self,part_id = self.name,cooperativity=self.cooperativity,protein=self.protein)
        #above is all the species with DNA bound to regulators. Now, we need to extract only the ones which
        #are transcribable

        if self.leak is not False:
            #this part takes care of the promoter not bound to anything
            species += mech_tx.update_species(dna = self.dna_to_bind, transcript = self.transcript, protein = self.get_protein_for_expression(), component = self, part_id = self.name+"_leak")
            #self.leak_complexes += []

        for bound_complex in bound_species: 
            species_inside = []
            for regulator in self.regulators:
                if(regulator in bound_complex):
                    species_inside += [regulator.name] 
            if(set(species_inside) in [set(a) for a in self.tx_capable_list]):
                #only the transcribable complexes in tx_capable_list get transcription reactions
                tx_capable_species = mech_tx.update_species(dna = bound_complex, transcript = self.transcript, \
                                                protein = self.get_protein_for_expression(), component = self, part_id = self.name)
                species +=tx_capable_species
                self.tx_capable_complexes +=[bound_complex]
            else:
                #in this case there's a combination of regulators which does not feature in tx_capable_list
                #this means: 
                # 1) this complex does nothing so don't add it
                # 2) we said we wanted leak, so then you should add this, but with the "_leak" parameters
                #                                                         (that happens in update_reactions)
                if(self.leak is not False):
                    leak_species = mech_tx.update_species(dna = bound_complex, transcript = self.transcript, \
                                        protein = self.get_protein_for_expression(), component = self, part_id = self.name+"_leak")
                    species += leak_species
                    self.leak_complexes += [bound_complex]
        species+=bound_species  
        return species

    def update_reactions(self):
        
        reactions = []
        mech_tx = self.get_mechanism("transcription")
        mech_b = self.get_mechanism('binding')

        if self.leak is not False:
            #once again the DNA not bound to anything gets special treatment
            reactions += mech_tx.update_reactions(dna = self.dna_to_bind, component = self, part_id = self.name+"_leak", \
                                                            transcript = self.transcript, protein = self.protein)
        #the binding reactions happen no matter what
        reactions += mech_b.update_reactions(self.regulators,self.dna_to_bind,component = self,\
                                                        part_id = self.name,cooperativity=self.cooperativity, protein = self.protein)
        if((self.tx_capable_complexes is None) or self.tx_capable_complexes == []):
            species = self.update_species()
            #this could mean we haven't run update_species() yet
            
            if(self.tx_capable_complexes == []):
                if(self.leak_complexes is None or self.leak_complexes == []):
                    #if it's still zero after running update_species then we could be in trouble
                    warn("nothing can transcribe from combinatorial promoter {}".format(self.name))
            
        if(len(self.tx_capable_complexes)>0):
            for specie in self.tx_capable_complexes:
                tx_partid = self.name
                for part in specie.species_set:
                    #construct the name of the promoter with regulators bound
                    if part.material_type == "dna":
                        #the DNA doesn't matter
                        pass
                    else:
                        #put in the regulators!
                        tx_partid += "_"+part.name
                if(tx_partid[0]=="_"):
                    #this will only happen if the name of the dna is ""
                    tx_partid = tx_partid[1:]
                #if it's bound to RNAP then it transcribes, right?
                tx_partid = tx_partid+"_RNAP"
                reactions += mech_tx.update_reactions(dna = specie, component = self, part_id = tx_partid, \
                                            transcript = self.transcript, protein = self.get_protein_for_expression())
        if(len(self.leak_complexes)>0):
            for specie in self.leak_complexes:
                #in this case every reaction uses the "promoter_leak" partid
                leak_partid = self.name+"_leak"
                reactions += mech_tx.update_reactions(dna = specie, component = self, part_id = leak_partid, \
                                            transcript = self.transcript, protein = self.get_protein_for_expression())

        return reactions
