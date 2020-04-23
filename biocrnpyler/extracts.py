# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from .component import Protein, ChemicalComplex
from .mechanism import Transcription_MM, Translation_MM, Degredation_mRNA_MM, EmptyMechanism, OneStepGeneExpression
from .mixture import Mixture
from .chemical_reaction_network import Species
       
#A Model for Gene Expression without any Machinery (eg Ribosomes, Polymerases, etc.)
class ExpressionExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], **kwargs):

        dummy_translation = EmptyMechanism(name = "dummy_translation", mechanism_type = "translation")
        mech_expression = OneStepGeneExpression()

        default_mechanisms = {
            "transcription": mech_expression,
            "translation": dummy_translation
        }

        default_components = []
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

#A Model for Transcription and Translation in Cell Extract with Ribosomes, Polymerases, and Endonucleases.
#This model does not include any energy
class TxTlExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                 rnap = "RNAP", ribosome = "Ribo", rnaase = "RNAase", **kwargs):
        
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)

        init = kwargs.get('init')
        if init:
            self.rnap.get_species().initial_concentration = init[rep(rnap)]
            self.rnaase.get_species().initial_concentration = init[repr(rnaase)]
            self.ribosome.get_species().initial_concentration = init[repr(ribosome)]

        mech_tx = Transcription_MM(rnap = self.rnap.get_species())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_species())
        mech_rna_deg = Degredation_mRNA_MM(nuclease = self.rnaase.get_species()) 


        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.rnaase]
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

# Below are unimplemented classes...
# Outline for BasicBuffer
# TODO
class BasicBuffer(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                atp = "ATP", ntp = "NTP", aa  = "AA", **kwargs):
    
        if isinstance(ntp, Species):
            self.ntp = ntp
        if isinstance(ntp, str):
            self.ntp = Protein(name=ntp)
        else:
            raise ValueError("ntp argument must be a str or chemical_reaction_network.specie")

        if isinstance(aa, Species):
            self.aa = aa
        elif isinstance(aa, str):
            self.aa = Protein(name=aa)
        else:
            raise ValueError("aa argument must be a str or chemical_reaction_network.specie")

        self.atp.get_species().initial_concentration = init["protein_RNAP"]
        self.ntp.get_species().initial_concentration = init["protein_RNAase"]
        self.aa.get_species().initial_concentration = init["protein_Ribo"]
        # mech_tx = Transcription_MM(rnap = self.atp.get_specie()) # TODO: Implement energy mechanisms here
        # mech_tl = Translation_MM(ribosome = self.ntp.get_specie())
        # mech_rna_deg = Degredation_mRNA_MM(nuclease = self.aa.get_specie())


        # default_mechanisms = {
        #     mech_tx.type: mech_tx,
        #     mech_tl.type: mech_tl,
        #     mech_rna_deg.type: mech_rna_deg
        # }

        default_components = [self.atp, self.ntp, self.aa]
        # Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        # components=components+default_components, **kwargs)
        raise NotImplementedError

class Extract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], 
                rnap = "RNAP", ribosome = "Ribo", rnaase = "RNAase", **kwargs):
        raise NotImplementedError

class CustomExtract(Mixture):
    def __init__(self, name = "", **kwargs):
        raise NotImplementedError


class EnergyBuffer(Mixture):
    def __init__(self, name = "", **kwargs):
        raise NotImplementedError
                