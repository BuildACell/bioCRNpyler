# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from .component import Protein, ChemicalComplex
from .mechanism import Transcription_MM, Translation_MM, Degredation_mRNA_MM
from .mixture import Mixture
from .chemical_reaction_network import Species
       

class BasicExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                 rnap = "RNAP", ribosome = "Ribo", rnaase = "RNAase", **kwargs):
        init = kwargs.get('init')
        parameter_warnings = kwargs.get('parameter_warnings')
        if parameter_warnings:
            warn('Parameter warnings have been set True. Verbose warnings regarding parameter files will be displayed.')
        else:
            parameter_warnings = False
            kwargs['parameter_warnings'] = parameter_warnings
        if not init and parameter_warnings:
            warn('Initial concentrations for extract species will all be set to zero.')
        if isinstance(rnap, Species):
            self.rnap = rnap
        elif isinstance(rnap, str):
            self.rnap = Protein(name=rnap)
        else:
            raise ValueError("rnap argument must be a str or chemical_reaction_network.species")

        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        if isinstance(ribosome, str):
            self.ribosome = Protein(name=ribosome)
        else:
            raise ValueError("rnap argument must be a str or chemical_reaction_network.species")

        if isinstance(rnaase, Species):
            self.rnaase = rnaase
        elif isinstance(rnaase, str):
            self.RNAase = Protein(name=rnaase)
        else:
            raise ValueError("rnaase argument must be a str or chemical_reaction_network.species")

        if init:
            self.rnap.get_species().initial_concentration = init["protein_RNAP"]
            self.RNAase.get_species().initial_concentration = init["protein_RNAase"]
            self.ribosome.get_specie().initial_concentration = init["protein_Ribo"]
        mech_tx = Transcription_MM(rnap = self.rnap.get_species())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_species())
        mech_rna_deg = Degredation_mRNA_MM(nuclease = self.RNAase.get_species())


        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_rna_deg.mechanism_type: mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.RNAase]
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

# Outline for BasicBuffer
class BasicBuffer(Mixture):
    def __init__(self, name="", mechanisms={}, components=[],
                atp = "ATP", ntp = "NTP", aa  = "AA", **kwargs):
        init = kwargs.get('init')
        parameter_warnings = kwargs.get('parameter_warnings')
        if parameter_warnings:
            warn('Parameter warnings have been set True. Verbose warnings regarding parameter files will be displayed.')
        else:
            parameter_warnings = False
            kwargs['parameter_warnings'] = parameter_warnings
        if not init:
            warn('Initial concentrations for extract species will all be set to zero.')
        if isinstance(atp, Species):
            self.atp = atp
        elif isinstance(atp, str):
            self.atp = Protein(name=atp)
        else:
            raise ValueError("atp argument must be a str or chemical_reaction_network.specie")

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

# You can add other extract models here based on the code above.