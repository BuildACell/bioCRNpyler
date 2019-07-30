# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from .component import Protein, Complex
from .mechanism import Transcription_MM, Translation_MM, Degredation_mRNA_MM
from .mixture import Mixture
from .chemical_reaction_network import Specie
       

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
        if not init:
            warn('Initial concentrations for extract species will all be set to zero.')
        if isinstance(rnap, Specie):
            self.rnap = rnap
        elif isinstance(rnap, str):
            self.rnap = Protein(name=rnap)
        else:
            raise ValueError("rnap parameter must be a str or chemical_reaction_network.specie")

        if isinstance(ribosome, Specie):
            self.ribosome = ribosome
        if isinstance(ribosome, str):
            self.ribosome = Protein(name=ribosome)
        else:
            raise ValueError("rnap parameter must be a str or chemical_reaction_network.specie")

        if isinstance(rnaase, Specie):
            self.rnaase = rnaase
        elif isinstance(rnaase, str):
            self.RNAase = Protein(name=rnaase)
        else:
            raise ValueError("rnaase parameter must be a str or chemical_reaction_network.specie")

        self.rnap.get_specie().initial_concentration = init["protein_RNAP"]
        self.RNAase.get_specie().initial_concentration = init["protein_RNAase"]
        self.ribosome.get_specie().initial_concentration = init["protein_Ribo"]
        mech_tx = Transcription_MM(rnap = self.rnap.get_specie())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_specie())
        mech_rna_deg = Degredation_mRNA_MM(nuclease = self.RNAase.get_specie())


        default_mechanisms = {
            mech_tx.type: mech_tx,
            mech_tl.type: mech_tl,
            mech_rna_deg.type: mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.RNAase]
        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, 
                        components=components+default_components, **kwargs)

