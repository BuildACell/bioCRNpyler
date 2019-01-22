# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from .component import Protein, Complex
from .mechanism import Transcription_MM, Translation_MM, Degredation_mRNA_MM
from .mixture import Mixture


class BasicExtract(Mixture):
    def __init__(self, name="", mechanisms={}, components=[], parameters={},
                 rnap_name="RNAP", ribosome_name="Ribo", rnaase_name="RNAase", **kwargs):

        self.rnap = Protein(name=rnap_name)
        self.ribosome = Complex(name=ribosome_name)
        self.RNAase = Protein(name=rnaase_name)

        mech_tx = Transcription_MM(rnap=self.rnap.get_specie())
        mech_tl = Translation_MM(ribosome=self.ribosome.get_specie())
        mech_rna_deg = Degredation_mRNA_MM(nuclease=self.RNAase.get_specie())

        default_mechanisms = {
            mech_tx.type: mech_tx,
            mech_tl.type: mech_tl,
            mech_rna_deg.type: mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.RNAase]

        Mixture.__init__(self, name=name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, components=components+default_components, parameters=parameters, **kwargs)

       
class txtl(Mixture):
    '''
    TODO : Implement this
    '''
    def __init__(self, Mixture):
        self.Mixture = Mixture
        self.name = Mixture.name

    def extract(self, name, parameters = {}):
        return

    def buffer(self, name = "", components = [], parameters = {}, **kwargs):
        self.name = name
        return 
   
    def new_tube(self, name = ""):
        return

    def add_dna(self, dna):
        return

    def combine_tubes(self):
        crn = self.compile_crn()
        return crn

