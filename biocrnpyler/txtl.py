from warnings import warn
import component as comp
import chemical_reaction_network as crn
import mechanism
import mixture

class BasicExtract(mixture.Mixture):
    def __init__(self, name="", mechanisms={}, components=[], parameters={},
                 rnap_name = "RNAP", ribosome_name = "Ribo", rnaase_name = "RNAase", **kwargs):

        self.rnap = comp.Protein(name=rnap_name)
        self.ribosome = comp.Complex(name=ribosome_name)
        self.RNAase = comp.Protein(name=rnaase_name)

        mech_tx = mechanism.Transcription_MM(rnap = self.rnap.get_specie())
        mech_tl = mechanism.Translation_MM(ribosome=self.ribosome.get_specie())
        mech_rna_deg = mechanism.Degredation_mRNA_MM(nuclease = self.RNAase.get_specie())

        default_mechanisms = {
            mech_tx.type:mech_tx,
            mech_tl.type:mech_tl,
            mech_rna_deg.type:mech_rna_deg
        }

        default_components = [self.rnap, self.ribosome, self.RNAase]

        mixture.Mixture.__init__(self, name = name, default_mechanisms=default_mechanisms, mechanisms=mechanisms, components=components+default_components, parameters=parameters, **kwargs)