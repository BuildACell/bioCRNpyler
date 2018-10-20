# import Component
import mechanism
import chemical_reaction_network as crn
from .component import Component, RNA
from .mechanism import Reversible_Bimolecular_Binding

class guideRNA(RNA):
    def __init__(self, guide_name, dCas9 = "dCas9", **keywords):
        self.gRNA = crn.specie(guide_name, type = "rna")

        if isinstance(dCas9, crn.specie):
            self.dCas = dCas9
        elif isinstance(dCas9, str):
            self.dCas = crn.specie(dCas9, type = "protein")
        else:
            raise ValueError("dCas9 parameter must be a chemical_reaction_network.specie or a string")


        self.default_mechanisms = {
            "dCas9_binding" : mechanism.Reversible_Bimolecular_Binding()
        }
        component.RNA.__init__(self, name = guide_name, **keywords)

    def update_species(self):
        species = [self.gRNA, self.dCas]
        species += self.mechanisms['dCas9_binding'].update_species(self.gRNA, self.dCas)
        return species

    def update_reactions(self):
        ku = self.get_parameter("ku", part_id=self.gRNA.name, mechanism=self.mechanisms['dCas9_binding'])
        kb = self.get_parameter("kb", part_id=self.gRNA.name, mechanism=self.mechanisms['dCas9_binding'])
        rxns = self.mechanisms['dCas9_binding'].update_reactions(self.gRNA, self.dCas, kb = kb, ku = ku)
        return rxns