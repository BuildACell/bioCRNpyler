import component
import mechanism
import chemical_reaction_network as crn

class guideRNA(component.RNA):
    def __init__(self, guide_name, dCas9 = "dCas9", **keywords):


        if isinstance(dCas9, crn.specie):
            self.dCas = dCas9
        elif isinstance(dCas9, str):
            self.dCas = crn.specie(dCas9, type = "protein")
        elif isinstance(dCas9, component.Component) and dCas9.get_specie()!= None:
            self.dCas = dCas9.get_specie()
        else:
            raise ValueError("dCas9 parameter must be a chemical_reaction_network.specie, Component with get_specie(), or a string")


        self.default_mechanisms = {
            "dCas9_binding" : mechanism.Reversible_Bimolecular_Binding(name = "dCas9_binding", type = "bimolecular binding")
        }

        component.RNA.__init__(self, name = guide_name, **keywords)
        self.gRNA = self.get_specie()

    def get_dCasComplex(self):
        binding_species = self.mechanisms['dCas9_binding'].update_species(self.gRNA, self.dCas)
        if len(binding_species) > 1:
            raise ValueError('dCas9_binding mechanisms '+self.mechanisms['dCas9_binding'].name+" returned multiple complexes. Unclear which is active." )
        else:
            return binding_species[0]

    def update_species(self):
        species = [self.gRNA, self.dCas]
        species += self.mechanisms['dCas9_binding'].update_species(self.gRNA, self.dCas)
        return species

    def update_reactions(self):
        ku = self.get_parameter("ku", part_id=self.gRNA.name, mechanism=self.mechanisms['dCas9_binding'])
        kb = self.get_parameter("kb", part_id=self.gRNA.name, mechanism=self.mechanisms['dCas9_binding'])
        rxns = self.mechanisms['dCas9_binding'].update_reactions(self.gRNA, self.dCas, kb = kb, ku = ku)
        return rxns