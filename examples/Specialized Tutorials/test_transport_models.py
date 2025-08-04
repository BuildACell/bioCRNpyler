from biocrnpyler import *
alphaHL_monomer = IntegralMembraneProtein('alphaHL_monomer',
                                          product='alphaHL',
                                          size=7)
alphaHL_channel = MembraneChannel(alphaHL_monomer.product,
                                  substrate='ATP')
compartment_internal = alphaHL_monomer.get_species().compartment

# ActivatedPromoter
activator = Species("T7RNAP",
                    material_type = "small_molecule",
                    compartment=compartment_internal)
hill_parameters = {"k":1.0, "n":4, "K":20, "kleak":0.0001}
P_activatable = ActivatablePromoter("P_activtable",
                                    activator = activator,
                                    leak = False,
                                    parameters = hill_parameters)
activatable_assembly = DNAassembly("activatable_assembly",
                                   promoter = P_activatable,
                                   rbs = "rbs",
                                   initial_concentration = 1*10**-3,
                                   protein= alphaHL_monomer.membrane_protein)
#Mechanisms
mech_int = Membrane_Protein_Integration() 
mech_tra = Simple_Transport()
all_mechanisms = {
    mech_tra.mechanism_type:mech_tra,
    mech_int.mechanism_type:mech_int
}
E = EnergyTxTlExtract(components=[activatable_assembly,
                                  alphaHL_monomer,
                                  alphaHL_channel,
                                  ],
                      mechanisms = all_mechanisms,
                      parameter_file = "all_parameters.txt")
CRN = E.compile_crn(compartment=compartment_internal)
print(CRN.species)
CRN.write_sbml_file("test_transport_models.xml")

# Load SBML using libsbml
import libsbml
model = libsbml.readSBML("test_transport_models.xml").getModel()
# assert 2 compartments, one Internal and one External
assert model.getNumCompartments() == 2
for compartment in model.getListOfCompartments():
    assert compartment.getName() in ["Internal", "External"]
    assert compartment.getId() in ["Internal", "External"]

# assert that there are 18 species in the model
assert model.getNumSpecies() == 18

