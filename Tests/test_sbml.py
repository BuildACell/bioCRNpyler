import pytest
from biocrnpyler.sbmlutil import *
from biocrnpyler.species import Species, Complex
from biocrnpyler.propensities import MassAction, HillPositive, HillNegative, ProportionalHillPositive, ProportionalHillNegative 
from biocrnpyler.parameter import ParameterEntry, ParameterKey
from biocrnpyler.chemical_reaction_network import Reaction, ChemicalReactionNetwork


def test_create_sbml_model():
    #tests creating an SBML model

    document, model = create_sbml_model()


def test_add_all_species():
    #Tests add_species via add_all_species

    document, model = create_sbml_model()

    #Some Species
    S1, S2, S3, S4 = Species("S1"), Species("S2"), Species("S3"), Species("S4")

    #These two Complexes push the naming convention to its limits
    C1 = Complex([Complex([S1, S2, S3]), S4])
    C2 = Complex([Complex([S1, S2]), S3, S4])

    species = [S1, S2, S3, S4, C1, C2]
    add_all_species(model, species)

    assert len(model.getListOfSpecies()) == len(species)
    assert test_validate_sbml(document) == 0



def test_add_reaction():
    #tests adding reactions of all types, with different stochiometries

    
    #create species
    S1, S2, S3 = Species("S1"), Species("S2"), Species("S3")
    species = [S1, S2, S3]
    #create some parameters
    key1 = ParameterKey(name = "k", mechanism = "m", part_id = "pid")
    k1 = ParameterEntry("k", 1.11, key1)
    k2 = ParameterEntry("k", 2.22)

    stochiometries = [([], [S1]), ([S1], []), ([S1], [S2]), (2*[S1], [S2, S3])]
    propensities = [MassAction(k_forward = 1.), MassAction(k_forward = 1, k_reverse = .1), MassAction(k_forward = k1), MassAction(k_forward = k1, k_reverse = k2),
                    HillPositive(k = 1, n = 2., K = 3., s1 = S1), HillPositive(k = k1, n = 2., K = k2, s1 = S1),
                    HillNegative(k = 1, n = 2., K = 3., s1 = S1), HillNegative(k = k1, n = k2, K = 3., s1 = S1),
                    ProportionalHillPositive(k = 1, n = 2, K = 3., s1 = S1, d = S2), ProportionalHillPositive(k = k1, n = 2, K = k2, s1 = S1, d = S2),
                    ProportionalHillNegative(k = 1, n = 2, K = 3., s1 = S1, d = S2), ProportionalHillNegative(k = k1, n = 2, K = k2, s1 = S1, d = S2)
                    ]

    rxn_num = 0
    for prop in propensities: #Cycle through different propensity types
        for inputs, outputs in stochiometries: #cycle through different stochiometries
            for stochastic in [True, False]: #Toggle Stochastic
                for for_bioscrape in [True, False]: #Toggle for_bioscrape
                    document, model = create_sbml_model()
                    add_all_species(model, species)
                    rxn = Reaction(inputs, outputs, propensity_type = prop) #create a reaction
                    args = {"stochastic":stochastic, "for_bioscrape":for_bioscrape} #args for printing
                    try:
                        sbml_reaction = add_reaction(model, rxn, f"r{rxn_num}", stochastic = stochastic, for_bioscrape = for_bioscrape) #add reaction to SBML Model
                        rxn_num += 1 #increment reaction id
                        crn_reactants = [str(s) for s in inputs]
                        crn_products = [str(s) for s in outputs]
                        sbml_products = [p.getSpecies() for p in sbml_reaction.getListOfProducts()]
                        sbml_reactants = [p.getSpecies() for p in sbml_reaction.getListOfReactants()]

                        #test that the reaction has the right inputs and outputs
                        assert set(crn_reactants) == set(sbml_reactants)
                        assert set(crn_products) == set(sbml_products)

                        if len(crn_reactants) > 0:
                            assert all(crn_reactants.count(s.getSpecies()) == int(s.getStoichiometry()) for s in sbml_reaction.getListOfReactants())

                        if len(crn_products) > 0:
                            assert all(crn_products.count(s.getSpecies()) == int(s.getStoichiometry()) for s in sbml_reaction.getListOfProducts())

                        # sbml_reaction.setReversible(False)
                        # print(sbml_reaction.isSetReversible())
                        # assert not sbml_reaction.isSetReversible()

                        #TODO is there a smart way to test that the rate formula is correct?

                        #test annotations
                        if for_bioscrape:
                            sbml_annotation = sbml_reaction.getAnnotation().toXMLString()
                            assert f"type={prop.name}" in sbml_annotation

                            #Test that the sbml annotation has keys for all species and parameters
                            for k in prop.propensity_dict["parameters"]:
                                #convert k_reverse and k_forward to just k
                                k = k.replace("_reverse", "").replace("_forward", "")
                                assert f"{k}=" in sbml_annotation

                            for s in prop.propensity_dict["species"]:
                                assert f"{s}=" in sbml_annotation
                        else:
                            sbml_annotation = sbml_reaction.getAnnotation()
                            assert sbml_annotation is None

                        #Make sure this whole thing can be written by libsml, for good measure!
                        sbml_string = libsbml.writeSBMLToString(document)
                        assert test_validate_sbml(document) == 0
                    except Exception as e: #collect errors to display with args
                        error_txt = f"Unexpected Error: in sbmlutil.add_reaction {rxn} with args {args}. \n {str(e)}."
                        raise Exception(error_txt)

def test_generate_sbml_model():

    #Test a non-reversible reaction
    s1 = Species("S1")
    s2 =  Species("S2")
    rx0 = Reaction.from_massaction(inputs=[s1], outputs=[s2], k_forward=0.1)
    crn = ChemicalReactionNetwork(species = [s1, s2], reactions = [rx0])

    # generate an sbml model
    document, model = crn.generate_sbml_model()
    # all species from the CRN are accounted for
    assert len(model.getListOfSpecies()) == len(crn.species)
    # all reactions from the CRN are accounted for
    assert len(model.getListOfReactions()) == len(crn.reactions)

    #reversible needs to be off!
    # assert not model.getListOfReactions()[0].isSetReversible()

    # test a reversible reaction
    rx1 = Reaction.from_massaction(inputs=[s1], outputs=[s2], k_forward=0.1, k_reverse=0.1)
    rxn_list = [rx1]
    crn = ChemicalReactionNetwork(species=[s1, s2], reactions=rxn_list)

    # generate an sbml model
    document, model = crn.generate_sbml_model()
    # all species from the CRN are accounted for
    assert len(model.getListOfSpecies())==len(crn.species)
    # all reactions from the CRN are accounted for
    assert len(model.getListOfReactions())== 2
    # although  sbml represents a reverisble reaction with reversible flag
    # BioCRNpyler always creates two reactions, because this is correct
    # for stochastic simulation with SBML.
    sbml_rxn = model.getListOfReactions()
    # assert not sbml_rxn[0].isSetReversible()
    # assert not sbml_rxn[1].isSetReversible()

    #Test propagation for for_bioscrape keyword
    # generate an sbml model with for_bioscrape = True
    document, model = crn.generate_sbml_model(for_bioscrape = True)
    for r in model.getListOfReactions():
        assert r.getAnnotation() is not None
    assert test_validate_sbml(document) == 0
    # generate an sbml model with for_bioscrape = False
    document, model = crn.generate_sbml_model(for_bioscrape = False)
    for r in model.getListOfReactions():
        assert r.getAnnotation() is None
    assert test_validate_sbml(document) == 0


def test_generate_sbml_model_parameter_names():

    s1 = Species("S1")
    s2 =  Species("S2")
    
    #The correct parameter formating should be: "name_partid_mechanism"
    key0 = ParameterKey(mechanism = "m", part_id = "p", name = "n")
    k0 = ParameterEntry("n", 1.0, parameter_key = key0)

    #Tests 3 Parameter entries with seemingly redundant keys.
    #In SBML parameter id, None will be kept blank, resulting in
    #different numbers of underscores between the V's
    k1 = ParameterEntry("v", 1.0, parameter_key = None)

    key2 = ParameterKey(mechanism = None, part_id = "v", name = "v")
    k2 = ParameterEntry("v", 2.0, parameter_key = key2)

    key3 = ParameterKey(mechanism = "v", part_id = None, name = "v")
    k3 = ParameterEntry("v", 2.0, parameter_key = key3)

    key4 = ParameterKey(mechanism = "v", part_id = "v", name = "v")
    k4 = ParameterEntry("v", 2.0, parameter_key = key4)

    rx0 = Reaction.from_massaction(inputs = [], outputs = [s1], k_forward = k0, k_reverse = k1) #Throw a duplicate key for good measure!
    rx1 = Reaction.from_massaction(inputs=[s1], outputs=[s2], k_forward=k1, k_reverse=k2)
    rx2 = Reaction.from_massaction(inputs=2*[s1], outputs=2*[s2], k_forward=k3, k_reverse=k4)
    rxn_list = [rx0, rx1, rx2]
    crn = ChemicalReactionNetwork(species=[s1, s2], reactions=rxn_list)

    document, model = crn.generate_sbml_model()
    assert test_validate_sbml(document) == 0
    correct_ids = set(["v_v_", "v__v", "v_v_v", "v__", "n_p_m"])
    ids = set([p.getId() for p in model.getListOfParameters()])
    assert ids == correct_ids
