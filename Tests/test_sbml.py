import pytest
import warnings
from unittest import TestCase
from biocrnpyler.sbmlutil import *
from biocrnpyler.species import Species, Complex
from biocrnpyler.compartments import Compartment
from biocrnpyler.propensities import MassAction, HillPositive, HillNegative, ProportionalHillPositive, ProportionalHillNegative, GeneralPropensity
from biocrnpyler.parameter import ParameterEntry, ParameterKey
from biocrnpyler.chemical_reaction_network import Reaction, ChemicalReactionNetwork


class TestSBML(TestCase):
    def test_create_sbml_model(self):
        # tests creating an SBML model

        document, model = create_sbml_model()

    def test_add_all_species(self):
        # Tests add_species via add_all_species

        document, model = create_sbml_model()

        # Some Species
        S1, S2, S3, S4 = Species("S1"), Species(
            "S2"), Species("S3"), Species("S4")

        # These two Complexes push the naming convention to its limits
        C1 = Complex([Complex([S1, S2, S3]), S4])
        C2 = Complex([Complex([S1, S2]), S3, S4])

        species = [S1, S2, S3, S4, C1, C2]

        # initial conditions
        init_cond_dict = {S1: 1, S2: 2, S3: 3, S4: 4}

        add_all_species(model, species, init_cond_dict)
        # Test get species by name:
        self.assertEqual(str(getSpeciesByName(model, model.getSpecies(0).getName()).getName()), "S1")
        with self.assertRaisesRegex(ValueError, '"name" must be a string. Received 24.'):
            getSpeciesByName(model, 24)

        self.assertEqual(len(model.getListOfSpecies()), len(species))
        self.assertEqual(validate_sbml(document), 0)

    def test_add_all_compartments(self):
        document, model = create_sbml_model()

        # Some Species
        S1 = Species("S1", compartment="S1_compartment")
        S2 = Species("S2", compartment="S2_compartment")
        S3_compartment = Compartment("S3_compartment", size=2e-4)
        S3 = Species("S3", compartment=S3_compartment)
        S4_compartment = Compartment(
            "S4_compartment", size=4e-4, spatial_dimensions=2, unit = "litre")
        S4 = Species("S4", compartment=S4_compartment, material_type="protein")

        species = [S1, S2, S3, S4]
        list_of_compartments = [S1.compartment,
                                S2.compartment, S3_compartment, S4_compartment]
        list_of_sizes = [1e-6, 1e-6, 2e-4, 4e-4]
        list_of_dimensions = [3, 3, 3, 2]
        add_all_compartments(model, list_of_compartments)
        # initial conditions
        init_cond_dict = {S1: 1, S2: 2, S3: 3, S4: 4}
        add_all_species(model, species, init_cond_dict)

        self.assertEqual(len(model.getListOfCompartments()),
                         len(list_of_compartments))
        for compartment, size, dim in zip(model.getListOfCompartments(), list_of_sizes, list_of_dimensions):
            self.assertEqual(compartment.getSize(), size)
            self.assertEqual(compartment.getSpatialDimensions(), dim)
        self.assertEqual(len(model.getListOfSpecies()), len(species))
        self.assertEqual(validate_sbml(document), 0)

    def test_add_reaction(self):
        """
        Test adding reaction to SBML for combinatorially created list of models
        """
        # create species
        S1, S2, S3 = Species("S1"), Species("S2"), Species("S3")
        species = [S1, S2, S3]
        # create some parameters
        key1 = ParameterKey(name="k", mechanism="m", part_id="pid")
        k1 = ParameterEntry("k1", 1.11, key1)
        k2 = ParameterEntry("k2", 2.22)
        stochiometries = [([], [S1]), ([S1], []),
                          ([S1], [S2]), (2*[S1], [S2, S3])]
        propensities = [MassAction(k_forward=1.), MassAction(k_forward=1, k_reverse=.1), MassAction(k_forward=k1), MassAction(k_forward=k1, k_reverse=k2),
                        HillPositive(k=1, n=2., K=3., s1=S1), HillPositive(
                            k=k1, n=2., K=k2, s1=S1),
                        HillNegative(k=1, n=2., K=3., s1=S1), HillNegative(
                            k=k1, n=k2, K=3., s1=S1),
                        ProportionalHillPositive(k=1, n=2, K=3., s1=S1, d=S2), ProportionalHillPositive(
                            k=k1, n=2, K=k2, s1=S1, d=S2),
                        ProportionalHillNegative(k=1, n=2, K=3., s1=S1, d=S2), ProportionalHillNegative(
                            k=k1, n=2, K=k2, s1=S1, d=S2),
                        GeneralPropensity('k1*2 - k2/S1^2', propensity_species=[S1], propensity_parameters=[k1, k2]), GeneralPropensity(
                            'S1^2 + S2^2 + S3^2', propensity_species=[S1, S2, S3], propensity_parameters=[])
                        ]

        for prop in propensities:  # Cycle through different propensity types
            for inputs, outputs in stochiometries:  # cycle through different stochiometries
                for stochastic in [True, False]:  # Toggle Stochastic
                    rxn_num = 0
                    model_id = f"{prop.name}_model_with_stochastic_{stochastic}"
                    document, model = create_sbml_model(model_id=model_id)
                    add_all_species(model, species, {})
                    # create a reaction
                    rxn = Reaction(inputs, outputs, propensity_type=prop)
                    try:
                        # add reaction to SBML Model
                        sbml_reaction = add_reaction(
                            model, rxn, f"r{rxn_num}", stochastic=stochastic)
                        rxn_num += 1  # increment reaction id
                        crn_reactants = [str(s) for s in inputs]
                        crn_products = [str(s) for s in outputs]
                        sbml_products = [p.getSpecies()
                                         for p in sbml_reaction.getListOfProducts()]
                        sbml_reactants = [p.getSpecies()
                                          for p in sbml_reaction.getListOfReactants()]

                        # test that the reaction has the right inputs and outputs
                        assert set(crn_reactants) == set(sbml_reactants)
                        assert set(crn_products) == set(sbml_products)

                        if len(crn_reactants) > 0:
                            assert all(crn_reactants.count(s.getSpecies()) == int(
                                s.getStoichiometry()) for s in sbml_reaction.getListOfReactants())

                        if len(crn_products) > 0:
                            assert all(crn_products.count(s.getSpecies()) == int(
                                s.getStoichiometry()) for s in sbml_reaction.getListOfProducts())
                        assert not sbml_reaction.getReversible()
                        # TODO is there a smart way to test that the rate formula is correct?
                        # Validate the SBML model
                        assert validate_sbml(document) == 0
                    except Exception as e:  # collect errors to display
                        error_txt = f"Unexpected Error: in sbmlutil.add_reaction {rxn} for {model_id}. \n {str(e)}."
                        raise Exception(error_txt)

        k1_err = ParameterEntry("k1*", 1.11, key1)
        error_propensities = [MassAction(k_forward = k1_err), 
                            HillPositive(k=k1_err, n=2., K=3., s1=S1),
                            GeneralPropensity('k1*2 - k2/S1^2*', propensity_species=[S1], propensity_parameters=[k1, k2])]
        for prop in error_propensities:
            with pytest.raises(ValueError, match = "Could not write the rate law for reaction to SBML. Check the propensity functions of reactions."):
                rxn = Reaction([],[S1], propensity_type = prop)
                sbml_reaction = add_reaction(model, rxn, "r_err")
    def test_add_reaction_for_bioscrape(self):
        """
        Generates models for bioscrape including the particular annotations needed.
        """
        S1, S2, S3 = Species("S1"), Species("S2"), Species("S3")
        species = [S1, S2, S3]
        for_bioscrape = True  # Toggle for_bioscrape
        propensities = [MassAction(k_forward=1, k_reverse=.1), HillPositive(k=1, n=2., K=3., s1=S1),
                        HillNegative(k=1, n=2., K=3., s1=S1),
                        ProportionalHillPositive(k=1, n=2, K=3., s1=S1, d=S2),
                        ProportionalHillNegative(k=1, n=2, K=3., s1=S1, d=S2)]
        for prop in propensities:
            rxn_num = 0
            model_id = f"{prop.name}_model_with_for_bioscrape_{for_bioscrape}"
            document, model = create_sbml_model(model_id=model_id)
            add_all_species(model, species, {})
            # create a reaction
            rxn = Reaction([S1], [S2, S3], propensity_type=prop)
            try:
                # add reaction to SBML Model
                sbml_reaction = add_reaction(
                    model, rxn, f"r{rxn_num}", for_bioscrape=for_bioscrape)
                rxn_num += 1  # increment reaction id
                crn_reactants = [str(S1)]
                crn_products = [str(S2), str(S3)]
                sbml_products = [p.getSpecies()
                                 for p in sbml_reaction.getListOfProducts()]
                sbml_reactants = [p.getSpecies()
                                  for p in sbml_reaction.getListOfReactants()]

                # test that the reaction has the right inputs and outputs
                assert set(crn_reactants) == set(sbml_reactants)
                assert set(crn_products) == set(sbml_products)

                if len(crn_reactants) > 0:
                    assert all(crn_reactants.count(s.getSpecies()) == int(
                        s.getStoichiometry()) for s in sbml_reaction.getListOfReactants())

                if len(crn_products) > 0:
                    assert all(crn_products.count(s.getSpecies()) == int(
                        s.getStoichiometry()) for s in sbml_reaction.getListOfProducts())

                assert not sbml_reaction.getReversible()
                # Check annotations:
                sbml_annotation = sbml_reaction.getAnnotationString()
                check_var = f"type={prop.name}" in str(sbml_annotation)
                assert check_var == True
                # Test that the sbml annotation has keys for all species and parameters
                for k in prop.propensity_dict["parameters"]:
                    # convert k_reverse and k_forward to just k
                    k = k.replace("_reverse", "").replace("_forward", "")
                    check_var = f"{k}=" in sbml_annotation
                    assert check_var == True

                    # be sure that "k=" only shows up once in the annotation
                    assert sbml_annotation.count(f"{k}=") == 1

                for s in prop.propensity_dict["species"]:
                    check_var = f"{s}=" in sbml_annotation
                    assert check_var == True
                # TODO is there a smart way to test that the rate formula is correct?
                assert validate_sbml(document) == 0  # Validate the SBML model
            except Exception as e:  # collect errors
                error_txt = f"Unexpected Error: in sbmlutil.add_reaction {rxn} for {model_id}. \n {str(e)}."
                raise Exception(error_txt)


def test_generate_sbml_model():

    # Test a non-reversible reaction
    s1 = Species("S1")
    s2 = Species("S2")
    rx0 = Reaction.from_massaction(inputs=[s1], outputs=[s2], k_forward=0.1)
    crn = ChemicalReactionNetwork(species=[s1, s2], reactions=[rx0])

    # generate an sbml model
    document, model = crn.generate_sbml_model()
    # all species from the CRN are accounted for
    assert len(model.getListOfSpecies()) == len(crn.species)
    # all reactions from the CRN are accounted for
    assert len(model.getListOfReactions()) == len(crn.reactions)

    # reversible needs to be off!
    # assert not model.getListOfReactions()[0].isSetReversible()

    # test a reversible reaction
    rx1 = Reaction.from_massaction(
        inputs=[s1], outputs=[s2], k_forward=0.1, k_reverse=0.1)
    rxn_list = [rx1]
    crn = ChemicalReactionNetwork(species=[s1, s2], reactions=rxn_list)

    # generate an sbml model
    document, model = crn.generate_sbml_model()
    # all species from the CRN are accounted for
    assert len(model.getListOfSpecies()) == len(crn.species)
    # all reactions from the CRN are accounted for
    assert len(model.getListOfReactions()) == 2
    # although  sbml represents a reverisble reaction with reversible flag
    # BioCRNpyler always creates two reactions, because this is correct
    # for stochastic simulation with SBML.
    sbml_rxn = model.getListOfReactions()
    # assert not sbml_rxn[0].isSetReversible()
    # assert not sbml_rxn[1].isSetReversible()

    # Test propagation for for_bioscrape keyword
    # generate an sbml model with for_bioscrape = True
    document, model = crn.generate_sbml_model(for_bioscrape=True)
    for r in model.getListOfReactions():
        assert r.getAnnotation() is not None
    assert validate_sbml(document) == 0
    # generate an sbml model with for_bioscrape = False
    document, model = crn.generate_sbml_model(for_bioscrape=False)
    for r in model.getListOfReactions():
        assert r.getAnnotation() is None
    assert validate_sbml(document) == 0


def test_generate_sbml_model_parameter_names():

    s1 = Species("S1")
    s2 = Species("S2")

    # The correct parameter formating should be: "name_partid_mechanism"
    key0 = ParameterKey(mechanism="m", part_id="p", name="n")
    k0 = ParameterEntry("n", 1.0, parameter_key=key0)

    # Tests 3 Parameter entries with seemingly redundant keys.
    # In SBML parameter id, None will be kept blank, resulting in
    # different numbers of underscores between the V's
    k1 = ParameterEntry("v", 1.0, parameter_key=None)

    key2 = ParameterKey(mechanism=None, part_id="v", name="v")
    k2 = ParameterEntry("v", 2.0, parameter_key=key2)

    key3 = ParameterKey(mechanism="v", part_id=None, name="v")
    k3 = ParameterEntry("v", 2.0, parameter_key=key3)

    key4 = ParameterKey(mechanism="v", part_id="v", name="v")
    k4 = ParameterEntry("v", 2.0, parameter_key=key4)

    # Throw a duplicate key for good measure!
    rx0 = Reaction.from_massaction(
        inputs=[], outputs=[s1], k_forward=k0, k_reverse=k1)
    rx1 = Reaction.from_massaction(
        inputs=[s1], outputs=[s2], k_forward=k1, k_reverse=k2)
    rx2 = Reaction.from_massaction(
        inputs=2*[s1], outputs=2*[s2], k_forward=k3, k_reverse=k4)
    rxn_list = [rx0, rx1, rx2]
    crn = ChemicalReactionNetwork(species=[s1, s2], reactions=rxn_list)

    document, model = crn.generate_sbml_model()
    assert validate_sbml(document) == 0
    correct_ids = set(["v_v_", "v__v", "v_v_v", "v__", "n_p_m"])
    ids = set([p.getId() for p in model.getListOfParameters()])
    assert ids == correct_ids


def test_sbml_basics():
    """
    This test aims to test all Python libSBML functions and their return values. 
    It will catch errors if certain libSBML functions are not working as expected. 
    For example: At times, if the reaction rate law string is mis-formatted, it is 
    ignored by the SBML writing code and not caught by any other tests. Similarly, this 
    test will catch modifier reference, parameter values etc. 
    """
    def check(value, message):
        """If 'value' is None, prints an error message constructed using
        'message' and then exits with status code 1 (for libsbml). If 'value' is an integer,
        it assumes it is a libSBML return status code.  If the code value is
        LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
        prints an error message constructed using 'message' along with text from
        libSBML explaining the meaning of the code, and exits with status code 1.
        """
        if value == None:
            raise SystemExit(
                'LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value == libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                    + 'LibSBML returned error code ' + str(value) + ': "' \
                    + libsbml.OperationReturnValue_toString(value).strip() + '"'
                raise SystemExit(err_msg)
        else:
            return

    from biocrnpyler.sbmlutil import _create_global_parameter, _create_local_parameter
    from biocrnpyler.sbmlutil import _create_modifiers, _create_products, _create_reactants
    # generate an sbml model
    document, model = create_sbml_model()

    # If the document/model creation failed the following two lines
    # will catch it:
    check(document, 'create SBMLDocument object')
    check(model, 'create SBML Model object')
    S1 = Species("S1")
    S2 = Species("S2")
    compartment = add_compartment(model, S1.compartment)
    check(compartment, 'add compartment')
    sbml_species = add_species(
        model, compartment, S1, initial_concentration=10)
    sbml_species2 = add_species(
        model, compartment, S2, initial_concentration=10)
    check(sbml_species.getId(), 'get ID for SBML species')
    check(sbml_species.getInitialConcentration(),
          'get concentration for SBML species')
    check(sbml_species2.getId(), 'get ID for SBML species')
    check(sbml_species2.getInitialConcentration(),
          'get concentration for SBML species')

    prop_hill = HillPositive(k=1, s1=S2, K=10, n=2)
    crn_rxn = Reaction([S1], [], propensity_type=prop_hill)
    check(model.createReaction(), 'create a new SBML reaction')
    rx1 = model.getReaction(0)
    check(rx1, 'get the created reaction')
    check(rx1.setId('r1'), 'set ID for reaction')
    check(rx1.setReversible(False), 'set reversible attribute for reaction')
    _create_reactants(crn_rxn.inputs, rx1, model)
    _create_products(crn_rxn.outputs, rx1, model)
    _create_modifiers(crn_rxn, rx1, model)

    check(rx1.getModifier(0), 'create species modifier')
    modifier = rx1.getModifier(0)
    check(modifier.getSpecies(), 'get the species id for the modifier reference')
    rateLaw = rx1.createKineticLaw()
    check(rateLaw, 'create kineticLaw for reaction')
    local_param = _create_local_parameter(rateLaw, 'k_local', 10)
    check(local_param, 'create local parameter in SBML model')
    global_param = _create_global_parameter(model, 'k_global', 10)
    check(global_param, 'create global parameter in SBML model')
    with pytest.raises(ValueError, match = "Units for a parameter must be passed as strings."):
        global_param = _create_global_parameter(model, 'k_global', value = 10, p_unit = 24)
    with pytest.warns(Warning, match = "The string identifier for the unit 1_s is not supported by BioCRNpyler. Add this to the dictionary in biocrnpyler/units.py if you want this unit."):
        global_param = _create_global_parameter(model, 'k_global', value = 10, p_unit = "1_s")
    check(rateLaw.setFormula('k_local*10 - k_global/2'),
          'set rate formula for reaction')
    validator = validateSBML(ucheck=False)
    validation_result = validator.validate(document, print_results=True)
    if validation_result > 0:
        raise Exception('Invalid SBML model.')
