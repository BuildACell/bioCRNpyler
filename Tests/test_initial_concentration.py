#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from unittest import TestCase
from biocrnpyler import Species, Mixture, ExpressionExtract, Component, ParameterKey, ChemicalReactionNetwork, ChemicalComplex, parameter_to_value



def test_initial_condition_in_crn():
    S = Species("S")
    S2 = Species("S2")
    S3 = Species("S3")

    CRN = ChemicalReactionNetwork(species = [S, S2, S3], reactions = [])

    #No Initial Concentration Dict
    assert S not in CRN.initial_concentration_dict

    #test setter
    CRN.initial_concentration_dict = {S:10, S2:10}
    assert S in CRN.initial_concentration_dict
    assert CRN.initial_concentration_dict[S] == 10

    #Test setter overwriting
    CRN.initial_concentration_dict = {S:11, S3:10}
    assert S2 in CRN.initial_concentration_dict
    assert CRN.initial_concentration_dict[S] == 11
    assert CRN.initial_concentration_dict[S2] == 10
    assert CRN.initial_concentration_dict[S3] == 10

    #Test setter reset
    CRN.initial_concentration_dict = None
    assert S not in CRN.initial_concentration_dict
    assert S2 not in CRN.initial_concentration_dict
    assert S3 not in CRN.initial_concentration_dict


def test_species_initial_condition_in_component():

    S1 = Species("S1")
    S2 = Species("S2")

    key1 = ParameterKey(mechanism = "initial concentration", part_id = None, name = str(S1))
    key2 = ParameterKey(mechanism = "initial concentration", part_id = None, name = "C")

    C = ChemicalComplex([S1, S2], name = "C", parameters = {key1:10, key2:2.2})
    S3 = C.get_species()

    M = Mixture(name = "M", components = [C])

    #Initial condition found under the Species name
    assert parameter_to_value(M.get_initial_concentration(S1, C)[S1]) == 10

    #Initial condition defaults to 0
    assert parameter_to_value(M.get_initial_concentration(S2, C)[S2]) == 0

    #Initial condition found under the Component name
    assert parameter_to_value(M.get_initial_concentration(S3, C)[S3]) == 2.2


def test_species_initial_condition_in_mixture():

    S1 = Species("S1")
    S2 = Species("S2")

    C = ChemicalComplex([S1, S2], name = "C")
    S3 = C.get_species()

    mixture_name = "M"
    key1 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S1))
    key2 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = C.name)

    M = Mixture(name = mixture_name, components = [C], parameters = {key1:10, key2:2.2})

    #Initial condition found under the Species name
    assert parameter_to_value(M.get_initial_concentration(S1)[S1]) == 10

    #Initial condition defaults to 0
    assert parameter_to_value(M.get_initial_concentration(S2)[S2]) == 0

    #Initial condition found under the Component name
    assert parameter_to_value(M.get_initial_concentration(S3, C)[S3]) == 2.2

def test_species_initial_condition_defaulting():

    S1 = Species("S1")
    S2 = Species("S2")

    C = ChemicalComplex([S1, S2], name = "C")
    S3 = C.get_species()

    mixture_name = "M"
    key1 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S1))
    key2 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S2))
    key3 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S3))

    C = ChemicalComplex([S1, S2], name = "C", parameters = {key1:.11, key2:.22}, initial_concentration = .33)
    M = Mixture(name = mixture_name, components = [C], parameters = {key1:1.1, key2:2.2, key3:3.3})

    #Initial condition found under the Species name, in the Component, not Mixture
    assert parameter_to_value(M.get_initial_concentration(S1, C)[S1]) == .11
    assert parameter_to_value(M.get_initial_concentration(S2, C)[S2]) == .22
    assert parameter_to_value(M.get_initial_concentration(S3, C)[S3]) == .33


def test_species_initial_condition_at_compile():

    S1 = Species("S1")
    S2 = Species("S2")

    C = ChemicalComplex([S1, S2], name = "C")
    S3 = C.get_species()

    mixture_name = "M"
    key1 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S1))
    key2 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S2))
    key3 = ParameterKey(mechanism = "initial concentration", part_id = mixture_name, name = str(S3))

    C = ChemicalComplex([S1, S2], name = "C", parameters = {key1:.11, key2:.22}, initial_concentration = .33)
    M = ExpressionExtract(name = mixture_name, components = [C], parameters = {key1:1.1, key2:2.2, key3:3.3, "kb":1.0, "ku":1.0})

    #Initial condition found under the Species name, in the Component, not Mixture
    CRN = M.compile_crn()
    #Initial condition found under the Species name, in the Component, not Mixture
    assert parameter_to_value(CRN.initial_concentration_dict[S1]) == .11
    assert parameter_to_value(CRN.initial_concentration_dict[S2]) == .22
    assert parameter_to_value(CRN.initial_concentration_dict[S3]) == .33

    #These new conditions should overide the old ones
    CRN2 = M.compile_crn(initial_concentration_dict = {S1:1.1, S2:2.2})
    assert parameter_to_value(CRN2.initial_concentration_dict[S1]) == 1.1
    assert parameter_to_value(CRN2.initial_concentration_dict[S2]) == 2.2
    assert parameter_to_value(CRN2.initial_concentration_dict[S3]) == .33





	

