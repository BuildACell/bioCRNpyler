#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from unittest import TestCase
from biocrnpyler import Species, Mixture, Component, ParameterKey, ChemicalReactionNetwork, ChemicalComplex



def test_initial_condition_in_crn():
    S = Species("S")

    CRN = ChemicalReactionNetwork(species = [S], reactions = [])

    assert S not in CRN.initial_condition_dict

    CRN.initial_condition_dict = {S:10}
    assert S in CRN.initial_condition_dict
    assert CRN.initial_condition_dict[S] == 10


def test_species_initial_condition_in_component():

    S1 = Species("S1")
    S2 = Species("S2")

    key1 = ParameterKey(mechanism = "initial concentration", part_id = None, name = str(S1))
    key2 = ParameterKey(mechanism = "initial concentration", part_id = None, name = "C")

    C = ChemicalComplex([S1, S2], name = "C", parameters = {key1:10, key2:2.2})
    S3 = C.get_species()

    M = Mixture(name = "M", components = [C])

    #Initial condition found under the Species name
    assert M.get_initial_condition(S1, C)[S1].value == 10

    #Initial condition defaults to 0
    assert M.get_initial_condition(S2, C)[S2] == 0

    #Initial condition found under the Component name
    assert M.get_initial_condition(S3, C)[S3].value == 2.2


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
    assert M.get_initial_condition(S1)[S1].value == 10

    #Initial condition defaults to 0
    assert M.get_initial_condition(S2)[S2] == 0

    #Initial condition found under the Component name
    assert M.get_initial_condition(S3, C)[S3].value == 2.2

def test_sppecies_initial_condition_defaulting():

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
    assert M.get_initial_condition(S1, C)[S1].value == .11
    assert M.get_initial_condition(S2, C)[S2].value == .22
    assert M.get_initial_condition(S3, C)[S3].value == .33





	

