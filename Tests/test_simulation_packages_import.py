#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest


def test_bioscrape_import_simulate():
    from biocrnpyler import ChemicalReactionNetwork

    CRN = ChemicalReactionNetwork(species=[],reactions=[])
    with pytest.warns(None) as record:
        sim_result = CRN.simulate_with_bioscrape(timepoints=[], initial_condition_dict = {})

    assert sim_result is None

    # only one warning was triggered
    assert len(record) == 1
    # check the warning message
    assert str(record[0].message) == "bioscrape was not found, please install bioscrape"


def test_bioscrape_import_simulate_via_sbml():
    from biocrnpyler import ChemicalReactionNetwork

    CRN = ChemicalReactionNetwork(species=[],reactions=[])
    with pytest.warns(None) as record:
        sim_result, bioscrape_model  = CRN.simulate_with_bioscrape_via_sbml(timepoints=[], initial_condition_dict={})

    assert sim_result is None
    assert bioscrape_model is None

    # only one warning was triggered
    assert len(record) == 1
    # check the warning message
    assert str(record[0].message) == "bioscrape was not found, please install bioscrape"


def test_libroadrunner_import():
    from biocrnpyler import ChemicalReactionNetwork

    CRN = ChemicalReactionNetwork(species=[],reactions=[])
    with pytest.warns(None) as record:
        sim_results = CRN.runsim_roadrunner(timepoints=[], filename=None)

    assert sim_results is None

    # only one warning was triggered
    assert len(record) == 1
    # check the warning message
    assert str(record[0].message) == "libroadrunner was not found, please install libroadrunner"

