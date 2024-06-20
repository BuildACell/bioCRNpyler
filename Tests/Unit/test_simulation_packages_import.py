#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import warnings
import pytest


def test_bioscrape_import_simulate():
    from biocrnpyler import ChemicalReactionNetwork, Species
    try:
        import numpy as np

        X = Species("X")
        CRN = ChemicalReactionNetwork(species=[X],reactions=[])
        with pytest.warns() as record:
            sim_result = CRN.simulate_with_bioscrape(timepoints=np.linspace(0, 10, 100), initial_condition_dict = {str(X):1})


        #In the case bioscrape is not imported
        if sim_result is None:
            # check the warning message
            assert "simulate_with_bioscrape is depricated and will cease working in a future release." in str(record[0].message)
        #In the case bioscrape is imported
        else:
            assert str(X) in sim_result
    except ModuleNotFoundError:
        print('test skipped')


def test_bioscrape_import_simulate_via_sbml():
    from biocrnpyler import ChemicalReactionNetwork, Species
    try:
        import numpy as np

        X = Species("X")
        CRN = ChemicalReactionNetwork(species=[X],reactions=[])
        with warnings.catch_warnings(record=True) as record:
            sim_result, bioscrape_model = CRN.simulate_with_bioscrape_via_sbml(timepoints=np.linspace(0, 10, 100), initial_condition_dict={str(X):1}, return_model=True)


        #In the case bioscrape is not imported
        if sim_result is None and bioscrape_model is None:
            # only one warning was triggered
            assert len(record) == 1
            # check the warning message
            assert str(record[0].message) == "bioscrape was not found, please install bioscrape"

        #In the case bioscrape is imported
        else:
            assert str(X) in sim_result
            assert bioscrape_model is not None
    except ModuleNotFoundError:
        print('test skipped')


def test_libroadrunner_import():
    from biocrnpyler import ChemicalReactionNetwork
    CRN = ChemicalReactionNetwork(species=[], reactions=[])
    try:
        import roadrunner
    except ModuleNotFoundError:
        # libroadrunner is not installed, let's check if it triggers a warning inside simulate_with_roadrunner()
        with pytest.warns(UserWarning, match='libroadrunner was not found, please install libroadrunner'):
            CRN.simulate_with_roadrunner(timepoints=list(range(0,10)))
    else:
        # no exception was triggered, we can simulate the CRN with roadrunner
        sim_results = CRN.simulate_with_roadrunner(timepoints=list(range(0,10)))
        assert sim_results is not None




