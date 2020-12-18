#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import RBS, Species


def test_promoter_from_rbs():
    test_name = 'test_rbs'
    test_length = 10
    rbs = RBS(name=test_name, length=test_length)

    assert rbs.assembly is None
    assert rbs.transcript is None
    assert rbs.protein is None

    test_assembly = 'test_assembly'
    test_transcript = 'test_transcript'
    test_protein = 'test_protein'
    rbs2= RBS.from_rbs(name=rbs,assembly=test_assembly, transcript=test_transcript, protein=test_protein)

    assert rbs2.assembly == test_assembly
    assert rbs2.transcript == test_transcript
    assert rbs2.protein == test_protein
    # check that original values are kept
    assert rbs2.name == test_name
    assert rbs2.length == test_length

    with pytest.raises(TypeError, match='RBS can be initialized from string or another RBS!'):
        rbs.from_rbs(name=Species(name=test_name), assembly=test_assembly, transcript=test_transcript, protein=test_protein)
