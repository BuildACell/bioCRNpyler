#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest
from biocrnpyler import Promoter, Species


def test_promoter_from_promoter():
    test_name = 'test_promoter'
    test_length = 10
    prom = Promoter(name=test_name, length=test_length)

    assert prom.assembly is None
    assert prom.transcript is None
    assert prom.protein is None

    test_assembly = 'test_assembly'
    test_transcript = 'test_transcript'
    test_protein = 'test_protein'
    prom2= Promoter.from_promoter(name=prom,assembly=test_assembly, transcript=test_transcript, protein=test_protein)

    assert prom2.assembly == test_assembly
    assert prom2.transcript == test_transcript
    assert prom2.protein == test_protein
    # check that original values are kept
    assert prom2.name == test_name
    assert prom2.length == test_length

    with pytest.raises(TypeError, match='Promoter can be initialized from string or another promoter!'):
        Promoter.from_promoter(name=Species(name=test_name), assembly=test_assembly, transcript=test_transcript, protein=test_protein)
