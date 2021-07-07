import pytest
from biocrnpyler import UserDefined,Origin,DNABindingSite

def test_user_defined():

    part = UserDefined("part 1")
    assert part.name == "part 1"

    assert part.update_species()==[]
    assert part.update_reactions()==[]
def test_origin():
    part = Origin("part 1")
    assert part.name == "part 1"

    assert part.update_species()==[]
    assert part.update_reactions()==[]
def test_dnabindingsite():
    part = DNABindingSite("part 1",["testprotein"])
    assert part.name == "part 1"