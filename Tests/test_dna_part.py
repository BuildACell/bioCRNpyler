import pytest
from biocrnpyler import DNA_part

def test_instantiation():

    part = DNA_part("part 1")
    assert part.name == "part 1"

    with pytest.raises(AttributeError, match=r"DNA_part should not recieve initial_concentration keyword. Pass this into the DNAassembly or DNA_construct instead*"):
        part2 = DNA_part("part 2", initial_concentration = 10)


    part3 = DNA_part("part3", sequence = "abc", assembly = "a1")

    assert part3.sequence == "abc"
    assert part3.assembly == "a1"

def test_equality_without_ordered_polymer():

    part1 = DNA_part("part")
    part1b = DNA_part("part")
    part2 = DNA_part("part2")
    assert part1 == part1b
    assert part1 != part2

    part1a = DNA_part("part", assembly = "a1")
    assert part1 != part1a
    part1a2 = DNA_part("part", assembly = "a1")
    assert part1a == part1a2
    part1b = DNA_part("part", assembly = "b1")
    assert part1a != part1b
