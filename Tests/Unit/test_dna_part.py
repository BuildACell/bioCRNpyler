import pytest
from biocrnpyler import DNA_part, Operator, UserDefined

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

    partfor = DNA_part("part",direction="forward")
    x = part1.set_dir("forward")
    assert partfor == x

def test_misc_parts():
    p1 = Operator("op1",binder=["b1"])

    assert(p1.name=="op1")
    assert(len(p1.binder)==1)
    assert(p1.binder[0].name=="b1")

    u1 = UserDefined("op2",dpl_type="Origin")

    assert(u1.name=="op2")
    assert(u1.dpl_type=="Origin")