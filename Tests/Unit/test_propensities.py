#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from biocrnpyler import MassAction, ProportionalHillNegative, ProportionalHillPositive, HillPositive, HillNegative, Hill
from biocrnpyler import GeneralPropensity
from biocrnpyler import ParameterEntry
from biocrnpyler import Species
import pytest
from biocrnpyler import sbmlutil


def test_massaction_forward_rate():

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=0)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=-1)


def test_massaction_reserve_rate():
    with pytest.raises(TypeError, match=r"missing 1 required positional argument: 'k_forward'"):
        MassAction(k_reverse=0.1)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=1, k_reverse=0)

    with pytest.raises(ValueError, match=r"Propensity parameters must be Parameters or floats with positive values.*"):
        MassAction(k_forward=1, k_reverse=-1)


def test_massaction_is_reverable():
    mak = MassAction(k_forward=1, k_reverse=0.1)
    assert mak.is_reversible


def test_massaction_from_parameters():
    mak = MassAction(k_forward=ParameterEntry("k", .1))
    assert mak.k_forward == .1

    mak = MassAction(k_forward=ParameterEntry("k", .1), k_reverse=ParameterEntry("k", .01))
    assert mak.k_reverse == .01


def test_hill_negative_init():
    s1 = Species('s1')
    hillneg = HillNegative(k=1, s1=s1, K=5, n=2)
    assert s1 == hillneg.s1

    assert hillneg.is_reversible is False

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        HillNegative(k=1, s1='s1', K=5, n=2)

    with pytest.raises(TypeError):
        d = Species('d')
        HillNegative(k=1, s1=s1, K=5, n=2, d=d)


def test_hill_negative_rate_formula():
    in_s1 = Species('s1')
    in_k = 1
    in_K = 5
    in_n = 2
    philpos = HillNegative(k=in_k, s1=in_s1, K=in_K, n=in_n)

    assert philpos._get_rate_formula(philpos.propensity_dict) == f"{in_k} / ( {in_K}^{in_n} + {in_s1}^{in_n} )"


def test_hill_positive_init():
    s1 = Species('s1')
    philpos = HillPositive(k=1, s1=s1, K=5, n=2)
    assert s1 == philpos.s1

    assert philpos.is_reversible is False

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        HillPositive(k=1, s1='s1', K=5, n=2)

    with pytest.raises(TypeError):
        d = Species('d')
        HillPositive(k=1, s1=s1, K=5, n=2, d=d)


def test_hill_positive_rate_formula():
    in_s1 = Species('s1')
    in_k = 1
    in_K = 5
    in_n = 2
    philpos = HillPositive(k=in_k, s1=in_s1, K=in_K, n=in_n)

    assert philpos._get_rate_formula(philpos.propensity_dict) == f"{in_k}*{in_s1}^{in_n} / ( {in_K}^{in_n} + {in_s1}^{in_n} )"


def test_proportional_hill_positive_init():
    d = Species('d')
    s1 = Species('s1')
    philpos = ProportionalHillPositive(k=1, s1=s1, K=5, n=2, d=d)
    assert d == philpos.d
    assert s1 == philpos.s1

    assert philpos.is_reversible is False

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        ProportionalHillPositive(k=1, s1=s1, K=5, n=2, d='')

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        ProportionalHillPositive(k=1, s1='s1', K=5, n=2, d=d)


def test_proportional_hill_positive_rate_formula():
    in_d = Species('d')
    in_s1 = Species('s1')
    in_k = 1
    in_K = 5
    in_n = 2
    philpos = ProportionalHillPositive(k=in_k, s1=in_s1, K=in_K, n=in_n, d=in_d)

    assert philpos._get_rate_formula(philpos.propensity_dict) == f"{in_k}*{in_d}*{in_s1}^{in_n} / ( {in_K}^{in_n} + {in_s1}^{in_n} )"


def test_proportional_hill_negative_init():
    d = Species('d')
    s1 = Species('s1')
    philneg = ProportionalHillNegative(k=1, s1=s1, K=5, n=2, d=d)
    assert d == philneg.d
    assert s1 == philneg.s1

    assert philneg.is_reversible is False

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        ProportionalHillNegative(k=1, s1=s1, K=5, n=2, d='')

    with pytest.raises(TypeError, match='Propensity expected a Species'):
        ProportionalHillNegative(k=1, s1='s1', K=5, n=2, d=d)


def test_proportional_hill_negative_rate_formula():
    in_d = Species('d')
    in_s1 = Species('s1')
    in_k = 1
    in_K = 5
    in_n = 2
    philneg = ProportionalHillNegative(k=in_k, s1=in_s1, K=in_K, n=in_n, d=in_d)

    assert philneg._get_rate_formula(philneg.propensity_dict) == f"{in_k}*{in_d} / ( {in_K}^{in_n} + {in_s1}^{in_n} )"


def test_general_propensity():
    S1, S2, S3 = Species("S1"), Species("S2"), Species("S3")
    # create some parameters
    k1 = ParameterEntry("k1", 1.11)
    k2 = ParameterEntry("k2", 2.22)

    gn1 = GeneralPropensity('k1*2 - k2/S1^2', propensity_species=[S1], propensity_parameters=[k1, k2])
    assert str(S1) in gn1.propensity_dict['species']
    assert k1.parameter_name in gn1.propensity_dict['parameters']
    assert k2.parameter_name in gn1.propensity_dict['parameters']
    assert gn1.pretty_print() == 'k1*2 - k2/S1^2\n  k1=1.11\n  k2=2.22\n'

    gn2 = GeneralPropensity('S1^2 + S2^2 + S3^2', propensity_species=[S1, S2, S3], propensity_parameters=[])
    assert str(S1) in gn1.propensity_dict['species']
    assert k1.parameter_name not in gn2.propensity_dict['parameters']

    with pytest.raises(TypeError, match='propensity_species must be a list of Species!'):
        GeneralPropensity('S1^2 + S2^2 + S3^2', propensity_species=[k1], propensity_parameters=[])

    with pytest.raises(TypeError, match='propensity_parameter must be a list of ParameterEntry!'):
        GeneralPropensity('S1^2 + S2^2 + S3^2', propensity_species=[], propensity_parameters=[S2])

    test_formula = 'S3^2'
    with pytest.raises(ValueError, match=f'must be part of the formula'):
        GeneralPropensity(test_formula, propensity_species=[S1], propensity_parameters=[])

    test_formula = 'k2*S3^2'
    with pytest.raises(ValueError, match=f'must be part of the formula'):
        GeneralPropensity(test_formula, propensity_species=[S3], propensity_parameters=[k1])


def test_propensity_dict_massaction():
    k1 = ParameterEntry(parameter_value = '1', parameter_name = 'k1')
    k2 = ParameterEntry(parameter_value = '2', parameter_name = 'k2')

    #Should store the ParameterEntry in this case
    P1 = MassAction(k_forward = k1, k_reverse = k2)
    assert P1.propensity_dict["parameters"]["k_forward"] == k1
    assert P1.propensity_dict["parameters"]["k_reverse"] == k2

    #assert getters work (should return values instead of ParameterEntries)
    assert P1.k_forward == k1.value
    assert P1.k_reverse == k2.value

    #Should store a numerical value in this case
    P2 = MassAction(k_forward = k1.value, k_reverse = k2.value)
    assert P2.propensity_dict["parameters"]["k_forward"] == k1.value
    assert P2.propensity_dict["parameters"]["k_reverse"] == k2.value

    #assert getters work
    assert P2.k_forward == k1.value
    assert P2.k_reverse == k2.value


def test_propensity_dict_hill():
    d = Species('d')
    s1 = Species('s1')
    k = ParameterEntry(parameter_value = '1', parameter_name = 'k')
    K = ParameterEntry(parameter_value = '2', parameter_name = 'K')
    n = ParameterEntry(parameter_value = '3', parameter_name = 'n')

    #Should store the ParameterEntry in this case
    P1 = Hill(k = k, K = K, n = n, s1 = s1, d = d)
    assert P1.propensity_dict["parameters"]["k"] == k
    assert P1.propensity_dict["parameters"]["K"] == K
    assert P1.propensity_dict["parameters"]["n"] == n
    #assert the getters work (should return values instead of ParameterEntries)
    assert P1.k == k.value
    assert P1.K == K.value
    assert P1.n == n.value

    #Should store a numerical value in this case
    P2 = Hill(k = k.value, K = K.value, n = n.value, s1 = s1, d = d)
    assert P2.propensity_dict["parameters"]["k"] == k.value
    assert P2.propensity_dict["parameters"]["K"] == K.value
    assert P2.propensity_dict["parameters"]["n"] == n.value
    #assert the getters work
    assert P2.k == k.value
    assert P2.K == K.value
    assert P2.n == n.value