#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List
from .chemical_reaction_network import Species


class Propensity(object):
    def __init__(self):
        pass

    @staticmethod
    def is_valid_propensity(propensity_type) -> bool:
        """checks whether the given propensity_type is valid propensity

        It recursively checks all subclasses of Propensity until it finds the
        propensity type. Otherwise raise a Type error

        :param propensity_type: Propensity
        :returns bool
        """
        pass

    def check_species(self, species_list: List[Species]) -> bool:
        required_species = self.get_species()
        if required_species:
            found = False
            for species in species_list:
                if species in required_species:
                    found = True
            return  found
        else:
            # the current propensity doesn't require any species
            return True

    def get_species(self):
        """returns the instance variables that are species type"""
        return []


class GeneralPropensity(Propensity):
    def __init__(self, propensity_function: str):
        super(GeneralPropensity, self).__init__()
        self.propensity_function  = propensity_function


class MassAction(Propensity):
    def __init__(self, k: float, k_rev: float = None):
        super(MassAction, self).__init__()
        self.k = k
        self.k_rev = k_rev


class HillPositive(MassAction):
    def __init__(self, k: float, s1: Species, K: float, n: float):
        """ Hill positive propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*s1^n/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        """
        super(HillPositive, self).__init__(k=k)
        self.s1 = s1
        self.K = K


class HillNegative(MassAction):
    def __init__(self, k: float, s1: Species, K: float, n: float):
        """ Hill negative propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*1/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        """
        super(HillNegative, self).__init__(k=k)
        self.s1 = s1
        self.K = K
        self.n = n


class ProportionalHillPositive(HillPositive):
    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        """ proportional Hill positive propensity with the following formula

            p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        super(ProportionalHillPositive, self).__init__(k=k, s1=s1, K=K, n=n)
        self.d = d


class PropotionalHillNegative(HillNegative):
    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        """ proportional Hill positive propensity with the following formula

            p(s1, d; k, K, n) = k*d/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        super(PropotionalHillNegative, self).__init__(k=k, s1=s1, K=K, n=n)
        self.d = d
