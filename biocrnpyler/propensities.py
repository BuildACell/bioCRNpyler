#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List


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

    @staticmethod
    def get_available_propensities() -> List[str]:
        pass

    def check_species(self, species_list) -> bool:
        required_species = self.get_species()
        if required_species:
            found = False
            for species in species_list:
                if species in required_species:
                    found = True
            return found
        else:
            # the current propensity doesn't require any species
            return True

    def get_species(self):
        """returns the instance variables that are species type"""
        return []

    @property
    def is_reversible(self):
        raise NotImplementedError

    @property
    def k_forward(self):
        raise NotImplementedError

    @property
    def k_reverse(self):
        raise NotImplementedError


class GeneralPropensity(Propensity):
    def __init__(self, propensity_function: str):
        super(GeneralPropensity, self).__init__()
        self.propensity_function = propensity_function


class MassAction(Propensity):
    def __init__(self, k_forward: float, k_reverse: float = None):
        super(MassAction, self).__init__()
        self.k_forward = k_forward
        self.k_reverse = k_reverse

    @property
    def k_forward(self):
        return self._k_forward

    @k_forward.setter
    def k_forward(self, new_k_forward):
        if new_k_forward <= 0:
            raise ValueError(f'Forward reaction rate coefficient is negative! '
                             f'k={new_k_forward}')

        self._k_forward = new_k_forward

    @property
    def k_reverse(self):
        return self._k_reverse

    @k_reverse.setter
    def k_reverse(self, new_k_reverse):
        if new_k_reverse is not None and new_k_reverse <= 0:
            raise ValueError(f'Reverse reaction rate coefficient is negative! '
                             f'k={new_k_reverse}')
        self._k_reverse = new_k_reverse

    @property
    def is_reversible(self):
        if self.k_reverse is None:
            return False
        else:
            return True


class HillPositive(MassAction):
    def __init__(self, k_forward: float, s1, K: float, n: float):
        """ Hill positive propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*s1^n/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        """
        super(HillPositive, self).__init__(k_forward=k_forward)
        self.s1 = s1
        self.K = K


class HillNegative(MassAction):
    def __init__(self, k_forward: float, s1, K: float, n: float):
        """ Hill negative propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*1/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        """
        super(HillNegative, self).__init__(k_forward=k_forward)
        self.s1 = s1
        self.K = K
        self.n = n


class ProportionalHillPositive(HillPositive):
    def __init__(self, k_forward: float, s1, K: float, n: float, d):
        """ proportional Hill positive propensity with the following formula

            p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        super(ProportionalHillPositive, self).__init__(k_forward=k_forward, s1=s1, K=K, n=n)
        self.d = d


class ProportionalHillNegative(HillNegative):
    def __init__(self, k_forward: float, s1, K: float, n: float, d):
        """ proportional Hill positive propensity with the following formula

            p(s1, d; k, K, n) = k*d/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        super(ProportionalHillNegative, self).__init__(k_forward=k_forward, s1=s1, K=K, n=n)
        self.d = d
