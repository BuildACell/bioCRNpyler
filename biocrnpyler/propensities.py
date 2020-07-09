#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import Set, Union
from numbers import Number
from .parameter import Parameter


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
        for propensity in Propensity.get_available_propensities():
            if isinstance(propensity_type, propensity):
                return True
        return False

    @staticmethod
    def _all_subclasses(cls):
        """Returns a set of all subclasses of cls (recursively calculated)

        Source:
        https://stackoverflow.com/questions/3862310/how-to-find-all-the-subclasses-of-a-class-given-its-name
        :param cls: A class in the codebase, for example Propensity
        :return: set of all subclasses from cls
        """
        return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in Propensity._all_subclasses(c)])

    @staticmethod
    def get_available_propensities() -> Set:
        return Propensity._all_subclasses(Propensity)

    @staticmethod
    def _get_parameter_name(rate_coeff: Union[Number, Parameter]) -> str:
        """checks whether a rate_coeff is a Number, if not this function calls for its string represention

        :param rate_coeff: a number of a parameter object
        :return: either None or the parameter's name
        """
        if not isinstance(rate_coeff, Number):
            rate_coeff_name = str(rate_coeff)
        else:
            rate_coeff_name = 'k'
        return rate_coeff_name

    def pretty_print(self, **kwargs):
        raise NotImplementedError

    @property
    def is_reversible(self):
        raise NotImplementedError

    @property
    def k_forward(self):
        raise NotImplementedError

    @property
    def k_reverse(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError

    @property
    def species(self):
        """returns the instance variables that are species type"""
        raise NotImplementedError

    def create_kinetic_law(self, reaction, reverse_reaction, stochastic):
        raise NotImplementedError


class GeneralPropensity(Propensity):
    def __init__(self, propensity_function: str):
        super(GeneralPropensity, self).__init__()
        self.propensity_function = propensity_function

    def create_kinetic_law(self, reaction):
        raise NotImplementedError("SBML writing of general propensities not implemented")


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
        if not isinstance(self, MassAction):
            # only mass action object can have reverse rates
            raise NotImplementedError
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

    def pretty_print(self, **kwargs):
        if 'reaction_direction' in kwargs and kwargs['react_direction'] is 'reverse':
            if self.is_reversible:
                return f''
            else:
                return None
        else:
            return f''

    def __eq__(self, other):
        return self.k_forward == other.k_forward and self.k_reverse == other.k_reverse

    @property
    def species(self):
        return []

    def create_kinetic_law(self, reaction, stochastic, direction='forward'):
        # Create a kinetic law for the reaction
        ratelaw = reaction.createKineticLaw()

        annotation_dict = {"type": self}
        if direction == 'forward':
            rate_coeff = self.k_forward
        else:
            rate_coeff = self.k_reverse

        # get the name of the kinetic rate coefficient
        rate_coeff_name = Propensity._get_parameter_name(rate_coeff)

        param = ratelaw.createParameter()
        param.setId(rate_coeff_name)
        param.setConstant(True)
        param.setValue(rate_coeff)
        annotation_dict["k"] = rate_coeff

        # Create Rate-strings for massaction propensities
        ratestring = rate_coeff_name

        if stochastic:
            for i in range(stoichiometry):
                if i > 0:
                    ratestring += f" * ( {species_id} - {i} )"
                else:
                    ratestring += f" * {species_id}"

        else:
                ratestring += f" * {species_id}^{stoichiometry}"


        # Set the ratelaw to the ratestring
        math_ast = libsbml.parseL3Formula(ratestring)
        ratelaw.setMath(math_ast)

        if propensity_annotation:
            annotation_string = "<PropensityType>"
            for k in annotation_dict:
                annotation_string += " "+k + "=" + str(annotation_dict[k])
            annotation_string += "</PropensityType>"
            reaction.appendAnnotation(annotation_string)

        return ratelaw


class Hill(MassAction):
    def __init__(self, k_forward: float, s1, K: float, n: float):
        super(Hill, self).__init__(k_forward=k_forward)
        
        self.s1 = s1
        self.K = K
        self.n = n

    def __eq__(self, other):
        return super(Hill, self).__eq__(other) \
               and self.s1 is other.s1 \
               and self.K == other.K \
               and self.n == other.n

    @property
    def species(self):
        return [self.s1]
        
    def create_kinetic_law(self, reaction, direction='forward'):
        if direction != 'forward':
            raise ValueError('Only MassAction can have reverse reaction!')

        ratelaw = super(Hill, self).create_kinetic_law(reaction=reaction)

        param_n = ratelaw.createParameter()
        param_n.setId("n")
        param_n.setConstant(True)
        param_n.setValue(self.n)
    
        param_K = ratelaw.createParameter()
        param_K.setId("K")
        param_K.setConstant(True)
        param_K.setValue(self.K)

        annotation_dict["k"] = self.k_forward
        annotation_dict["K"] = self.K
        annotation_dict["n"] = self.n


        #Create ratestring for non-massaction propensities
    if propensity_type == "hillpositive":

        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()

        ratestring+=f"*{s_species_id}^n/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id

    elif propensity_type == "hillnegative":
        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()

        ratestring+=f"/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id

    elif propensity_type == "proportionalhillpositive":

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        d_species_id = getSpeciesByName(model,d).getId()

        ratestring+=f"*{d_species_id}*{s_species_id}^n/({s_species_id}^n + K)"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id

    elif propensity_type == "proportionalhillnegative":

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        d_species_id = getSpeciesByName(model,d).getId()

        ratestring+=f"*{d_species_id}/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id


class HillPositive(Hill):
    def __init__(self, k_forward: float, s1, K: float, n: float):
        """ Hill positive propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*s1^n/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        """
        super(HillPositive, self).__init__(k_forward=k_forward, s1=s1, K=K, n=n)

    def pretty_print(self,**kwargs):
        if 'reaction_direction' in kwargs and kwargs['reaction_direction'] is 'reverse':
            return super().pretty_print(**kwargs)

        return f'k_f({self.s1}) = {self.k_forward}*{self.s1}^{self.n}/({self.s1}^{self.n} + {self.K})'


class HillNegative(Hill):
    def __init__(self, k_forward: float, s1, K: float, n: float):
        """ Hill negative propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*1/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        """
        super(HillNegative, self).__init__(k_forward=k_forward, s1=s1, K=K, n=n)

    def pretty_print(self,**kwargs):
        if 'reaction_direction' in kwargs and kwargs['reaction_direction'] is 'reverse':
            return super().pretty_print(**kwargs)

        return f'k_f({self.s1}) = {self.k_forward}/({self.s1}^{self.n} + {self.K})'


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

    def pretty_print(self,**kwargs):
        if 'reaction_direction' in kwargs and kwargs['reaction_direction'] is 'reverse':
            return super().pretty_print(**kwargs)

        return '*'.join([self.d, super(ProportionalHillPositive, self).pretty_print(**kwargs)])

    def __eq__(self, other):
        return super(ProportionalHillPositive, self).__eq__(other) \
               and self.d is other.d

    @property
    def species(self):
        return [*super(ProportionalHillPositive, self).species, self.d]


class ProportionalHillNegative(HillNegative):
    def __init__(self, k_forward: float, s1, K: float, n: float, d):
        """ proportional Hill negative propensity with the following formula

            p(s1, d; k, K, n) = k*d/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        super(ProportionalHillNegative, self).__init__(k_forward=k_forward, s1=s1, K=K, n=n)
        self.d = d

    def pretty_print(self,**kwargs):
        # TODO no reverse for nonmassaction!!!
        if 'reaction_direction' in kwargs and kwargs['reaction_direction'] is 'reverse':
            return super().pretty_print(**kwargs)

        return '*'.join([self.d, super(ProportionalHillNegative, self).pretty_print(**kwargs)])

    def __eq__(self, other):
        return super(ProportionalHillNegative, self).__eq__(other) \
               and self.d is other.d

    @property
    def species(self):
        return [*super(ProportionalHillNegative, self).species, self.d]