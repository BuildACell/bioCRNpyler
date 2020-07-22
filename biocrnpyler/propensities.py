#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import Set, Union, Dict, List
from numbers import Real
from .parameter import Parameter, ModelParameter, ParameterEntry
from .sbmlutil import getSpeciesByName, _create_ratelaw_parameter
from .species import Species
import libsbml


class Propensity(object):
    def __init__(self):
        self.propensity_dict = {'species': {}, 'parameters': {}}

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
    def _get_parameter_name(rate_coeff: Union[Real, Parameter]) -> str:
        """checks whether a rate_coeff is a Real number, if not this function calls for its string represention

        :param rate_coeff: a number of a parameter object
        :return: either None or the parameter's name
        """
        if not isinstance(rate_coeff, Real):
            rate_coeff_name = str(rate_coeff)
        else:
            rate_coeff_name = 'k'
        return rate_coeff_name

    def _check_parameter(self, parameter, allow_None = False, positive = True):
        """
        A helper function used in setters to set parameters and do type checking
        """
        if isinstance(parameter, Parameter) and (parameter.value > 0 or not positive):
            return parameter
        elif (isinstance(parameter, float) or isinstance(parameter, int)) and (parameter > 0 or not positive):
            return parameter
        elif parameter is None and allow_None:
            return parameter
        else:
            if positive:
                raise ValueError(f"Propensity parameters must be Parameters or floats with positive values. Recieved {parameter}.")
            else:
                raise ValueError(f"Propensity parameters must be Parameters or floats. Recieved {parameter}.")

    def _check_species(self, species, allow_None = False):
        """
        A helper function used in setters to set species and do type checking
        """
        if isinstance(species, Species):
            return species
        elif species is None and allow_None:
            return species
        else:
            raise ValueError(f"Propensity expected a Species: recieved {species}.")


    def pretty_print(self, show_parameters = True, **kwargs):
        txt += self.pretty_print_rate(**kwargs)
        if show_parameters:
            txt += "\n"+self.pretty_print(parameters)
        return txt

    def pretty_print_rate(self, **kwargs):
        raise NotImplementedError("class Propensity is meant to be subclassed!")

    def pretty_print_parameters(self, show_keys = True, **kwargs):
        txt = ""
        for k in self.propensity_dict["parameters"]:
            p = self.propensity_dict["parameters"][k]
            if isinstance(p, Parameter):
                txt += f"\t{k}={p.value}"#p.pretty_print(**kwargs)+"\n"
                if isinstance(p, ModelParameter) and show_keys:
                    txt+=f"\n\t\tfound_key={p.found_key}. search_key={p.search_key}."
                txt+="\n"
            else:
                txt += f"\t{k}={p}\n"
        return txt

    @property
    def is_reversible(self):
        """
        By default, Propensities are assumed to NOT be reversible.
        """
        return False

    @property
    def k_forward(self):
        raise NotImplementedError

    @property
    def k_reverse(self):
        raise NotImplementedError

    def __eq__(self, other):
        if other.__class__ == self.__class__:
            return other.propensity_dict == self.propensity_dict

    @property
    def species(self) -> List:
        """returns the instance variables that are species type"""
        return self.propensity_dict['species'].values()

    def create_kinetic_law(self, reaction, reverse_reaction, stochastic):
        raise NotImplementedError("class Propensity is meant to be subclassed!")

    @classmethod
    def from_dict(cls, propensity_dict):
        merged = propensity_dict['parameters']
        merged.update(propensity_dict['species'])
        return cls(**merged)


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

        self.propensity_dict['parameters']['k_forward'] = k_forward
        if self.is_reversible:
            self.propensity_dict['parameters']['k_reverse'] = k_reverse

    @property
    def k_forward(self):
        if isinstance(self._k_forward, Parameter):
            return self._k_forward.value
        else:
            return self._k_forward

    @k_forward.setter
    def k_forward(self, new_k_forward):
        self._k_forward = self._check_parameter(new_k_forward)
        self.propensity_dict['parameters']['k_forward'] = self.k_forward

    @property
    def k_reverse(self):
        if isinstance(self._k_reverse, Parameter):
            return self._k_reverse.value
        else:
            return self._k_reverse

    @k_reverse.setter
    def k_reverse(self, new_k_reverse):
        self._k_reverse = self._check_parameter(new_k_reverse, allow_None = True)
        self.propensity_dict['parameters']['k_reverse'] = self.k_reverse

    @property
    def is_reversible(self):
        if self.k_reverse is None:
            return False
        else:
            return True

    def pretty_print_rate(self, **kwargs):
        if 'reaction_direction' in kwargs and kwargs['react_direction'] is 'reverse':
            if self.is_reversible:
                return f''
            else:
                return None
        else:
            return f''

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        # Create a kinetic law for the sbml_reaction
        ratelaw = sbml_reaction.createKineticLaw()

        # annotation_dict = {"type": self}
        if 'crn_reaction' in kwargs:
            crn_reaction = kwargs['crn_reaction']
        else:
            raise ValueError('crn_reaction reference is needed for Massaction kinetics!')

        # set up the forward sbml_reaction
        reactant_species = {}
        for w_species in crn_reaction.inputs:
            species_id = getSpeciesByName(model, str(w_species.species)).getId()
            reactant_species[species_id] = w_species

        # get the name of the kinetic rate coefficient
        rate_coeff_name = Propensity._get_parameter_name(self.k_forward)
        param = _create_ratelaw_parameter(ratelaw, rate_coeff_name, self.k_forward)

        # annotation_dict["k"] = rate_coeff

        rate_formula = self._get_rate_formula(rate_coeff_name, stochastic, reactant_species)

        # setup the reverse sbml_reaction
        if crn_reaction.is_reversible:
            sbml_reaction.setReversible(True)
            rate_formula += '-'
            reactant_species = {}
            for w_species in crn_reaction.outputs:
                species_id = getSpeciesByName(model, str(w_species.species)).getId()
                reactant_species[species_id] = w_species

            rate_coeff_name = Propensity._get_parameter_name(self.k_reverse)
            param = _create_ratelaw_parameter(ratelaw, rate_coeff_name, self.k_reverse)

            rate_formula += self._get_rate_formula(rate_coeff_name, stochastic, reactant_species)
        else:
            sbml_reaction.setReversible(False)

        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        ratelaw.setMath(math_ast)

        # if propensity_annotation:
        #     annotation_string = "<PropensityType>"
        #     for k in annotation_dict:
        #         annotation_string += " "+k + "=" + str(annotation_dict[k])
        #     annotation_string += "</PropensityType>"
        #     sbml_reaction.appendAnnotation(annotation_string)

        return ratelaw

    def _get_rate_formula(self, rate_coeff_name, stochastic, reactant_species) -> str:

        # Create Rate-strings for massaction propensities
        ratestring = rate_coeff_name

        for species_id, weighted_species in reactant_species.items():
            if stochastic:
                ratestring += '*'.join(f" ( {species_id} - {i} )" for i in range(weighted_species.stoichiometry))
            else:
                ratestring += f" * {species_id}^{weighted_species.stoichiometry}"

        return ratestring


class Hill(Propensity):
    def __init__(self, k: float, s1: Species, K: float, n: float, d:float):
        Propensity.__init__(self)
        self.s1 = s1
        self.K = K
        self.n = n
        self.d = d

    @property
    def k(self):
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k
    @k.setter
    def k(self, new_k):
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self.k

    @property
    def K(self):
        if isinstance(self._K, Parameter):
            return self._K.value
        else:
            return self._K

    @K.setter
    def K(self, new_K):
        self._K = self._check_parameter(new_K)
        self.propensity_dict['parameters']['K'] = self.K

    @property
    def n(self):
        if isinstance(self._n, Parameter):
            return self._n.value
        else:
            return self._n
    @n.setter
    def n(self, new_n):
        self._n = self._check_parameter(new_n)
        self.propensity_dict['parameters']['n'] = self.n

    @property
    def s1(self):
        return self._s1
    @s1.setter
    def s1(self, new_s1):
        self._s1 = self._check_species(new_s1)
        self.propensity_dict['species']['s1'] = self.s1

    @property
    def d(self):
        return self._d
    @d.setter
    def d(self, new_d):
        self._d = self._check_species(new_d, allow_None = True)
        self.propensity_dict['species']['d'] = self.d

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        raise NotImplementedError("Propensity class Hill is meant to be subclassed: try HillPositive, HillNegative, ProportionalHillPositive, or ProportionalHillNegative.")

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        raise NotImplementedError("Propensity class Hill is meant to be subclassed: try HillPositive, HillNegative, ProportionalHillPositive, or ProportionalHillNegative.")

    
    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        param_n = ratelaw.createParameter()
        param_n.setId("n")
        param_n.setConstant(True)
        param_n.setValue(self.n)
    
        param_K = ratelaw.createParameter()
        param_K.setId("K")
        param_K.setConstant(True)
        param_K.setValue(self.K)

        param_k = ratelaw.createParameter()
        parak_k.setId("k")
        param_K.setConstant(True)
        pram_K.setValue(self.k)

        s = str(self.s1).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()


        # annotation_dict["k"] = self.k_forward
        # annotation_dict["K"] = self.K
        # annotation_dict["n"] = self.n

        #
        #     #Create ratestring for non-massaction propensities
        # if propensity_type == "hillpositive":
        #
        #     ratestring+=f"*{s_species_id}^n/({s_species_id}^n+K)"
        #
        #     annotation_dict["s1"] = s_species_id
        #
        # elif propensity_type == "hillnegative":
        #     s = str(propensity_params['s1']).replace("'", "")
        #     s_species_id = getSpeciesByName(model,s).getId()
        #
        #     ratestring+=f"/({s_species_id}^n+K)"
        #
        #     annotation_dict["s1"] = s_species_id
        #
        # elif propensity_type == "proportionalhillpositive":
        #
        #     s = str(propensity_params['s1']).replace("'", "")
        #     d = str(propensity_params['d']).replace("'", "")
        #     s_species_id = getSpeciesByName(model,s).getId()
        #     d_species_id = getSpeciesByName(model,d).getId()
        #
        #     ratestring+=f"*{d_species_id}*{s_species_id}^n/({s_species_id}^n + K)"
        #
        #     annotation_dict["s1"] = s_species_id
        #     annotation_dict["d"] = d_species_id
        #
        # elif propensity_type == "proportionalhillnegative":
        #
        #     s = str(propensity_params['s1']).replace("'", "")
        #     d = str(propensity_params['d']).replace("'", "")
        #     s_species_id = getSpeciesByName(model,s).getId()
        #     d_species_id = getSpeciesByName(model,d).getId()
        #
        #     ratestring+=f"*{d_species_id}/({s_species_id}^n+K)"
        #
        #     annotation_dict["s1"] = s_species_id
        #     annotation_dict["d"] = d_species_id



class HillPositive(Hill):
    def __init__(self, k: float, s1:Species, K: float, n: float):
        """ Hill positive propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*s1^n/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        """
        Hill.__init__(self = self, k=k, s1=s1, K=K, n=n, d=None)

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f'k_f({self.s1.pretty_print(**kwargs)}) = k {self.s1.pretty_print(**kwargs)}^n / ({self.s1.pretty_print(**kwargs)}^n + K)'
        


class HillNegative(Hill):
    def __init__(self, k: float, s1, K: float, n: float):
        """ Hill negative propensity is a nonlinear propensity with the following formula

            p(s1; k, K, n) = k*1/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        """
        Hill.__init__(self = self, k=k, s1=s1, K=K, n=n, d=None)

    def pretty_print_rate(self, show_parameters = True, **kwargs):
        return f'k_f({self.s1.pretty_print(**kwargs)}) = k / ({self.s1.pretty_print(**kwargs)}^n + K)'


class ProportionalHillPositive(HillPositive):
    def __init__(self, k: float, s1:Species, K: float, n: float, d:Species):
        """ proportional Hill positive propensity with the following formula

            p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        Hill.__init__(self = self, k=k, s1=s1, K=K, n=n, d=d)

    def pretty_print_rate(self, show_parameters = True,  **kwargs):
        return f'k_f({self.s1.pretty_print(**kwargs)}, {self.d.pretty_print(**kwargs)}) = k {self.d.pretty_print(**kwarks)} {self.s1.pretty_print(**kwargs)}^n/({self.s1.pretty_print(**kwargs)}^n + K)'


class ProportionalHillNegative(HillNegative):
    def __init__(self, k_forward: float, s1:Species, K: float, n: float, d:Species):
        """ proportional Hill negative propensity with the following formula

            p(s1, d; k, K, n) = k*d/(s1^n + K)

        :param k_forward: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        Hill.__init__(self = self, k=k, s1=s1, K=K, n=n, d=d)

    def pretty_print(self, show_parameters = True, **kwargs):
        return f'k_f({self.s1.pretty_print(**kwargs)}, {self.d.pretty_print(**kwargs)}) = k {self.d.pretty_print(**kwargs)} /({self.s1.pretty_print(**kwargs)}^{self.n} + K)'
