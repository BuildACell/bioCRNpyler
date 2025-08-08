#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import numbers
from collections import defaultdict
from typing import List, Set, Union

import libsbml  # type: ignore

from .parameter import ModelParameter, Parameter, ParameterEntry
from ..utils.sbmlutil import (_create_global_parameter, _create_local_parameter,
                              getSpeciesByName)
from .species import Species


class Propensity(object):
    def __init__(self):
        self.propensity_dict = {'species': {}, 'parameters': {}}
        self.name = None

    @staticmethod
    def is_valid_propensity(propensity_type) -> bool:
        """checks whether the given propensity_type is valid propensity.

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
        """Returns a set of all subclasses of cls (recursively calculated).

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

    def _create_sbml_parameter(self, parameter_name, sbml_model, ratelaw, rename_dict=None):
        """Creates an sbml parameter for a parameter of the given name.

        if self.propensity_dict["parameter"]["parameter_name"] is a Parameter,
            creates a global parameter "Parameter.name_Parameter.part_id_Parameter.mechanism"
            where part_id and mechanism can be empty (but _ will always be incldued for uniqueness).
        if self.propensity_dict["parameter"]["parameter_name"] is a Number,
            creates a local parameter "parameter_name".
        rename_dict allows for param.name to be changed to rename_dict[param.name]
        """
        p = self.propensity_dict["parameters"][parameter_name]
        if isinstance(p, ParameterEntry):
            v = p.value
            p_unit = p.unit
            if p_unit == "":
                p_unit = None
            m = p.parameter_key.mechanism
            if m is None:
                m = ""
            pid = p.parameter_key.part_id
            if pid is None:
                pid = ""

            if rename_dict is None or p.parameter_name not in rename_dict:
                sbml_name = p.parameter_name+"_"+pid+"_"+m
            else:
                sbml_name = rename_dict[p.parameter_name]+"_"+pid+"_"+m

            return _create_global_parameter(sbml_model, sbml_name, v, p_unit)

        elif isinstance(p, int) or isinstance(p, float):
            v = p
            if rename_dict is None or parameter_name not in rename_dict:
                sbml_name = parameter_name
            else:
                sbml_name = rename_dict[parameter_name]

            return _create_local_parameter(ratelaw, sbml_name, v)

        else:
            raise TypeError(
                f"Invalid item in propensity_diction['parameter']: {p}. Only numbers of ParameterEntries accepted.")

    def _check_parameter(self, parameter, allow_None=False, positive=True):
        """A helper function used in setters to set parameters and do type checking."""
        if isinstance(parameter, Parameter) and (parameter.value > 0 or not positive):
            return parameter
        elif isinstance(parameter, numbers.Real) and (parameter > 0 or not positive):
            return parameter
        elif parameter is None and allow_None:
            return parameter
        else:
            if positive:
                raise ValueError(
                    f"Propensity parameters must be Parameters or floats with positive values. Recieved {type(parameter)}.")
            else:
                raise ValueError(
                    f"Propensity parameters must be Parameters or floats. Recieved {type(parameter)}.")

    def _check_species(self, species, allow_None=False):
        """A helper function used in setters to set species and do type checking."""
        if isinstance(species, Species):
            return species
        elif species is None and allow_None:
            return species
        else:
            raise TypeError(
                f"Propensity expected a Species: received {type(species)}.")

    def pretty_print(self, show_parameters=True, **kwargs):
        txt = self.pretty_print_rate(**kwargs)
        if show_parameters:
            txt += "\n"+self.pretty_print_parameters(**kwargs)
        return txt

    def pretty_print_rate(self, **kwargs):
        raise NotImplementedError(
            "class Propensity is meant to be subclassed!")

    def pretty_print_parameters(self, show_keys=True, **kwargs):
        txt = ""
        for k in self.propensity_dict["parameters"]:
            p = self.propensity_dict["parameters"][k]
            if isinstance(p, Parameter):
                txt += f"  {k}={p.value}"  # p.pretty_print(**kwargs)+"\n"
                if isinstance(p, ModelParameter) and show_keys:
                    txt += f"\n  found_key=(mech={p.found_key.mechanism}, partid={p.found_key.part_id}, name={p.found_key.name}).\n  search_key=(mech={p.search_key.mechanism}, partid={p.search_key.part_id}, name={p.search_key.name})."
                txt += "\n"
            elif p is not None:
                txt += f"  {k}={p}\n"
        return txt

    @property
    def is_reversible(self):
        """By default, Propensities are assumed to NOT be reversible."""
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
        """returns the instance variables that are species type."""
        return list(self.propensity_dict['species'].values())

    def create_kinetic_law(self, reaction, reverse_reaction, stochastic, **kwargs):
        raise NotImplementedError(
            "class Propensity is meant to be subclassed!")

    @classmethod
    def from_dict(cls, propensity_dict):
        merged = propensity_dict['parameters']
        merged.update(propensity_dict['species'])
        return cls(**merged)

    def _create_annotation(self, model, propensity_dict_in_sbml, **kwargs):
        """Create simulator specific annotations to write to Kinetic laws or any other
        part of the SBML model object.

        Annotations are used to take advantage of a simulator specific need/feature.
        """
        annotation_string = ''
        # Add your own simulator specific annotations here

        # For the bioscrape simulator:
        # Check if `for_bioscrape` keyword argument has been passed in **kwargs
        if 'for_bioscrape' in kwargs:
            for_bioscrape = kwargs.get('for_bioscrape')
        else:
            for_bioscrape = False

        if for_bioscrape:
            annotation_string = self._create_bioscrape_annotation(
                propensity_dict_in_sbml)

        return annotation_string

    def _create_bioscrape_annotation(self, propensity_dict_in_sbml):
        """Propensity Annotations are Used to take advantage of Bioscrape Propensity types
        for faster simulation."""
        annotation_dict = defaultdict()
        for param_name, param_value in propensity_dict_in_sbml['parameters'].items():
            annotation_dict[param_name] = param_value

        for species_name, species in propensity_dict_in_sbml['species'].items():
            annotation_dict[species_name] = species

        annotation_dict["type"] = self.name

        annotation_string = "<PropensityType>"
        for k in annotation_dict:
            annotation_string += " " + str(k) + "=" + str(annotation_dict[k])
        annotation_string += "</PropensityType>"

        # replace strings to match with bioscrape naming convention
        annotation_string = annotation_string.replace('k_forward', 'k', 1)
        # Bioscrape doesn't have the concept of a reversible reaction - so for both the reverse and forward cases
        # we just make the annotation parameter be called 'k'
        annotation_string = annotation_string.replace('k_reverse', 'k', 1)
        return annotation_string

    def _translate_propensity_dict_to_sbml(self, model, ratelaw):
        # get copy of the propensity_dict and fill with sbml names
        propensity_dict_in_sbml = copy.deepcopy(self.propensity_dict)
        for param_name in propensity_dict_in_sbml['parameters'].keys():
            parameter_in_sbml = self._create_sbml_parameter(
                param_name, model, ratelaw)
            propensity_dict_in_sbml['parameters'][param_name] = parameter_in_sbml.getId(
            )

        for species_name, species in propensity_dict_in_sbml['species'].items():
            propensity_dict_in_sbml['species'][species_name] = str(species)

        return propensity_dict_in_sbml


class GeneralPropensity(Propensity):
    def __init__(self, propensity_function: str, propensity_species: List[Species], propensity_parameters: List[ParameterEntry]):
        """A class to define a general propensity.

        :param propensity_function: valid propensity formula defined as a string
        :param propensity_species: list of species that are part of the propensity_function
        :param propensity_parameters: list of parameters that are part of the propensity_function
        """
        super(GeneralPropensity, self).__init__()
        self.propensity_function = propensity_function

        if len(propensity_species) > 0 and not all(isinstance(s, Species) for s in propensity_species):
            raise TypeError('propensity_species must be a list of Species!')

        if len(propensity_parameters) > 0 and not all(isinstance(s, ParameterEntry) for s in propensity_parameters):
            raise TypeError(
                'propensity_parameter must be a list of ParameterEntry!')

        for species in propensity_species:
            if str(species) not in self.propensity_function:
                raise ValueError(
                    f'species: {species} must be part of the formula: {self.propensity_function}')

            self.propensity_dict['species'].update({str(species): species})

        for parameter in propensity_parameters:
            if parameter.parameter_name not in self.propensity_function:
                raise ValueError(
                    f'species: {parameter.parameter_name} must be part of the formula: {self.propensity_function}')

            self.propensity_dict['parameters'].update(
                {parameter.parameter_name: parameter.value})

        self.name = 'general'

    def pretty_print_rate(self, **kwargs):
        return self.propensity_function

    def create_kinetic_law(self, model, sbml_reaction, **kwargs):
        """Creates KineticLaw object for SBML using the propensity_function string."""
        ratelaw = sbml_reaction.createKineticLaw()

        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw)

        # replacing the species defined in CRN with valid SBML names
        for species_in_crn, species_in_sbml in propensity_dict_in_sbml['species'].items():
            self.propensity_function = self.propensity_function.replace(
                species_in_crn, species_in_sbml)

        # replacing the parameters defined in CRN with valid SBML names
        for parameter_in_crn, parameter_in_sbml in propensity_dict_in_sbml['parameters'].items():
            self.propensity_function = self.propensity_function.replace(
                parameter_in_crn, parameter_in_sbml)

        math_ast = libsbml.parseL3Formula(self.propensity_function)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. Check the propensity functions of reactions.")
        return ratelaw


class MassAction(Propensity):
    def __init__(self, k_forward: Union[float, ParameterEntry], k_reverse: Union[float, ParameterEntry] = None):
        super(MassAction, self).__init__()
        self.k_forward = k_forward
        self.k_reverse = k_reverse
        self.name = 'massaction'

    @property
    def k_forward(self):
        if isinstance(self._k_forward, Parameter):
            return self._k_forward.value
        else:
            return self._k_forward

    @k_forward.setter
    def k_forward(self, new_k_forward):
        self._k_forward = self._check_parameter(new_k_forward)
        self.propensity_dict['parameters']['k_forward'] = self._k_forward

    @property
    def k_reverse(self):
        if isinstance(self._k_reverse, Parameter):
            return self._k_reverse.value
        else:
            return self._k_reverse

    @k_reverse.setter
    def k_reverse(self, new_k_reverse):
        self._k_reverse = self._check_parameter(new_k_reverse, allow_None=True)
        if self._k_reverse is not None:
            self.propensity_dict['parameters']['k_reverse'] = self._k_reverse

    @property
    def is_reversible(self):
        if self.k_reverse is None:
            return False
        else:
            return True

    def pretty_print_rate(self, **kwargs):
        crn_reaction = kwargs["reaction"]
        reactant_species = {}
        for w_species in crn_reaction.inputs:
            reactant_species[str(w_species.species)] = w_species
        txt = " Kf=" + \
            self._get_rate_formula(
                "k_forward", kwargs["stochastic"], reactant_species)
        if self.is_reversible:
            reactant_species = {}
            for w_species in crn_reaction.outputs:
                reactant_species[str(w_species.species)] = w_species
            txt += "\n Kr=" + \
                self._get_rate_formula(
                    "k_reverse", kwargs["stochastic"], reactant_species)
        return txt

    def create_kinetic_law(self, model, sbml_reaction, stochastic, reverse_reaction=False, **kwargs):

        if 'crn_reaction' in kwargs:
            crn_reaction = kwargs['crn_reaction']
        else:
            raise ValueError(
                'crn_reaction reference is needed for Massaction kinetics!')

        # create a kinetic law for the sbml_reaction
        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw)

        # set up the forward sbml_reaction
        if not reverse_reaction:
            reactant_species = {}
            for w_species in crn_reaction.inputs:
                species_id = str(w_species.species)
                reactant_species[species_id] = w_species
            param = propensity_dict_in_sbml['parameters']['k_forward']
            # remove the other parameter from the propensities
            propensity_dict_in_sbml['parameters'].pop('k_reverse', None)
            # if k_reverse is a local parameter, remove it
            ratelaw.removeLocalParameter("k_reverse")
        # set up a reverse reaction
        elif reverse_reaction:
            reactant_species = {}
            for w_species in crn_reaction.outputs:
                species_id = str(w_species.species)
                reactant_species[species_id] = w_species
            param = propensity_dict_in_sbml['parameters']['k_reverse']
            # remove the other parameter from the propensities
            propensity_dict_in_sbml['parameters'].pop('k_forward', None)
            # if k_forward is a local parameter, remove it
            ratelaw.removeLocalParameter("k_forward")

        rate_formula = self._get_rate_formula(
            param, stochastic, reactant_species)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. Check the propensity functions of reactions.")
        annotation_string = self._create_annotation(
            model, propensity_dict_in_sbml=propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        return ratelaw

    def _get_rate_formula(self, rate_coeff_name, stochastic, reactant_species) -> str:

        # Create Rate-strings for massaction propensities
        ratestring = rate_coeff_name

        for species_id, weighted_species in reactant_species.items():
            if stochastic:
                ratestring += '*'
                ratestring += f"{species_id}"
                ratestring += '*'
                ratestring += '*'.join(f" ( {species_id} - {i} )" for i in range(
                    1, weighted_species.stoichiometry))

                # Remove trailing *
                if ratestring[len(ratestring)-1] == "*":
                    ratestring = ratestring[:-1]
            else:
                if weighted_species.stoichiometry > 1:
                    ratestring += f" * {species_id}^{weighted_species.stoichiometry}"
                else:
                    ratestring += f" * {species_id}"
        return ratestring


class Hill(Propensity):
    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        Propensity.__init__(self)
        self.k = k
        self.s1 = s1
        self.K = K
        self.n = n
        if d is not None:
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
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def K(self):
        if isinstance(self._K, Parameter):
            return self._K.value
        else:
            return self._K

    @K.setter
    def K(self, new_K):
        self._K = self._check_parameter(new_K)
        self.propensity_dict['parameters']['K'] = self._K

    @property
    def n(self):
        if isinstance(self._n, Parameter):
            return self._n.value
        else:
            return self._n

    @n.setter
    def n(self, new_n):
        self._n = self._check_parameter(new_n)
        self.propensity_dict['parameters']['n'] = self._n

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
        self._d = self._check_species(new_d, allow_None=True)
        self.propensity_dict['species']['d'] = self._d

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        raise NotImplementedError(
            "Propensity class Hill is meant to be subclassed: try HillPositive, HillNegative, ProportionalHillPositive, or ProportionalHillNegative.")

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """This code is reused in all Hill Propensity subclasses."""
        if 'reverse_reaction' in kwargs and kwargs['reverse_reaction'] is True:
            raise ValueError(
                'reverse reactions cannot exist for Hill type Propensities!')

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw)

        rate_formula = self._get_rate_formula(
            propensity_dict=propensity_dict_in_sbml)
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(
            model, propensity_dict_in_sbml, **kwargs)
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. Check the propensity functions of reactions.")

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        raise NotImplementedError(
            'Hill does not have a rate formula! Check out the subclasses.')


class HillPositive(Hill):
    def __init__(self, k: float, s1: Species, K: float, n: float):
        """ Hill positive propensity is a nonlinear propensity with the following formula.

            p(s1; k, K, n) = k*s1^n/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        """
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=None)
        self.name = 'hillpositive'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        return f' Kf = k {self.s1.pretty_print(**kwargs)}^n / ( K^n + {self.s1.pretty_print(**kwargs)}^n )'

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        rate_formula = f"{k}*{s1}^{n} / ( {K}^{n} + {s1}^{n} )"
        return rate_formula


class HillNegative(Hill):
    def __init__(self, k: float, s1: Species, K: float, n: float):
        """ Hill negative propensity is a nonlinear propensity with the following formula.

            p(s1; k, K, n) = k*1/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        """
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=None)
        self.name = 'hillnegative'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        return f' Kf = k / ( 1 + ({self.s1.pretty_print(**kwargs)}/K)^n )'

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        rate_formula = f"{k} / ( 1 + ({s1}/{K})^{n} )"
        return rate_formula


class ProportionalHillPositive(HillPositive):
    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        """ proportional Hill positive propensity with the following formula.

            p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=d)
        self.name = 'proportionalhillpositive'

    def pretty_print_rate(self, show_parameters=True,  **kwargs):
        return f' Kf = k {self.d.pretty_print(**kwargs)} {self.s1.pretty_print(**kwargs)}^n / ( K^n + {self.s1.pretty_print(**kwargs)}^n )'

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        d = propensity_dict['species']['d']
        return f"{k}*{d}*{s1}^{n} / ( {K}^{n} + {s1}^{n} )"


class ProportionalHillNegative(HillNegative):
    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        """ proportional Hill negative propensity with the following formula.

            p(s1, d; k, K, n) = k*d/(s1^n + K)

        :param k: rate constant (float)
        :param s1: species (chemical_reaction_network.species)
        :param K: dissociation constant (float)
        :param n: cooperativity (float)
        :param d: species (chemical_reaction_network.species)
        """
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=d)
        self.name = 'proportionalhillnegative'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        return f' Kf = k {self.d.pretty_print(**kwargs)} / ( 1 + ({self.s1.pretty_print(**kwargs)}/K)^{self.n} )'

    def _get_rate_formula(self, propensity_dict):
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        d = propensity_dict['species']['d']
        return f"{k}*{d} / ( 1 + ({s1}/{K})^{n} )"
