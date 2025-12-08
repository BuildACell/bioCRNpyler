#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import numbers
from collections import defaultdict
from typing import List, Set, Union

import libsbml  # type: ignore

from ..utils.sbmlutil import _create_global_parameter, _create_local_parameter
from .parameter import ModelParameter, Parameter, ParameterEntry
from .species import Species


class Propensity(object):
    """Base class for reaction propensity functions in BioCRNpyler.

    Propensities define the rate laws for chemical reactions in a CRN.
    Different propensity types implement different kinetic models such as
    mass action, Hill functions, and custom formulas. Propensities can be
    deterministic (ODE) or stochastic (Gillespie).

    Attributes
    ----------
    propensity_dict : dict
        Dictionary with 'species' and 'parameters' keys storing the
        species and parameters used in the propensity function.
    name : str or None
        Name identifier for the propensity type.

    See Also
    --------
    MassAction : Mass action kinetics propensity.
    GeneralPropensity : Custom formula propensity.
    Hill : Base class for Hill-type propensities.

    Notes
    -----
    This is an abstract base class that should be subclassed to implement
    specific propensity types. Subclasses must implement:

    - `create_kinetic_law`: Generate SBML kinetic law
    - `pretty_print_rate`: Human-readable rate formula

    The `propensity_dict` structure:

    - 'species': {<name>: Species object, ...}
    - 'parameters': {<name>: Parameter or number, ...}

    """

    def __init__(self):
        """Initialize a Propensity object.

        Creates an empty propensity dictionary with 'species' and
        'parameters' keys.

        """
        self.propensity_dict = {'species': {}, 'parameters': {}}
        self.name = None

    @staticmethod
    def is_valid_propensity(propensity_type) -> bool:
        """Check if an object is a valid Propensity subclass instance.

        Recursively checks all subclasses of Propensity to determine if
        the given object is a valid propensity type.

        Parameters
        ----------
        propensity_type : object
            Object to check for Propensity validity.

        Returns
        -------
        bool
            True if `propensity_type` is an instance of Propensity or any
            of its subclasses, False otherwise.

        Examples
        --------
        >>> prop = bcp.MassAction(k_forward=100.0)
        >>> bcp.Propensity.is_valid_propensity(prop)
        True
        >>> bcp.Propensity.is_valid_propensity("not a propensity")
        False

        """
        for propensity in Propensity.get_available_propensities():
            if isinstance(propensity_type, propensity):
                return True
        return False

    @staticmethod
    def _all_subclasses(cls):
        """Recursively find all subclasses of a class.

        Parameters
        ----------
        cls : type
            A class to find subclasses of.

        Returns
        -------
        set
            Set of all subclasses of `cls`, including nested subclasses.

        Notes
        -----
        Source: https://stackoverflow.com/questions/3862310/
        how-to-find-all-the-subclasses-of-a-class-given-its-name

        """
        return set(cls.__subclasses__()).union(
            [
                s
                for c in cls.__subclasses__()
                for s in Propensity._all_subclasses(c)
            ]
        )

    @staticmethod
    def get_available_propensities() -> Set:
        """Get all available propensity subclasses.

        Returns
        -------
        set
            Set of all Propensity subclass types available in BioCRNpyler.

        Examples
        --------
        >>> propensities = bcp.Propensity.get_available_propensities()
        >>> bcp.MassAction in propensities
        True

        """
        return Propensity._all_subclasses(Propensity)

    def _create_sbml_parameter(
        self, parameter_name, sbml_model, ratelaw, rename_dict=None
    ):
        """Create an SBML parameter for the propensity.

        Creates either a global or local SBML parameter depending on
        whether the parameter is a ParameterEntry or a number.

        Parameters
        ----------
        parameter_name : str
            Name of the parameter in propensity_dict['parameters'].
        sbml_model : libsbml.Model
            SBML model object to add global parameters to.
        ratelaw : libsbml.KineticLaw
            SBML kinetic law object to add local parameters to.
        rename_dict : dict, optional
            Dictionary to rename parameters. Maps original parameter names
            to new names.

        Returns
        -------
        libsbml.Parameter or libsbml.LocalParameter
            Created SBML parameter object.

        Raises
        ------
        TypeError
            If parameter is not a ParameterEntry, int, or float.

        Notes
        -----
        If parameter is a ParameterEntry, creates a global parameter named
        '<parameter_name>_<part_id>_<mechanism>'. If parameter is a number,
        creates a local parameter with the given name.

        """
        p = self.propensity_dict['parameters'][parameter_name]
        if isinstance(p, ParameterEntry):
            v = p.value
            p_unit = p.unit
            if p_unit == '':
                p_unit = None
            m = p.parameter_key.mechanism
            if m is None:
                m = ''
            pid = p.parameter_key.part_id
            if pid is None:
                pid = ''

            if rename_dict is None or p.parameter_name not in rename_dict:
                sbml_name = p.parameter_name + '_' + pid + '_' + m
            else:
                sbml_name = (
                    rename_dict[p.parameter_name] + '_' + pid + '_' + m
                )

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
                f"Invalid item in propensity_diction['parameter']: {p}. "
                "Only numbers of ParameterEntries accepted."
            )

    def _check_parameter(self, parameter, allow_None=False, positive=True):
        """Validate parameter type and value.

        Parameters
        ----------
        parameter : Parameter, float, or None
            Parameter to validate.
        allow_None : bool, default=False
            If True, allows None as a valid value.
        positive : bool, default=True
            If True, requires parameter value to be positive.

        Returns
        -------
        Parameter, float, or None
            The validated parameter.

        Raises
        ------
        ValueError
            If parameter is invalid type or has invalid value.

        """
        if isinstance(parameter, Parameter) and (
            parameter.value > 0 or not positive
        ):
            return parameter
        elif isinstance(parameter, numbers.Real) and (
            parameter > 0 or not positive
        ):
            return parameter
        elif parameter is None and allow_None:
            return parameter
        else:
            if positive:
                raise ValueError(
                    f"Propensity parameters must be Parameters or floats "
                    f"with positive values. Recieved {type(parameter)}."
                )
            else:
                raise ValueError(
                    "Propensity parameters must be Parameters or floats. "
                    f"Recieved {type(parameter)}."
                )

    def _check_species(self, species, allow_None=False):
        """Validate species type.

        Parameters
        ----------
        species : Species or None
            Species to validate.
        allow_None : bool, default=False
            If True, allows None as a valid value.

        Returns
        -------
        Species or None
            The validated species.

        Raises
        ------
        TypeError
            If species is not a Species object and not None when allowed.

        """
        if isinstance(species, Species):
            return species
        elif species is None and allow_None:
            return species
        else:
            raise TypeError(
                f"Propensity expected a Species: received {type(species)}."
            )

    def pretty_print(self, show_parameters=True, **kwargs):
        """Generate human-readable string representation of propensity.

        Parameters
        ----------
        show_parameters : bool, default=True
            If True, includes parameter values in output.
        **kwargs
            Additional keyword arguments passed to formatting methods.

        Returns
        -------
        str
            Formatted string showing rate formula and optionally
            parameters.

        """
        txt = self.pretty_print_rate(**kwargs)
        if show_parameters:
            txt += '\n' + self.pretty_print_parameters(**kwargs)
        return txt

    def pretty_print_rate(self, **kwargs):
        """Generate human-readable rate formula string.

        Parameters
        ----------
        **kwargs
            Formatting keyword arguments (e.g., reaction, stochastic).

        Returns
        -------
        str
            Formatted rate formula string.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses.

        """
        raise NotImplementedError(
            "class Propensity is meant to be subclassed!"
        )

    def pretty_print_parameters(self, show_keys=True, **kwargs):
        """Generate formatted string of all propensity parameters.

        Parameters
        ----------
        show_keys : bool, default=True
            If True, shows search and found keys for ModelParameter
            objects (useful for debugging parameter lookup).
        **kwargs
            Additional formatting keyword arguments.

        Returns
        -------
        str
            Formatted string listing all parameters and their values.

        """
        txt = ''
        for k in self.propensity_dict['parameters']:
            p = self.propensity_dict['parameters'][k]
            if isinstance(p, Parameter):
                txt += f"  {k}={p.value}"  # p.pretty_print(**kwargs)+"\n"
                if isinstance(p, ModelParameter) and show_keys:
                    txt += (
                        f"\n  found_key=(mech={p.found_key.mechanism}, "
                        f"partid={p.found_key.part_id}, "
                        f"name={p.found_key.name})."
                        f"\n  search_key=(mech={p.search_key.mechanism}, "
                        f"partid={p.search_key.part_id}, "
                        f"name={p.search_key.name})."
                    )
                txt += '\n'
            elif p is not None:
                txt += f"  {k}={p}\n"
        return txt

    @property
    def is_reversible(self):
        """bool: Whether the propensity represents a reversible reaction.

        Default is False. Subclasses override this for reversible kinetics.

        """
        return False

    @property
    def k_forward(self):
        """Float : Forward rate constant.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses that use rate constants.

        """
        raise NotImplementedError

    @property
    def k_reverse(self):
        """Float or None: Reverse rate constant for reversible reactions.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses that use rate constants.

        """
        raise NotImplementedError

    def __eq__(self, other):
        """Test equality between propensities.

        Parameters
        ----------
        other : Propensity
            Other propensity to compare with.

        Returns
        -------
        bool
            True if propensities have the same class and propensity_dict.

        """
        if other.__class__ == self.__class__:
            return other.propensity_dict == self.propensity_dict

    @property
    def species(self) -> List:
        """List of Species : All species used in the propensity function."""
        return list(self.propensity_dict['species'].values())

    def create_kinetic_law(
        self, reaction, reverse_reaction, stochastic, **kwargs
    ):
        """Create SBML kinetic law for a reaction.

        Parameters
        ----------
        reaction : libsbml.Reaction
            SBML reaction object to add kinetic law to.
        reverse_reaction : bool
            If True, creates kinetic law for reverse direction.
        stochastic : bool
            If True, uses stochastic propensity formulation.
        **kwargs
            Additional arguments (e.g., crn_reaction, model).

        Returns
        -------
        libsbml.KineticLaw
            Created SBML kinetic law object.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses.

        """
        raise NotImplementedError(
            "class Propensity is meant to be subclassed!"
        )

    @classmethod
    def from_dict(cls, propensity_dict):
        """Create a propensity from a dictionary.

        Parameters
        ----------
        propensity_dict : dict
            Dictionary with 'parameters' and 'species' keys containing
            parameter and species values.

        Returns
        -------
        Propensity
            New instance of the propensity class.

        """
        merged = propensity_dict['parameters']
        merged.update(propensity_dict['species'])
        return cls(**merged)

    def _create_annotation(self, model, propensity_dict_in_sbml, **kwargs):
        """Create simulator-specific annotations for SBML model.

        Annotations enable simulator-specific features and optimizations.
        Currently supports bioscrape annotations.

        Parameters
        ----------
        model : libsbml.Model
            SBML model object.
        propensity_dict_in_sbml : dict
            Propensity dictionary with SBML identifiers.
        **kwargs
            Additional arguments. If 'for_bioscrape' is True, creates
            bioscrape annotations.

        Returns
        -------
        str
            XML annotation string to append to SBML reaction.

        """
        annotation_string = ''
        # Add your own simulator specific annotations here

        # For the bioscrape simulator:
        # Check if `for_bioscrape` keyword argument was passed in **kwargs
        if 'for_bioscrape' in kwargs:
            for_bioscrape = kwargs.get('for_bioscrape')
        else:
            for_bioscrape = False

        if for_bioscrape:
            annotation_string = self._create_bioscrape_annotation(
                propensity_dict_in_sbml
            )

        return annotation_string

    def _create_bioscrape_annotation(self, propensity_dict_in_sbml):
        """Generate bioscrape-specific propensity type annotation.

        Creates XML annotation that bioscrape uses to optimize simulation
        by identifying propensity types.

        Parameters
        ----------
        propensity_dict_in_sbml : dict
            Propensity dictionary with SBML identifiers for species and
            parameters.

        Returns
        -------
        str
            XML annotation string with bioscrape PropensityType tag.

        Notes
        -----
        Bioscrape uses propensity type annotations for faster simulation.
        The annotation includes parameter and species names with their
        SBML identifiers, plus the propensity type name.

        """
        annotation_dict = defaultdict()
        for param_name, param_value in propensity_dict_in_sbml[
            'parameters'
        ].items():
            annotation_dict[param_name] = param_value

        for species_name, species in propensity_dict_in_sbml[
            'species'
        ].items():
            annotation_dict[species_name] = species

        annotation_dict['type'] = self.name

        annotation_string = '<PropensityType>'
        for k in annotation_dict:
            annotation_string += ' ' + str(k) + '=' + str(annotation_dict[k])
        annotation_string += '</PropensityType>'

        # replace strings to match with bioscrape naming convention
        annotation_string = annotation_string.replace('k_forward', 'k', 1)
        # Bioscrape doesn't have the concept of a reversible reaction
        # - so for both the reverse and forward cases.  We just make
        # the annotation parameter be called 'k'.
        annotation_string = annotation_string.replace('k_reverse', 'k', 1)
        return annotation_string

    def _translate_propensity_dict_to_sbml(self, model, ratelaw):
        """Translate internal propensity representation to SBML format.

        Converts propensity_dict with Species objects and Parameters to
        SBML identifiers that can be used in SBML rate formulas.

        Parameters
        ----------
        model : libsbml.Model
            SBML model object to add global parameters to.
        ratelaw : libsbml.KineticLaw
            SBML kinetic law object to add local parameters to.

        Returns
        -------
        dict
            Copy of propensity_dict with species and parameters replaced
            by their SBML identifier strings.

        """
        # get copy of the propensity_dict and fill with sbml names
        propensity_dict_in_sbml = copy.deepcopy(self.propensity_dict)
        for param_name in propensity_dict_in_sbml['parameters'].keys():
            parameter_in_sbml = self._create_sbml_parameter(
                param_name, model, ratelaw
            )
            propensity_dict_in_sbml['parameters'][param_name] = (
                parameter_in_sbml.getId()
            )

        for species_name, species in propensity_dict_in_sbml[
            'species'
        ].items():
            propensity_dict_in_sbml['species'][species_name] = str(species)

        return propensity_dict_in_sbml


class GeneralPropensity(Propensity):
    """Propensity with user-defined formula string.

    A `GeneralPropensity` allows specification of arbitrary kinetic rate
    laws using a formula string. The formula can reference species and
    parameters that are validated and tracked.

    Parameters
    ----------
    propensity_function : str
        Valid mathematical formula as a string (e.g., 'k*S1*S2'). Must
        contain all referenced species and parameters.
    propensity_species : list of Species
        List of Species objects used in the formula. Each species must
        appear in the propensity_function string.
    propensity_parameters : list of ParameterEntry
        List of ParameterEntry objects used in the formula. Each
        parameter name must appear in the propensity_function string.

    Attributes
    ----------
    propensity_function : str
        The mathematical formula defining the rate law.
    name : str
        Set to 'general' for this propensity type.

    Raises
    ------
    TypeError
        If propensity_species or propensity_parameters contain invalid
        types.
    ValueError
        If species or parameters in lists do not appear in the formula.

    See Also
    --------
    MassAction : Mass action kinetics propensity.
    Hill : Hill-type propensities.

    Notes
    -----
    The propensity_function string must be a valid mathematical expression
    that can be parsed by libsbml.parseL3Formula(). It can include:

    - Arithmetic operators: `+, -, *, /, ^`
    - Mathematical functions: exp, log, sin, cos, etc.
    - Species names (as strings matching their representation)
    - Parameter names (matching ParameterEntry.parameter_name)

    Examples
    --------
    Create a custom Michaelis-Menten propensity:

    >>> S = bcp.Species('S')
    >>> E = bcp.Species('E')
    >>> kcat = bcp.ParameterEntry('kcat', 0.1)
    >>> Km = bcp.ParameterEntry('Km', 10.0)
    >>> prop = bcp.GeneralPropensity(
    ...     propensity_function='kcat*E*S/(Km + S)',
    ...     propensity_species=[S, E],
    ...     propensity_parameters=[kcat, Km]
    ... )

    Create a custom regulatory function:

    >>> X = bcp.Species('X')
    >>> Y = bcp.Species('Y')
    >>> k1 = bcp.ParameterEntry('k1', 1.0)
    >>> k2 = bcp.ParameterEntry('k2', 0.5)
    >>> prop = bcp.GeneralPropensity(
    ...     propensity_function='k1*X + k2*Y^2',
    ...     propensity_species=[X, Y],
    ...     propensity_parameters=[k1, k2]
    ... )

    """

    def __init__(
        self,
        propensity_function: str,
        propensity_species: List[Species],
        propensity_parameters: List[ParameterEntry],
    ):
        super(GeneralPropensity, self).__init__()
        self.propensity_function = propensity_function

        if len(propensity_species) > 0 and not all(
            isinstance(s, Species) for s in propensity_species
        ):
            raise TypeError("propensity_species must be a list of Species!")

        if len(propensity_parameters) > 0 and not all(
            isinstance(s, ParameterEntry) for s in propensity_parameters
        ):
            raise TypeError(
                "propensity_parameter must be a list of ParameterEntry!"
            )

        for species in propensity_species:
            if str(species) not in self.propensity_function:
                raise ValueError(
                    f"species: {species} must be part of the formula: "
                    f"{self.propensity_function}"
                )

            self.propensity_dict['species'].update({str(species): species})

        for parameter in propensity_parameters:
            if parameter.parameter_name not in self.propensity_function:
                raise ValueError(
                    f"species: {parameter.parameter_name} must be part of "
                    f"the formula: {self.propensity_function}"
                )

            self.propensity_dict['parameters'].update(
                {parameter.parameter_name: parameter.value}
            )

        self.name = 'general'

    def pretty_print_rate(self, **kwargs):
        """Return the propensity function formula string.

        Returns
        -------
        str
            The propensity_function formula.

        """
        return self.propensity_function

    def create_kinetic_law(self, model, sbml_reaction, **kwargs):
        """Create SBML kinetic law using the propensity_function string.

        Translates species and parameter names to SBML identifiers and
        creates the kinetic law from the formula string.

        Parameters
        ----------
        model : libsbml.Model
            SBML model object.
        sbml_reaction : libsbml.Reaction
            SBML reaction to add kinetic law to.
        **kwargs
            Additional keyword arguments (unused for GeneralPropensity).

        Returns
        -------
        libsbml.KineticLaw
            Created SBML kinetic law object.

        Raises
        ------
        ValueError
            If the propensity_function cannot be parsed as valid SBML
            formula.

        """
        ratelaw = sbml_reaction.createKineticLaw()

        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw
        )

        # replacing the species defined in CRN with valid SBML names
        for species_in_crn, species_in_sbml in propensity_dict_in_sbml[
            'species'
        ].items():
            self.propensity_function = self.propensity_function.replace(
                species_in_crn, species_in_sbml
            )

        # replacing the parameters defined in CRN with valid SBML names
        for parameter_in_crn, parameter_in_sbml in propensity_dict_in_sbml[
            'parameters'
        ].items():
            self.propensity_function = self.propensity_function.replace(
                parameter_in_crn, parameter_in_sbml
            )

        math_ast = libsbml.parseL3Formula(self.propensity_function)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. "
                "Check the propensity functions of reactions."
            )
        return ratelaw


class MassAction(Propensity):
    r"""Mass action kinetics propensity.

    Implements mass action rate laws for chemical reactions. Supports both
    irreversible and reversible reactions. Propensities are computed
    differently for deterministic (ODE) and stochastic (Gillespie)
    simulations.

    Parameters
    ----------
    k_forward : float or ParameterEntry
        Forward reaction rate constant. Must be positive.
    k_reverse : float, ParameterEntry, or None, optional
        Reverse reaction rate constant. If None, reaction is irreversible.
        If provided, must be positive.

    Attributes
    ----------
    name : str
        Set to 'massaction' for this propensity type.

    See Also
    --------
    GeneralPropensity : Custom formula propensity.
    Hill : Hill-type propensities.

    Notes
    -----
    Deterministic (ODE) propensity: For reaction A + B --> C with
    rate constant $k$:

        $$ 'rate' = k [A] [B] $$

    where [A] and [B] are the concentrations of A and B.

    Stochastic (Gillespie) propensity: For reaction A + B --> C with
    rate constant $k$:
    $$
        'propensity' &= k \cdot A \cdot (B-1) &\quad& \text{if $A = B$} \\
        'propensity' &= k \cdot A \cdot B     &\quad& \text{otherwise}
    $$
    where $A$ and $B$ represent molecular counts.

    The stochastic formulation accounts for combinatorics of molecule
    selection. For stoichiometric coefficient $n > 1$:

        $$ 'factor' = S ... (S-1) ... (S-n+1) $$

    If `k_reverse` is provided, the reaction is reversible:

        $$ A + B <--> C $$

    Two kinetic laws are created: one for forward, one for reverse.

    Examples
    --------
    Create an irreversible mass action propensity:

    >>> prop = bcp.MassAction(k_forward=100.0)

    Create a reversible mass action propensity:

    >>> prop = bcp.MassAction(k_forward=100.0, k_reverse=0.01)
    >>> prop.is_reversible
    True

    Use with ParameterEntry objects:

    >>> kb = bcp.ParameterEntry('kb', 100.0, unit='1/(nM*s)')
    >>> ku = bcp.ParameterEntry('ku', 0.01, unit='1/s')
    >>> prop = bcp.MassAction(k_forward=kb, k_reverse=ku)

    """

    def __init__(
        self,
        k_forward: Union[float, ParameterEntry],
        k_reverse: Union[float, ParameterEntry] = None,
    ):
        super(MassAction, self).__init__()
        self.k_forward = k_forward
        self.k_reverse = k_reverse
        self.name = 'massaction'

    @property
    def k_forward(self):
        """float: Forward rate constant value."""
        if isinstance(self._k_forward, Parameter):
            return self._k_forward.value
        else:
            return self._k_forward

    @k_forward.setter
    def k_forward(self, new_k_forward):
        """Set the forward rate constant.

        Parameters
        ----------
        new_k_forward : float or ParameterEntry
            New forward rate constant. Must be positive.

        """
        self._k_forward = self._check_parameter(new_k_forward)
        self.propensity_dict['parameters']['k_forward'] = self._k_forward

    @property
    def k_reverse(self):
        """float: Reverse rate constant value, None if irreversible."""
        if isinstance(self._k_reverse, Parameter):
            return self._k_reverse.value
        else:
            return self._k_reverse

    @k_reverse.setter
    def k_reverse(self, new_k_reverse):
        """Set the reverse rate constant.

        Parameters
        ----------
        new_k_reverse : float, ParameterEntry, or None
            New reverse rate constant. If None, reaction is irreversible.
            If provided, must be positive.

        """
        self._k_reverse = self._check_parameter(
            new_k_reverse, allow_None=True
        )
        if self._k_reverse is not None:
            self.propensity_dict['parameters']['k_reverse'] = self._k_reverse

    @property
    def is_reversible(self):
        """bool: True if k_reverse is not None, False otherwise."""
        if self.k_reverse is None:
            return False
        else:
            return True

    def pretty_print_rate(self, **kwargs):
        """Generate human-readable rate formula string.

        Parameters
        ----------
        **kwargs
            Must include 'reaction' (CRN Reaction object) and 'stochastic'
            (bool) keys.

        Returns
        -------
        str
            Formatted rate formula showing forward rate and optionally
            reverse rate.

        """
        crn_reaction = kwargs['reaction']
        reactant_species = {}
        for w_species in crn_reaction.inputs:
            reactant_species[str(w_species.species)] = w_species
        txt = ' Kf=' + self._get_rate_formula(
            'k_forward', kwargs['stochastic'], reactant_species
        )
        if self.is_reversible:
            reactant_species = {}
            for w_species in crn_reaction.outputs:
                reactant_species[str(w_species.species)] = w_species
            txt += '\n Kr=' + self._get_rate_formula(
                'k_reverse', kwargs['stochastic'], reactant_species
            )
        return txt

    def create_kinetic_law(
        self,
        model,
        sbml_reaction,
        stochastic,
        crn_reaction=None,
        reverse_reaction=False,
        **kwargs,
    ):
        """Create SBML kinetic law for mass action reaction.

        Generates SBML kinetic law with proper mass action formula for
        either forward or reverse direction.

        Parameters
        ----------
        model : libsbml.Model
            SBML model object.
        sbml_reaction : libsbml.Reaction
            SBML reaction to add kinetic law to.
        stochastic : bool
            If True, uses stochastic mass action formula accounting for
            combinatorics.
        crn_reaction : Reaction
            Mass action reaction to use for the kinetic law.
        reverse_reaction : bool, default=False
            If True, creates kinetic law for reverse reaction.
        **kwargs
            Must include 'crn_reaction' (CRN Reaction object).

        Returns
        -------
        libsbml.KineticLaw
            Created SBML kinetic law object.

        Raises
        ------
        ValueError
            If crn_reaction not provided or if rate formula is invalid.

        """
        if crn_reaction is None:
            raise ValueError(
                "crn_reaction reference is needed for Massaction kinetics!"
            )

        # create a kinetic law for the sbml_reaction
        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw
        )

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
            ratelaw.removeLocalParameter('k_reverse')
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
            ratelaw.removeLocalParameter('k_forward')

        rate_formula = self._get_rate_formula(
            param, stochastic, reactant_species
        )
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. "
                "Check the propensity functions of reactions."
            )
        annotation_string = self._create_annotation(
            model, propensity_dict_in_sbml=propensity_dict_in_sbml, **kwargs
        )
        sbml_reaction.appendAnnotation(annotation_string)
        return ratelaw

    def _get_rate_formula(
        self, rate_coeff_name, stochastic, reactant_species
    ) -> str:
        r"""Generate mass action rate formula string.

        Creates the mathematical formula for mass action kinetics,
        accounting for stoichiometry and stochastic vs deterministic
        simulation.

        Parameters
        ----------
        rate_coeff_name : str
            Name of the rate constant parameter (SBML identifier).
        stochastic : bool
            If True, uses stochastic formula with combinatorics. If False,
            uses deterministic ODE formula.
        reactant_species : dict
            Dictionary mapping species_id (str) to WeightedSpecies objects
            with stoichiometry information.

        Returns
        -------
        str
            Rate formula string suitable for SBML parseL3Formula().

        Notes
        -----
        For deterministic: rate = $k * ['A']^n * ['B']^m$
        For stochastic: rate = $k * A * (A-1) * \ldots * (A-n+1) * B * \ldots$

        """
        # Create Rate-strings for massaction propensities
        ratestring = rate_coeff_name

        for species_id, weighted_species in reactant_species.items():
            if stochastic:
                ratestring += '*'
                ratestring += f"{species_id}"
                ratestring += '*'
                ratestring += '*'.join(
                    f" ( {species_id} - {i} )"
                    for i in range(1, weighted_species.stoichiometry)
                )

                # Remove trailing *
                if ratestring[len(ratestring) - 1] == '*':
                    ratestring = ratestring[:-1]
            else:
                if weighted_species.stoichiometry > 1:
                    ratestring += (
                        f" * {species_id}^"
                        + f"{weighted_species.stoichiometry}"
                    )
                else:
                    ratestring += f" * {species_id}"
        return ratestring


class Hill(Propensity):
    """Base class for Hill-type propensities.

    Hill propensities implement cooperative binding kinetics with
    sigmoidal response curves. This base class provides common
    functionality for positive and negative Hill functions.

    Parameters
    ----------
    k : float or ParameterEntry
        Maximum rate constant. Must be positive.
    s1 : Species
        Input species that drives the Hill function.
    K : float or ParameterEntry
        Half-saturation (dissociation) constant. Must be positive.
    n : float or ParameterEntry
        Hill coefficient (cooperativity). Must be positive. Values > 1
        indicate positive cooperativity, < 1 negative cooperativity.
    d : Species or None
        Optional species for proportional Hill functions. If provided,
        rate is proportional to this species concentration.

    Attributes
    ----------
    k : float
        Maximum rate constant value.
    K : float
        Half-saturation constant value.
    n : float
        Hill coefficient value.
    s1 : Species
        Input species.
    d : Species or None
        Proportional species (None for non-proportional Hill).

    See Also
    --------
    HillPositive : Positive Hill function.
    HillNegative : Negative Hill function (repression).
    ProportionalHillPositive : Proportional positive Hill.
    ProportionalHillNegative : Proportional negative Hill.

    Notes
    -----
    This is an abstract base class. Use the specific subclasses:

    - `HillPositive`: Activation, $k ['s1']^n / (K^n + ['s1']^n)$
    - `HillNegative`: Repression, $k / (1 + (['s1']/K)^n)$
    - `ProportionalHillPositive`: $k ['d'] ['s1']^n / (K^n + ['s1']^n)$
    - `ProportionalHillNegative`: $k ['d'] / (1 + (['s1']/K)^n)$

    Hill functions are not reversible - k_reverse is not supported.

    """

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
        """float: Maximum rate constant value."""
        if isinstance(self._k, Parameter):
            return self._k.value
        else:
            return self._k

    @k.setter
    def k(self, new_k):
        """Set the maximum rate constant.

        Parameters
        ----------
        new_k : float or ParameterEntry
            New maximum rate constant. Must be positive.

        """
        self._k = self._check_parameter(new_k)
        self.propensity_dict['parameters']['k'] = self._k

    @property
    def K(self):
        """float: Half-saturation (dissociation) constant value."""
        if isinstance(self._K, Parameter):
            return self._K.value
        else:
            return self._K

    @K.setter
    def K(self, new_K):
        """Set the half-saturation constant.

        Parameters
        ----------
        new_K : float or ParameterEntry
            New half-saturation constant. Must be positive.

        """
        self._K = self._check_parameter(new_K)
        self.propensity_dict['parameters']['K'] = self._K

    @property
    def n(self):
        """float: Hill coefficient (cooperativity) value."""
        if isinstance(self._n, Parameter):
            return self._n.value
        else:
            return self._n

    @n.setter
    def n(self, new_n):
        """Set the Hill coefficient.

        Parameters
        ----------
        new_n : float or ParameterEntry
            New Hill coefficient. Must be positive.

        """
        self._n = self._check_parameter(new_n)
        self.propensity_dict['parameters']['n'] = self._n

    @property
    def s1(self):
        """Species: Input species driving the Hill function."""
        return self._s1

    @s1.setter
    def s1(self, new_s1):
        """Set the input species.

        Parameters
        ----------
        new_s1 : Species
            New input species.

        """
        self._s1 = self._check_species(new_s1)
        self.propensity_dict['species']['s1'] = self.s1

    @property
    def d(self):
        """Species: Proportional species (None if not proportional)."""
        return self._d

    @d.setter
    def d(self, new_d):
        """Set the proportional species.

        Parameters
        ----------
        new_d : Species or None
            New proportional species. None for non-proportional Hill.

        """
        self._d = self._check_species(new_d, allow_None=True)
        self.propensity_dict['species']['d'] = self._d

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        """Generate human-readable rate formula string.

        Raises
        ------
        NotImplementedError
            Hill base class doesn't have a rate formula. Use subclasses
            HillPositive, HillNegative, ProportionalHillPositive, or
            ProportionalHillNegative.

        """
        raise NotImplementedError(
            "Propensity class Hill is meant to be subclassed: "
            "try HillPositive, HillNegative, ProportionalHillPositive, "
            "or ProportionalHillNegative."
        )

    def create_kinetic_law(self, model, sbml_reaction, stochastic, **kwargs):
        """Create SBML kinetic law for Hill propensity.

        This method is shared by all Hill subclasses.

        Parameters
        ----------
        model : libsbml.Model
            SBML model object.
        sbml_reaction : libsbml.Reaction
            SBML reaction to add kinetic law to.
        stochastic : bool
            If True, uses stochastic formulation (same as deterministic
            for Hill functions).
        **kwargs
            Additional arguments. 'reverse_reaction' is not supported.

        Returns
        -------
        libsbml.KineticLaw
            Created SBML kinetic law object.

        Raises
        ------
        ValueError
            If reverse_reaction=True (Hill propensities cannot be
            reversible) or if rate formula is invalid.

        """
        # This code is reused in all Hill Propensity subclasses.
        if (
            'reverse_reaction' in kwargs
            and kwargs['reverse_reaction'] is True
        ):
            raise ValueError(
                "reverse reactions cannot exist for Hill type Propensities!"
            )

        ratelaw = sbml_reaction.createKineticLaw()

        # translate the internal representation of a propensity to SBML format
        propensity_dict_in_sbml = self._translate_propensity_dict_to_sbml(
            model=model, ratelaw=ratelaw
        )

        rate_formula = self._get_rate_formula(
            propensity_dict=propensity_dict_in_sbml
        )
        # attach simulator specific annotations to the SBML model, if needed
        annotation_string = self._create_annotation(
            model, propensity_dict_in_sbml, **kwargs
        )
        sbml_reaction.appendAnnotation(annotation_string)
        # Set the ratelaw to the rateformula
        math_ast = libsbml.parseL3Formula(rate_formula)
        flag = ratelaw.setMath(math_ast)
        if not flag == libsbml.LIBSBML_OPERATION_SUCCESS or math_ast is None:
            raise ValueError(
                "Could not write the rate law for reaction to SBML. "
                "Check the propensity functions of reactions."
            )

        return ratelaw

    def _get_rate_formula(self, propensity_dict):
        """Generate Hill rate formula string.

        Parameters
        ----------
        propensity_dict : dict
            Propensity dictionary with SBML identifiers.

        Returns
        -------
        str
            Rate formula string.

        Raises
        ------
        NotImplementedError
            Hill base class doesn't have a rate formula. Use subclasses.

        """
        raise NotImplementedError(
            "Hill does not have a rate formula! Check out the subclasses."
        )


class HillPositive(Hill):
    r"""Positive Hill function propensity (activation).

    Implements an activating Hill function with sigmoidal dose-response
    curve. As s1 increases, the rate approaches $k$.

    Parameters
    ----------
    k : float or ParameterEntry
        Maximum rate constant. Must be positive.
    s1 : Species
        Input species that activates the reaction.
    K : float or ParameterEntry
        Half-saturation constant (concentration at half-maximal rate).
        Must be positive.
    n : float or ParameterEntry
        Hill coefficient (cooperativity). Must be positive. Values > 1
        indicate ultrasensitivity.

    Attributes
    ----------
    name : str
        Set to 'hillpositive' for this propensity type.

    See Also
    --------
    HillNegative : Repressive Hill function.
    ProportionalHillPositive : Proportional positive Hill.

    Notes
    -----
    The following formula is implemented:

        $$ p(s_1; k, K, n) = \frac{k s_1^n}{K^n + s_1^n}, $$

    leading to the following behaviors:

    - When $s_1 = 0$: rate ≈ 0
    - When $s_1 = K$: rate = $k$/2
    - When $s1 >> K$: rate --> $k$
    - Larger $n$ gives sharper (more switch-like) response

    Examples
    --------
    Create a Hill activation propensity:

    >>> X = bcp.Species('X')
    >>> prop = bcp.HillPositive(k=10.0, s1=X, K=50.0, n=2.0)

    Use with parameter objects:

    >>> kmax = bcp.ParameterEntry('kmax', 10.0)
    >>> Kd = bcp.ParameterEntry('Kd', 50.0)
    >>> hill_n = bcp.ParameterEntry('n', 2.0)
    >>> prop = bcp.HillPositive(k=kmax, s1=X, K=Kd, n=hill_n)

    """

    def __init__(self, k: float, s1: Species, K: float, n: float):
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=None)
        self.name = 'hillpositive'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        """Generate human-readable rate formula string.

        Returns
        -------
        str
            Formatted Hill positive formula.

        """
        return (
            f" Kf = k {self.s1.pretty_print(**kwargs)}^n / "
            + f"( K^n + {self.s1.pretty_print(**kwargs)}^n )"
        )

    def _get_rate_formula(self, propensity_dict):
        """Generate SBML-compatible Hill positive rate formula.

        Parameters
        ----------
        propensity_dict : dict
            Propensity dictionary with SBML identifiers.

        Returns
        -------
        str
            SBML rate formula string.

        """
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        rate_formula = f"{k}*{s1}^{n} / ( {K}^{n} + {s1}^{n} )"
        return rate_formula


class HillNegative(Hill):
    r"""Negative Hill function propensity (repression).

    Implements a repressive Hill function. As s1 increases, the rate
    decreases from k toward zero.

    Parameters
    ----------
    k : float or ParameterEntry
        Maximum rate constant (when s1=0). Must be positive.
    s1 : Species
        Input species that represses the reaction.
    K : float or ParameterEntry
        Half-saturation constant (concentration at half-maximal
        repression). Must be positive.
    n : float or ParameterEntry
        Hill coefficient (cooperativity). Must be positive. Values > 1
        indicate ultrasensitive repression.

    Attributes
    ----------
    name : str
        Set to 'hillnegative' for this propensity type.

    See Also
    --------
    HillPositive : Activating Hill function.
    ProportionalHillNegative : Proportional negative Hill.

    Notes
    -----
    The following mathematical formula is implemented:

        $$ p(s_1; k, K, n) = \frac{k}{1 + (s_1/K)^n} $$

    leading to the following behavior:

    - When s1 = 0: rate = $k$
    - When s1 = $K$: rate = $k$/2
    - When s1 >> $K$: rate --> 0
    - Larger $n$ gives sharper (more switch-like) repression

    Examples
    --------
    Create a Hill repression propensity:

    >>> R = bcp.Species('R')  # Repressor
    >>> prop = bcp.HillNegative(k=10.0, s1=R, K=50.0, n=2.0)

    Model transcriptional repression:

    >>> repressor = bcp.Species('repressor')
    >>> kmax = bcp.ParameterEntry('kmax', 1.0)
    >>> Ki = bcp.ParameterEntry('Ki', 100.0)
    >>> prop = bcp.HillNegative(k=kmax, s1=repressor, K=Ki, n=2.0)

    """

    def __init__(self, k: float, s1: Species, K: float, n: float):
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=None)
        self.name = 'hillnegative'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        """Generate human-readable rate formula string.

        Returns
        -------
        str
            Formatted Hill negative formula.

        """
        return f" Kf = k / ( 1 + ({self.s1.pretty_print(**kwargs)}/K)^n )"

    def _get_rate_formula(self, propensity_dict):
        """Generate SBML-compatible Hill negative rate formula.

        Parameters
        ----------
        propensity_dict : dict
            Propensity dictionary with SBML identifiers.

        Returns
        -------
        str
            SBML rate formula string.

        """
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        rate_formula = f"{k} / ( 1 + ({s1}/{K})^{n} )"
        return rate_formula


class ProportionalHillPositive(HillPositive):
    r"""Proportional positive Hill function propensity.

    Implements a positive Hill function with rate proportional to a
    species concentration. Commonly used for regulated production where
    the rate depends on both an activator and a template/enzyme.

    Parameters
    ----------
    k : float or ParameterEntry
        Maximum rate constant per unit of d. Must be positive.
    s1 : Species
        Input species that activates the reaction (e.g., transcription
        factor).
    K : float or ParameterEntry
        Half-saturation constant for s1. Must be positive.
    n : float or ParameterEntry
        Hill coefficient (cooperativity). Must be positive.
    d : Species
        Proportional species (e.g., DNA template, enzyme). Rate scales
        linearly with this species concentration.

    Attributes
    ----------
    name : str
        Set to 'proportionalhillpositive' for this propensity type.

    See Also
    --------
    HillPositive : Non-proportional positive Hill.
    ProportionalHillNegative : Proportional negative Hill.

    Notes
    -----
    The following mathematical formula: is used for the popensity:

        $$ p(s_1, d; k, K, n) = \frac{k d s_1^n}{K^n + s_1^n} $$

    This is commonly used for transcription, where

    - d = DNA template concentration
    - s1 = transcription factor concentration
    - Rate is proportional to both template and TF activation

    This results in the following behaviors:

    - When d = 0: rate = 0 (no template/enzyme)
    - When s1 = 0: rate ≈ 0 (no activation)
    - When s1 >> $K$: rate --> $k$ d (fully activated, proportional to d)

    Examples
    --------
    Model regulated transcription:

    >>> TF = bcp.Species('TF')  # Transcription factor
    >>> DNA = bcp.Species('DNA')  # DNA template
    >>> prop = bcp.ProportionalHillPositive(
    ...     k=0.1, s1=TF, K=50.0, n=2.0, d=DNA)

    Model enzyme with allosteric activator:

    >>> activator = bcp.Species('activator')
    >>> enzyme = bcp.Species('enzyme')
    >>> kcat = bcp.ParameterEntry('kcat', 10.0)
    >>> Ka = bcp.ParameterEntry('Ka', 100.0)
    >>> prop = bcp.ProportionalHillPositive(
    ...     k=kcat, s1=activator, K=Ka, n=2.0, d=enzyme)

    """

    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=d)
        self.name = 'proportionalhillpositive'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        """Generate human-readable rate formula string.

        Returns
        -------
        str
            Formatted proportional Hill positive formula.

        """
        return (
            f" Kf = k {self.d.pretty_print(**kwargs)} "
            + f"{self.s1.pretty_print(**kwargs)}^n / "
            + f"( K^n + {self.s1.pretty_print(**kwargs)}^n )"
        )

    def _get_rate_formula(self, propensity_dict):
        """Generate SBML-compatible proportional Hill positive formula.

        Parameters
        ----------
        propensity_dict : dict
            Propensity dictionary with SBML identifiers.

        Returns
        -------
        str
            SBML rate formula string.

        """
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        d = propensity_dict['species']['d']
        return f"{k}*{d}*{s1}^{n} / ( {K}^{n} + {s1}^{n} )"


class ProportionalHillNegative(HillNegative):
    r"""Proportional negative Hill function propensity.

    Implements a repressive Hill function with rate proportional to a
    species concentration. Commonly used for regulated production where
    a repressor inhibits production from a template/enzyme.

    Parameters
    ----------
    k : float or ParameterEntry
        Maximum rate constant per unit of d (when s1=0). Must be positive.
    s1 : Species
        Input species that represses the reaction (e.g., repressor).
    K : float or ParameterEntry
        Half-saturation constant for repression by s1. Must be positive.
    n : float or ParameterEntry
        Hill coefficient. Must be positive.
    d : Species
        Proportional species (e.g., DNA template, enzyme). Rate scales
        linearly with this species concentration.

    Attributes
    ----------
    name : str
        Set to 'proportionalhillnegative' for this propensity type.

    See Also
    --------
    HillNegative : Non-proportional negative Hill.
    ProportionalHillPositive : Proportional positive Hill.

    Notes
    -----
    The following mathematical formula: is used:

        $$ p(s_1, d; k, K, n) = \frac{k d}{1 + (s_1/K)^n} $$

    This is commonly used for repressed transcription where

    - d = DNA template concentration
    - s1 = repressor concentration
    - Rate is proportional to template but repressed by s1

    and resulting in the following behaviors:

    - When d = 0: rate = 0 (no template/enzyme)
    - When s1 = 0: rate = k*d (fully derepressed)
    - When s1 >> K: rate --> 0 (fully repressed)

    Examples
    --------
    Model repressed transcription:

    >>> repressor = bcp.Species('repressor')
    >>> DNA = bcp.Species('DNA')
    >>> prop = bcp.ProportionalHillNegative(
    ...     k=0.1, s1=repressor, K=50.0, n=2.0, d=DNA)

    Model enzyme with allosteric inhibitor:

    >>> inhibitor = bcp.Species('inhibitor')
    >>> enzyme = bcp.Species('enzyme')
    >>> kcat = bcp.ParameterEntry('kcat', 10.0)
    >>> Ki = bcp.ParameterEntry('Ki', 100.0)
    >>> prop = bcp.ProportionalHillNegative(
    ...     k=kcat, s1=inhibitor, K=Ki, n=2.0, d=enzyme)

    """

    def __init__(self, k: float, s1: Species, K: float, n: float, d: Species):
        Hill.__init__(self=self, k=k, s1=s1, K=K, n=n, d=d)
        self.name = 'proportionalhillnegative'

    def pretty_print_rate(self, show_parameters=True, **kwargs):
        """Generate human-readable rate formula string.

        Returns
        -------
        str
            Formatted proportional Hill negative formula.

        """
        return (
            f" Kf = k {self.d.pretty_print(**kwargs)} / "
            + f"( 1 + ({self.s1.pretty_print(**kwargs)}/K)^{self.n} )"
        )

    def _get_rate_formula(self, propensity_dict):
        """Generate SBML-compatible proportional Hill negative formula.

        Parameters
        ----------
        propensity_dict : dict
            Propensity dictionary with SBML identifiers.

        Returns
        -------
        str
            SBML rate formula string.

        """
        k = propensity_dict['parameters']['k']
        n = propensity_dict['parameters']['n']
        K = propensity_dict['parameters']['K']
        s1 = propensity_dict['species']['s1']
        d = propensity_dict['species']['d']
        return f"{k}*{d} / ( 1 + ({s1}/{K})^{n} )"
