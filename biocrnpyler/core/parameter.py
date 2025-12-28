# parameter.py - parameter processing
# RMM, 19 Aug 2018
#
# This file contains the Parameter class that is used for representing
# parameters, as well as utility functions for manipulating
# parameters.
#
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

"""Parameter processing module.

*Parameter Value Defaulting*

Not all parameters need to have the required headings.  The only two required
columns are 'param_val' and 'param_name'.  BioCRNpyler uses a form of
parameter name defaulting discussed below to find default parameters if no
exact match is in the config file. This makes it easy to set default
parameters for things like 'ku' and 'ktx' to quickly build models.

*Parameters inside BioCRNpyler*

Inside of BioCRNpyler, parameters are stored as a dictionary key value pair:
`(mechanism_name, part_id, param_name) --> param_val`. If that particular
parameter key cannot be found, the software will default to the following
keys: `(mechanism_type, part_id, param_name)` >> `(part_id, param_name)` >>
`(mechanism_name, param_name)` >> `(mechanism_type, param_name)` >>
`(param_name)` and give a warning.  As a note, `mechanism_name` refers to the
`.name` variable of a `Mechanism`, and `mechanism_type` refers to the `.type`
variable of a `Mechanism`.  Either of these can be used as a
`mechanism_id`. This allows for models to be constructed easily using default
parameter values and for parameters to be shared between different mechanisms
and/or components.

Units are read directly read from the column labeled "units" in the
parameter file.

*Initial Conditions are also Parameters*

The initial condition of any `Species` (or `Component`) will also be looked up
as a parameters automatically.  Initial conditions can be customized in
through the `custom_initial_conditions` keyword in the `Mixture` constructor.
The `custom_initial_conditions` keyword will take precedent over parameter
initial conditions.

"""

import csv
import numbers
import re
from collections import namedtuple  # Used for the parameter keys
from typing import Dict, List, Union
from warnings import warn

from ..utils.units import biocrnpyler_supported_units


class ParameterKey(namedtuple('ParameterKey', 'mechanism part_id name')):
    """Named tuple defining a parameter key.

    Parameters
    ----------
    mechanism : str, Mechanism, or None
        Mechanism identifier. Can be a string (used as both name and type), a
        Mechanism object (uses .name and .mechanism_type), or None.
    part_id : str or None
        Part/component identifier for the parameter.
    name : str
        Name of the parameter. Must start with a letter and contain at
        least one character.

    """


class Parameter(object):
    """Base class for representing parameters in BioCRNpyler.

    Parameters represent kinetic constants, initial concentrations, and
    other numerical values used in chemical reaction networks. This class
    provides validation for parameter names, values, and units.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter. Must start with a letter and contain at
        least one character.
    parameter_value : float or str
        Value of the parameter. Can be a number or a string in formats:
        '1.00', '1e4', or '2/5' (rational). Strings are automatically
        converted to numerical values.
    unit : str, optional
        Unit of the parameter (e.g., '1/s', 'nM', 'molecules'). If None,
        defaults to empty string.

    Attributes
    ----------
    parameter_name : str
        The validated name of the parameter.
    value : float
        The numerical value of the parameter.
    unit : str
        The unit string associated with the parameter.

    See Also
    --------
    ParameterEntry : Parameter with database lookup keys.
    ModelParameter : Parameter with search and found keys for defaulting.
    ParameterDatabase : Database for storing and retrieving parameters.

    Notes
    -----
    This is the base class for all parameter types in BioCRNpyler. In
    practice, subclasses `ParameterEntry` and `ModelParameter` are more
    commonly used as they support the parameter lookup and defaulting
    system.

    Parameter values provided as strings are automatically converted:

    - '1e4' --> 10000.0
    - '2/5' --> 0.4
    - '1.23' --> 1.23

    Examples
    --------
    Create a basic parameter:

    >>> param = bcp.Parameter('kb', 100.0, unit='1/s')
    >>> param.parameter_name
    'kb'
    >>> param.value
    100.0

    Create a parameter from a string value:

    >>> param = bcp.Parameter('ku', '1e-4', unit='1/s')
    >>> param.value
    0.0001

    Create a parameter from a rational string:

    >>> param = bcp.Parameter('ratio', '3/4')
    >>> param.value
    0.75

    """

    def __init__(
        self,
        parameter_name: str,
        parameter_value: Union[str, numbers.Real],
        unit=None,
    ):
        """Initialize a Parameter object.

        See class docstring for parameter descriptions.

        """
        self.parameter_name = parameter_name
        self.value = parameter_value
        self.unit = unit

    @property
    def parameter_name(self) -> str:
        """str: The name of the parameter."""
        return self._parameter_name

    @parameter_name.setter
    def parameter_name(self, new_parameter_name: str):
        """Set the parameter name with validation.

        Parameters
        ----------
        new_parameter_name : str
            New name for the parameter. Must be a string starting with a
            letter (not a number).

        Raises
        ------
        ValueError
            If `new_parameter_name` is not a string, or if it doesn't start
            with a letter.

        """
        if not isinstance(new_parameter_name, str):
            raise ValueError(
                f"parameter_name must be a string: "
                f"received {type(new_parameter_name)}."
            )
        if not re.search('^[a-z]+', new_parameter_name, re.IGNORECASE):
            raise ValueError(
                "parameter_name should be at least one character and "
                "cannot start with a number!"
            )

        self._parameter_name = new_parameter_name

    @property
    def value(self) -> numbers.Real:
        """float: The numerical value of the parameter."""
        return self._value

    @value.setter
    def value(self, new_parameter_value: Union[str, numbers.Real]):
        """Set the parameter value with validation and conversion.

        Parameters
        ----------
        new_parameter_value : float or str
            New value for the parameter. If a string, must be in format
            '1.00', '1e4', or '2/5'. Strings are automatically converted
            to float.

        Raises
        ------
        ValueError
            If `new_parameter_value` is not a number or valid string format.

        """
        if not (
            isinstance(new_parameter_value, numbers.Real)
            or isinstance(new_parameter_value, str)
        ):
            raise ValueError(
                f"parameter_value must be a float or int: "
                f"received {type(new_parameter_value)}."
            )
        if isinstance(new_parameter_value, str):
            if (
                re.search('[a-df-z]', new_parameter_value, re.I)
                or re.search(
                    '(^[1-9]+/[1-9]+)|(^[1-9]+e-?[0-9]+)|(^.?[0-9])',
                    new_parameter_value,
                    re.I,
                )
                is None
            ):
                raise ValueError(
                    f"No valid parameter value! Accepted formats: "
                    f"1.00 or 1e4 or 2/5, we got {new_parameter_value}"
                )

            self._value = Parameter._convert_rational(new_parameter_value)
        else:
            self._value = new_parameter_value

    @property
    def unit(self) -> str:
        """str: The unit string for the parameter."""
        return self._unit

    @unit.setter
    def unit(self, new_unit: str):
        """Set the parameter unit.

        Parameters
        ----------
        new_unit : str or None
            Unit string for the parameter. If None, sets to empty string.

        Raises
        ------
        ValueError
            If `new_unit` is not a string or None.

        """
        if new_unit is None:
            self._unit = ''
        elif not isinstance(new_unit, str):
            raise ValueError(
                f"All units must be strings. Recieved {new_unit}."
            )
        else:
            unit_str = new_unit.replace('/', ' per_')
            supported_units = biocrnpyler_supported_units()
            for unit in unit_str.split():
                if unit not in supported_units:
                    raise ValueError(f"unknown unit '{unit}'")
            self._unit = unit_str

    @staticmethod
    def _convert_rational(p_value: str) -> numbers.Real:
        """Convert a string parameter value to a numerical value.

        Handles rational fractions (e.g., '2/5') and standard float
        strings (e.g., '1.23' or '1e4').

        Parameters
        ----------
        p_value : str
            String representation of parameter value. Can be a fraction
            ('2/5') or standard float format ('1.23', '1e4').

        Returns
        -------
        float
            Numerical value of the parameter.

        """
        if '/' in p_value:
            nom, denom = p_value.split('/')
            return float(nom) / float(denom)
        else:
            return float(p_value)

    def __eq__(self, other):
        """Test equality between parameters or parameter and number.

        Parameters
        ----------
        other : Parameter or float
            Object to compare with. Can be another Parameter object or a
            numerical value.

        Returns
        -------
        bool
            True if values are equal, False otherwise.

        Raises
        ------
        TypeError
            If `other` cannot be compared (not a Parameter or number).

        """
        if isinstance(other, Parameter):
            return self.value == other.value
        else:
            try:
                return self.value == float(other)
            except TypeError:
                raise TypeError(
                    f"Cannot compare parameter {self} with {other}."
                )

    def __str__(self):
        """Return string representation of the parameter.

        Returns
        -------
        str
            String in format 'Parameter <name> = <value>'.

        """
        return f"Parameter {self.parameter_name} = {self.value}"

    def __hash__(self):
        """Return hash value for the parameter.

        Returns
        -------
        int
            Hash value based on parameter name, value, and unit.

        """
        return (
            hash(self._parameter_name) + hash(self._value) + hash(self._unit)
        )


class ParameterEntry(Parameter):
    """Parameter with database lookup key and metadata.

    A `ParameterEntry` extends `Parameter` with a lookup key for database
    storage and retrieval, plus additional metadata about the parameter's
    origin and context.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter.
    parameter_value : float or str
        Value of the parameter.
    parameter_key : dict, ParameterKey, str, or None, optional
        Lookup key for the parameter database.
    parameter_info : dict, optional
        Additional metadata about the parameter (e.g., source file,
        comments). If dict contains 'unit' key, it will update the
        parameter's unit.
    **kwargs
        Additional keyword arguments passed to Parameter constructor,
        including 'unit'.

    Attributes
    ----------
    parameter_key : ParameterKey
        The lookup key as a named tuple (mechanism, part_id, name).
    parameter_info : dict
        Dictionary of additional parameter metadata.

    See Also
    --------
    Parameter : Base parameter class.
    ModelParameter : Parameter with search and found keys.
    ParameterDatabase : Database for storing parameter entries.

    Notes
    -----
    The `parameter_key` value can be any of the following:

        - dict: {'mechanism': ..., 'part_id': ..., 'name': ...}
        - ParameterKey namedtuple: (mechanism, part_id, name)
        - str: parameter name (other fields set to None)
        - None: creates key with all fields None except name

    using the following conventions:

    - mechanism: str or None (mechanism name or type)
    - part_id: str or None (component/part identifier)
    - name: str (parameter name)

    These keys enable flexible parameter lookup with defaulting behavior
    in the ParameterDatabase.

    Examples
    --------
    Create a parameter entry with full key:

    >>> entry = bcp.ParameterEntry(
    ...     'kb',
    ...     100.0,
    ...     parameter_key={'mechanism': 'binding', 'part_id': 'promoter1'},
    ...     unit='1/s'
    ... )

    Create a parameter entry with just a name:

    >>> entry = bcp.ParameterEntry('ku', 0.01, parameter_key='ku')
    >>> entry.parameter_key
    ParameterKey(mechanism=None, part_id=None, name='ku')

    """

    def __init__(
        self,
        parameter_name: str,
        parameter_value: Union[str, numbers.Real],
        parameter_key=None,
        parameter_info=None,
        **kwargs,
    ):
        """Initialize a ParameterEntry object.

        See class docstring for parameter descriptions.

        """
        Parameter.__init__(self, parameter_name, parameter_value, **kwargs)

        self.parameter_key = parameter_key
        self.parameter_info = parameter_info

    # Helper function to create ParameterKeys
    @staticmethod
    def create_parameter_key(
        new_key: Union[Dict, ParameterKey, str], parameter_name=None
    ) -> ParameterKey:
        """Convert various input types to a ParameterKey namedtuple.

        Parameters
        ----------
        new_key : dict, ParameterKey, tuple, str, or None
            Input to convert to ParameterKey:

            - dict: Must have keys matching ParameterKey fields
            - ParameterKey: Returned as-is
            - 3-tuple: Converted to ParameterKey with proper field mapping
            - str: Used as 'name', other fields set to None
            - None: All fields set to None (requires parameter_name)

        parameter_name : str, optional
            Parameter name to use if not specified in new_key. Overrides
            `name` field if provided in dict.

        Returns
        -------
        ParameterKey
            Named tuple with fields (mechanism, part_id, name).

        Raises
        ------
        ValueError
            If `new_key` is not a valid type or format.

        Examples
        --------
        >>> key = bcp.ParameterEntry.create_parameter_key('kb')
        >>> key
        ParameterKey(mechanism=None, part_id=None, name='kb')

        >>> key = bcp.ParameterEntry.create_parameter_key(
        ...     {'mechanism': 'transcription', 'part_id': 'prom1'},
        ...     parameter_name='ktx'
        ... )
        >>> key
        ParameterKey(mechanism='transcription', part_id='prom1', name='ktx')

        """
        # New Key can be a named_tuple
        if isinstance(new_key, dict):
            new_key = dict(new_key)
            if parameter_name is not None:
                new_key['name'] = parameter_name
            for k in ParameterKey._fields:
                if k not in new_key:
                    new_key[k] = None
            return ParameterKey(
                **new_key
            )  # automatically unpack the keywords
        elif isinstance(new_key, ParameterKey):
            return new_key
        elif isinstance(new_key, tuple) and len(list(new_key)) == len(
            ParameterKey._fields
        ):
            # make a dictionary assuming correct ordering
            keywords = {
                ParameterKey._fields[i]: new_key[i]
                for i in range(len(ParameterKey._fields))
            }
            # automatically unpack the keywords
            return ParameterKey(**keywords)
        elif isinstance(new_key, str):
            return ParameterKey(mechanism=None, part_id=None, name=new_key)
        elif new_key is None and parameter_name is not None:
            return ParameterKey(
                mechanism=None, part_id=None, name=parameter_name
            )
        else:
            raise ValueError(
                f"parameter_key must be None, a dictionary, a ParameterKey, "
                f"a {len(ParameterKey._fields)}-tuple, or a string "
                f"(parameter name): received {new_key}."
            )

    @property
    def parameter_key(self) -> ParameterKey:
        """ParameterKey: The database lookup key for this parameter."""
        return self._parameter_key

    @parameter_key.setter
    def parameter_key(self, parameter_key: Union[Dict, ParameterKey, str]):
        """Set the parameter lookup key.

        Parameters
        ----------
        parameter_key : dict, ParameterKey, str, or None
            New parameter key. Automatically converted to ParameterKey
            namedtuple using `create_parameter_key`.

        """
        self._parameter_key = self.create_parameter_key(
            parameter_key, self.parameter_name
        )

    @property
    def parameter_info(self) -> Dict:
        """dict: Additional metadata about the parameter."""
        return self._parameter_info

    @parameter_info.setter
    def parameter_info(self, parameter_info: Dict):
        """Set parameter metadata.

        Parameters
        ----------
        parameter_info : dict or None
            Dictionary of additional parameter information. If dict
            contains 'unit' key, updates the parameter's unit attribute.

        Raises
        ------
        ValueError
            If `parameter_info` is not None or a dict, or if 'unit' in
            dict conflicts with existing unit.

        """
        if parameter_info is None:
            self._parameter_info = {}
        elif isinstance(parameter_info, dict):
            self._parameter_info = dict(parameter_info)

            # Update the units attribute, if necessary
            if (
                'unit' in parameter_info
                and self.unit != ''
                and self.unit != parameter_info['unit']
            ):
                raise ValueError(
                    f"Recieved multiple parameter units through constructor "
                    f"{self.unit} and parameter_info dictionary "
                    f"{parameter_info['unit']}."
                )
            elif 'unit' in parameter_info:
                self.unit = parameter_info['unit']

        else:
            raise ValueError(
                f"parameter_info must be None or a dictionary: "
                f"received {parameter_info}."
            )

    def get_sbml_id(self):
        """Generate SBML-compatible identifier for the parameter.

        Constructs an identifier string from the parameter key fields,
        formatted as: '<name>_<part_id>_<mechanism>'.

        Returns
        -------
        str
            SBML-compatible identifier string.

        """
        sbml_id = self.parameter_key.name + '_'
        if self.parameter_key.part_id is not None:
            sbml_id += self.parameter_key.part_id
        sbml_id += '_'
        if self.parameter_key.mechanism is not None:
            sbml_id += self.parameter_key.mechanism
        return sbml_id

    def __str__(self):
        """Return string representation of the parameter entry.

        Returns
        -------
        str
            String in format 'ParameterEntry(<key>) = <value>'.

        """
        return f"ParameterEntry({self.parameter_key}) = {self.value}"


class ModelParameter(ParameterEntry):
    """Parameter with search and found keys for defaulting behavior.

    A `ModelParameter` extends `ParameterEntry` with information about how
    the parameter was looked up in the database. It tracks both the
    original search key and the actual key where the parameter was found,
    enabling parameter defaulting and debugging.

    Parameters
    ----------
    parameter_name : str
        Name of the parameter.
    parameter_value : float or str
        Value of the parameter.
    search_key : dict, ParameterKey, tuple, or str
        The key originally searched for in the database. Usually includes
        mechanism, part_id, and name.
    found_key : dict, ParameterKey, tuple, or str
        The key where the parameter was actually found after defaulting.
        May have fewer fields than search_key.
    unit : str, optional
        Unit of the parameter.
    parameter_key : dict, ParameterKey, str, or None, optional
        Database lookup key (inherited from ParameterEntry).
    parameter_info : dict, optional
        Additional metadata (inherited from ParameterEntry).
    **kwargs
        Additional keyword arguments passed to ParameterEntry constructor.

    Attributes
    ----------
    search_key : ParameterKey
        The original lookup key as a named tuple.
    found_key : ParameterKey
        The key where parameter was found as a named tuple.

    See Also
    --------
    Parameter : Base parameter class.
    ParameterEntry : Parameter with database key.
    ParameterDatabase : Database with parameter defaulting.

    Notes
    -----
    The parameter defaulting hierarchy is:

        1. (mechanism_name, part_id, param_name)
        2. (mechanism_type, part_id, param_name)
        3. (None, part_id, param_name)
        4. (mechanism_name, None, param_name)
        5. (mechanism_type, None, param_name)
        6. (None, None, param_name)

    The `search_key` shows what was requested, while `found_key` shows
    which level of defaulting was used. This information is useful for
    debugging parameter lookups.

    Examples
    --------
    Create a model parameter showing search and found keys:

    >>> model_param = bcp.ModelParameter(
    ...     'kb',
    ...     100.0,
    ...     search_key={'mechanism': 'binding', 'part_id': 'prom1'},
    ...     found_key={'mechanism': 'binding', 'part_id': None},
    ...     unit='1/s'
    ... )
    >>> # Shows parameter was found using mechanism-level default

    """

    def __init__(
        self,
        parameter_name: str,
        parameter_value: Union[str, numbers.Real],
        search_key,
        found_key,
        unit=None,
        parameter_key=None,
        parameter_info=None,
        **kwargs,
    ):
        """Initialize a ModelParameter object.

        See class docstring for parameter descriptions.

        """
        ParameterEntry.__init__(
            self,
            parameter_name,
            parameter_value,
            unit=unit,
            parameter_key=parameter_key,
            parameter_info=parameter_info,
            **kwargs,
        )
        self.search_key = search_key
        self.found_key = found_key

    @property
    def search_key(self):
        """ParameterKey: The key originally searched for in database."""
        return self._search_key

    @search_key.setter
    def search_key(self, search_key):
        """Set the search key.

        Parameters
        ----------
        search_key : dict, ParameterKey, tuple, or str
            The key that was searched for. Automatically converted to
            ParameterKey namedtuple.

        """
        self._search_key = self.create_parameter_key(
            search_key, self.parameter_name
        )

    @property
    def found_key(self):
        """ParameterKey: The key where parameter was actually found."""
        return self._found_key

    @found_key.setter
    def found_key(self, found_key):
        """Set the found key.

        Parameters
        ----------
        found_key : dict, ParameterKey, tuple, or str
            The key where the parameter was found after defaulting.
            Automatically converted to ParameterKey namedtuple.

        """
        self._found_key = self.create_parameter_key(
            found_key, self.parameter_name
        )

    def __str__(self):
        """Return string representation of the model parameter.

        Returns
        -------
        str
            String showing parameter key, value, and search key in format
            'ModelParameter(<key>) = <value> search_key=<search_key>'.

        """
        return (
            f"ModelParameter({self.parameter_key}) = "
            + f"{self.value}\n    search_key={self.search_key}"
        )


class ParameterDatabase(object):
    """Database for storing and retrieving parameters with defaulting.

    A `ParameterDatabase` stores parameters with flexible lookup keys that
    enable parameter defaulting based on mechanism, part_id, and parameter
    name. Parameters can be loaded from dictionaries, files, or other
    databases.

    Parameters
    ----------
    parameter_dictionary : dict, optional
        Dictionary of parameters to load. Keys should be ParameterKey-like
        (dict, tuple, or str) and values should be numerical or Parameter
        objects.
    parameter_file : str or list of str, optional
        Path(s) to parameter file(s) to load. Files must be tab-separated
        (.tsv, .txt) or comma-separated (.csv).
    overwrite_parameters : bool, default=False
        If True, allows overwriting existing parameters when loading. If
        False, raises ValueError if duplicate keys are encountered.

    Attributes
    ----------
    parameters : dict
        Internal dictionary mapping ParameterKey to ParameterEntry objects.
        Access via indexing or iteration rather than directly.

    See Also
    --------
    Parameter : Base parameter class.
    ParameterEntry : Parameter with database key.
    ModelParameter : Parameter with search and found keys.

    Notes
    -----
    When searching for a parameter with `find_parameter(mechanism,
    part_id, param_name)`, the database searches in this order:

    1. (mechanism.name, part_id, param_name)
    2. (mechanism.type, part_id, param_name)
    3. (None, part_id, param_name)
    4. (mechanism.name, None, param_name)
    5. (mechanism.type, None, param_name)
    6. (None, None, param_name)

    This enables flexible parameter specification where specific parameters
    override more general ones.

    Parameter files should have these columns (column names are flexible):

    - 'param_name' or 'parameter_name' (required)
    - 'param_val' or 'value' (required)
    - 'mechanism_id' or 'mechanism' (optional)
    - 'part_id' or 'part' (optional)
    - 'units' or 'unit' (optional)

    Additional columns are stored in parameter_info.

    Parameter files are searched for in the following directories:

    1. The current directory
    2. All directories listed in the 'BCP_PATH' environment variable
    3. The BioCRNpyler source code directory

    The directories are search in this order and the first parameter file that
    is found is returned.  For files in the BioCRNpyler source code directory,
    common filename patterns are of the form '<type>/<name>_parameters.tsv'
    where <type> is 'components', 'mechanisms', or 'mixtures'.

    Examples
    --------
    Create a parameter database from a dictionary:

    >>> params = {
    ...     'kb': 100.0,
    ...     'ku': 0.01,
    ...     ('transcription', None, 'ktx'): 0.05
    ... }
    >>> db = bcp.ParameterDatabase(parameter_dictionary=params)

    Load parameters from a file:

    >>> db = bcp.ParameterDatabase(
    ...     parameter_file='mixtures/pure_parameters.tsv')

    Look up a parameter with defaulting:

    >>> param = db.find_parameter('transcription', 'promoter1', 'ktx')
    >>> param.value
    0.05

    Add a new parameter:

    >>> db.add_parameter('kcat', 10.0,
    ...     parameter_key={'mechanism': 'catalysis', 'part_id': 'enzyme1'})

    """

    def __init__(
        self,
        parameter_dictionary=None,
        parameter_file=None,
        overwrite_parameters=False,
    ):
        """Initialize a ParameterDatabase object.

        See class docstring for parameter descriptions.

        """
        self.parameters = {}  # create an emtpy dictionary to get parameters.

        if isinstance(parameter_file, str):
            self.load_parameters_from_file(
                parameter_file, overwrite_parameters=overwrite_parameters
            )
        elif isinstance(parameter_file, list):
            for p in parameter_file:
                if isinstance(p, str):
                    self.load_parameters_from_file(
                        p, overwrite_parameters=overwrite_parameters
                    )
                else:
                    raise ValueError(
                        "parameter_file must be a string or list of strings "
                        "representing file names and paths."
                    )
        elif parameter_file is not None:
            raise ValueError(
                "parameter_file must be a string representing a file name "
                "and path."
            )

        if isinstance(parameter_dictionary, dict):
            self.load_parameters_from_dictionary(
                parameter_dictionary,
                overwrite_parameters=overwrite_parameters,
            )
        elif parameter_dictionary is not None:
            raise ValueError(
                "parameter_dictionary must be None or a dictionary!"
            )

    # To check if a key or ParameterEntry is in a the ParameterDatabase
    def __contains__(self, val):
        """Check if a key or ParameterEntry is in the database.

        Parameters
        ----------
        val : ParameterEntry, dict, ParameterKey, tuple, or str
            Value to check. Can be a ParameterEntry object or any valid
            parameter key format.

        Returns
        -------
        bool
            True if the key exists in the database (and ParameterEntry
            values match if val is a ParameterEntry), False otherwise.

        """
        if isinstance(val, ParameterEntry):
            key = val.parameter_key
            if key in self.parameters and self.parameters[key] == val:
                return True
            else:
                return False
        else:
            try:
                key = ParameterEntry.create_parameter_key(val)
                return key in self.parameters
            except ValueError:
                return False

    # Ability to loop through parameters eg
    # for entry in ParameterDatabase: ...
    def __iter__(self):
        """Initialize iterator over parameter entries.

        Returns
        -------
        ParameterDatabase
            Self with iterator state initialized.

        """
        self.keys = list(self.parameters.keys())
        self.current_key_ind = 0
        return self

    def __next__(self):
        """Get next parameter entry in iteration.

        Returns
        -------
        ParameterEntry
            Next parameter entry in the database.

        Raises
        ------
        StopIteration
            When all parameters have been iterated.

        """
        if self.current_key_ind < len(self.keys):
            key = self.keys[self.current_key_ind]
            entry = self.parameters[key]
            self.current_key_ind += 1
            return entry
        else:
            raise StopIteration

    # Length method
    def __len__(self):
        """Return number of parameters in the database.

        Returns
        -------
        int
            Number of parameter entries stored.

        """
        return len(self.parameters)

    # Gets a parameter from the database
    # Only returns exact matches.
    def __getitem__(self, key):
        """Get a parameter by exact key match.

        Parameters
        ----------
        key : dict, ParameterKey, tuple, or str
            Parameter key to look up. No defaulting is performed.

        Returns
        -------
        ParameterEntry
            The parameter entry with the exact matching key.

        Raises
        ------
        KeyError
            If the exact key is not found in the database.

        """
        param_key = ParameterEntry.create_parameter_key(key)
        return self.parameters[param_key]

    # Sets a parameter in the databases - useful for quickly changing
    # parameters, but add_parameter is recommended.
    def __setitem__(self, parameter_key, value):
        """Set a parameter value by key.

        Parameters
        ----------
        parameter_key : dict, ParameterKey, tuple, or str
            Key for the parameter.
        value : float, str, or ParameterEntry
            New value or ParameterEntry object. If ParameterEntry, its key
            must match parameter_key.

        Raises
        ------
        ValueError
            If value is ParameterEntry with mismatched key.

        Notes
        -----
        This method automatically overwrites existing parameters.
        For more control, use `add_parameter` instead.

        """
        key = ParameterEntry.create_parameter_key(parameter_key)

        if isinstance(value, ParameterEntry):
            if key != value.parameter_key:
                raise ValueError(
                    f"Parameter Key does not match: ParameterDatabase key "
                    f"{key} is not the same as ParameterEntry Key "
                    f"{value.parameter_key}."
                )
            self.parameters[key] = value
        else:
            self.add_parameter(
                key.name,
                value,
                parameter_key=key,
                parameter_origin='Set Manually',
                overwrite_parameters=True,
            )

    def __str__(self):
        """Return string representation of the parameter database.

        Returns
        -------
        str
            String listing all parameters in the database.

        """
        txt = 'ParameterDatabase:'
        param_txt = '\n'.join([repr(p) for p in self.parameters])
        return txt + param_txt

    def add_parameter(
        self,
        parameter_name: str,
        parameter_value: Union[str, numbers.Real],
        parameter_origin=None,
        parameter_key=None,
        parameter_info=None,
        overwrite_parameters=False,
    ):
        """Add a parameter to the database.

        Parameters
        ----------
        parameter_name : str
            Name of the parameter.
        parameter_value : float or str
            Value of the parameter. Strings are converted to float.
        parameter_origin : str, optional
            Description of where the parameter came from (e.g., filename).
            Stored in `parameter_info`.
        parameter_key : dict, ParameterKey, str, or None, optional
            Lookup key for the parameter. If None, creates key with only
            the parameter name.
        parameter_info : dict, optional
            Additional metadata about the parameter. `parameter_origin` is
            added to this dict if provided.
        overwrite_parameters : bool, default=False
            If True, allows overwriting existing parameters. If False,
            raises ValueError if key already exists.

        Raises
        ------
        ValueError
            If key already exists in database and
            `overwrite_parameters=False`.

        Examples
        --------
        >>> db = bcp.ParameterDatabase()
        >>> db.add_parameter('kb', 100.0)
        >>> db.add_parameter('ku', 0.01,
        ...     parameter_key={'mechanism': 'binding'})

        """
        # Put parameter origin into parameter_info
        if parameter_info is None:
            parameter_info = {}
        if 'parameter origin' not in parameter_info:
            parameter_info['parameter origin'] = parameter_origin

        # Create ParameterEntry
        param = ParameterEntry(
            parameter_name,
            parameter_value,
            parameter_key=parameter_key,
            parameter_info=parameter_info,
        )
        key = param.parameter_key

        # Update parameter dictionary
        if key in self.parameters and not overwrite_parameters:
            raise ValueError(
                f"Duplicate parameter detected. Parameter with key = "
                f"{key} is already in the ParameterDatabase. To Overwrite "
                f"existing parameters, use overwrite_parameters = True."
            )
        else:
            self.parameters[key] = param

    def load_parameters_from_dictionary(
        self,
        parameter_dictionary: Dict[ParameterKey, Union[str, numbers.Real]],
        overwrite_parameters=False,
    ) -> None:
        """Load parameters from a dictionary.

        Parameters
        ----------
        parameter_dictionary : dict
            Dictionary with keys as ParameterKey-like objects (dict, tuple,
            or str) and values as numerical values or strings.
        overwrite_parameters : bool, default=False
            If True, allows overwriting existing parameters. If False,
            raises ValueError if duplicate keys are encountered.

        Raises
        ------
        ValueError
            If duplicate keys exist and `overwrite_parameters=False`.

        Examples
        --------
        >>> db = bcp.ParameterDatabase()
        >>> params = {
        ...     'kb': 100.0,
        ...     ('binding', None, 'ku'): 0.01,
        ... }
        >>> db.load_parameters_from_dictionary(params)

        """
        for k in parameter_dictionary:
            key = ParameterEntry.create_parameter_key(k)
            self.add_parameter(
                key.name,
                parameter_dictionary[k],
                parameter_key={
                    'part_id': key.part_id,
                    'mechanism': key.mechanism,
                },
                parameter_origin='parameter_dictionary',
                overwrite_parameters=overwrite_parameters,
            )

    def load_parameters_from_database(
        self, parameter_database, overwrite_parameters=False
    ) -> None:
        """Load parameters from another `ParameterDatabase`.

        Parameters
        ----------
        parameter_database : ParameterDatabase
            Another ParameterDatabase instance to copy parameters from.
        overwrite_parameters : bool, default=False
            If True, allows overwriting existing parameters. If False,
            raises ValueError if duplicate keys are encountered.

        Raises
        ------
        TypeError
            If `parameter_database` is not a ParameterDatabase instance.
        ValueError
            If duplicate keys exist and `overwrite_parameters=False`.

        Examples
        --------
        >>> db1 = bcp.ParameterDatabase(parameter_dictionary={'kb': 100.0})
        >>> db2 = bcp.ParameterDatabase()
        >>> db2.load_parameters_from_database(db1)
        >>> db2['kb'].value
        100.0

        """
        if not isinstance(parameter_database, ParameterDatabase):
            raise TypeError(
                f"paramater_database must be a ParamaterDatabase: "
                f"recievied {parameter_database}."
            )

        for k in parameter_database:
            if k not in self.parameters or overwrite_parameters:
                self.parameters[k.parameter_key] = parameter_database[
                    k.parameter_key
                ]
            else:
                raise ValueError(
                    f"Duplicate parameter detected. Parameter with key = "
                    f"{k} is already in the ParameterDatabase. To Overwrite "
                    f"existing parameters, use overwrite_parameters = True."
                )

    def load_parameters_from_file(
        self, filename: str, overwrite_parameters=False
    ) -> None:
        """Load parameters from a file.

        Reads parameters from a CSV or TSV file and adds them to the
        database. The file must have 'param_name' and 'param_val' columns.
        Optional columns include 'mechanism', 'part_id', and 'units'.

        Parameters
        ----------
        filename : str
            Path to parameter file. Must be tab-separated (.tsv, .txt) or
            comma-separated (.csv). File is searched in current directory
            and BioCRNpyler package paths.
        overwrite_parameters : bool, default=False
            If True, allows overwriting existing parameters. If False,
            raises ValueError if duplicate keys are encountered.

        Raises
        ------
        ValueError
            If file cannot be found, has invalid format, or contains
            duplicate keys when `overwrite_parameters=False`.

        Notes
        -----
        Accepted column names (case-sensitive, first match used):

        - param_name: 'parameter_name', 'parameter', 'param', 'param_name'
        - param_val: 'val', 'value', 'param_val', 'parameter_value'
        - mechanism: 'mechanism', 'mechanism_id'
        - part_id: 'part_id', 'part'
        - unit: 'units', 'unit'

        Additional columns are stored in parameter_info dictionary.

        File Format Example (CSV)::

        .. code::

            mechanism,part_id,param_name,param_val,unit
            binding,,kb,100,1/s
            binding,,ku,0.01,1/s
            transcription,prom1,ktx,0.05,1/s

        Examples
        --------
        >>> db = bcp.ParameterDatabase(
        ...    parameter_file='mixtures/pure_parameters.tsv')

        Load multiple files:

        >>> db = bcp.ParameterDatabase(
        ...     parameter_file=[
        ...         'mixtures/pure_parameters.tsv',
        ...         'components/tetr_parameters.tsv'])

        """
        from ..utils.fileutil import find_file_in_bcp_path

        # Find the file in the current path
        if (_filename := find_file_in_bcp_path(filename)) is None:
            raise ValueError(f"can't find parameter file {filename}")
        else:
            filename = _filename

        # Figure out the format of the parameter file from the file extension
        with open(filename) as f:
            file_type = filename.split('.')[-1]
            if file_type in ['tsv', 'txt']:
                delimiter = '\t'
            elif file_type in ['csv']:
                delimiter = ','
            else:
                raise ValueError(
                    "Parameter files must be tab-seperated (.tsv or .txt) "
                    "or comma-seperated (.csv) files."
                )

            # Read the CSV file, filtering out comment lines
            csvreader = csv.DictReader(
                filter(lambda row: row[0] != '#', f), delimiter=delimiter
            )

            # Used for flexible column headings
            accepted_field_names = {
                'mechanism': ['mechanism', 'mechanism_id'],
                'param_name': [
                    'parameter_name',
                    'parameter',
                    'param',
                    'param_name',
                ],
                'part_id': ['part_id', 'part'],
                'param_val': ['val', 'value', 'param_val', 'parameter_value'],
                'unit': ['unit', 'units'],
            }

            field_names = self._get_field_names(
                csvreader.fieldnames, accepted_field_names
            )

            # Determine which columns are in the CSV
            if field_names['param_name'] is None:
                warn(
                    "No param_name column was found, "
                    "could not load parameter!"
                )
            if field_names['mechanism'] is None:
                no_mechism_column = True
            else:
                no_mechism_column = False

            if field_names['part_id'] is None:
                no_part_id_column = True
            else:
                no_part_id_column = False

            # Load all parameters
            for row in csvreader:
                param_value = row[field_names['param_val']]
                field_columns = [
                    field_names['param_name'],
                    field_names['part_id'],
                    field_names['mechanism'],
                    field_names['param_val'],
                ]
                parameter_info = {
                    k: row[k] for k in row if k not in field_columns
                }
                # TODO test all these cases!

                # Case 1: No Param Name so skip the row
                if (
                    row[field_names['param_name']] is None
                    or len(row[field_names['param_name']]) == 0
                ):
                    pass

                # Case 2: Just a Param Name
                elif no_mechism_column and no_part_id_column:
                    param_name = row[field_names['param_name']]
                    self.add_parameter(
                        param_name,
                        param_value,
                        parameter_origin=filename,
                        parameter_info=parameter_info,
                        overwrite_parameters=overwrite_parameters,
                    )

                # Case 3: Part_id and Param Name
                elif no_mechism_column and no_part_id_column is False:
                    param_name = row[field_names['param_name']]
                    part_id = row[field_names['part_id']]

                    if part_id is not None and len(part_id) > 0:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_key={'part_id': part_id},
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )
                    else:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )

                # Case 4: mechanism and param name
                elif no_part_id_column and no_mechism_column is False:
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if mech_name is not None and len(mech_name) > 0:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_key={'mechanism': mech_name},
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )
                    else:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )

                # Case 5: mechanism, part_id, and param name
                else:
                    part_id = row[field_names['part_id']]
                    mech_name = row[field_names['mechanism']]
                    param_name = row[field_names['param_name']]
                    if (
                        part_id is not None
                        and len(part_id) > 0
                        and mech_name is not None
                        and len(mech_name) > 0
                    ):
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_key={
                                'part_id': part_id,
                                'mechanism': mech_name,
                            },
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )

                    elif part_id is not None and len(part_id) > 0:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_key={'part_id': part_id},
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )

                    elif mech_name is not None and len(mech_name) > 0:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_key={'mechanism': mech_name},
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )
                    else:
                        self.add_parameter(
                            param_name,
                            param_value,
                            parameter_origin=filename,
                            parameter_info=parameter_info,
                            overwrite_parameters=overwrite_parameters,
                        )

    @staticmethod
    def _get_field_names(
        field_names: List[str], accepted_field_names: Dict[str, List[str]]
    ) -> Dict[str, str]:
        """Map parameter file column names to standard field names.

        Searches through column names in a parameter file to find which
        valid aliases are being used, and creates a mapping dictionary.

        Parameters
        ----------
        field_names : list of str
            List of column names found in the parameter file.
        accepted_field_names : dict
            Dictionary mapping standard field names to lists of valid
            aliases. Format: {'field': ['alias1', 'alias2', ...]}.

        Returns
        -------
        dict
            Dictionary mapping standard field names to the actual column
            names used in the file. Fields not found are set to None.

        Raises
        ------
        ValueError
            If `field_names` or `accepted_field_names` are invalid types or
            empty.

        Warns
        -----
        UserWarning
            If a standard field has no matching column in the file.

        """
        if not isinstance(field_names, list):
            raise ValueError("field_names must be a list of strings")
        if isinstance(field_names, list) and len(field_names) == 0:
            raise ValueError("field_names cannot be empty list!")
        if not isinstance(accepted_field_names, dict):
            raise ValueError("accepted_field_names must be a dictionary")
        if (
            isinstance(accepted_field_names, dict)
            and len(accepted_field_names) == 0
        ):
            raise ValueError(
                "accepted_field_names cannot be empty dictionary"
            )

        return_field_names = dict.fromkeys(accepted_field_names.keys())
        for accepted_name in accepted_field_names:
            # try to find an possible accepted names in the
            # field_names using a generator
            try:
                loc_gen = (
                    idx
                    for idx, name in enumerate(
                        accepted_field_names[accepted_name]
                    )
                    if name in field_names
                )
                loc_idx = next(loc_gen)
            except StopIteration:
                # we have reached the end of the possible names
                return_field_names[accepted_name] = None
                warn(
                    f"parameter file contains no {accepted_name} column! "
                    "Please add a column named "
                    f"{accepted_field_names[accepted_name]}."
                )
            else:
                return_field_names[accepted_name] = accepted_field_names[
                    accepted_name
                ][loc_idx]

        return return_field_names

    def find_parameter(self, mechanism, part_id, param_name):
        """Search for a parameter with automatic defaulting.

        Searches the database for the best matching parameter using a
        hierarchical defaulting system. If an exact match is not found,
        progressively more general keys are tried.

        Parameters
        ----------
        mechanism : str, Mechanism, or None
            Mechanism identifier. Can be a string (used as both name and
            type), a Mechanism object (uses .name and .mechanism_type), or
            None.
        part_id : str or None
            Part/component identifier for the parameter.
        param_name : str
            Name of the parameter to find.

        Returns
        -------
        ModelParameter or None
            ModelParameter object with search_key and found_key attributes
            showing how the parameter was found. Returns None if no match
            found at any defaulting level.

        Raises
        ------
        ValueError
            If `mechanism` is not a string, Mechanism object, or None.

        Notes
        -----
        The method searches for parameters in this order:

        1. (mechanism.name, part_id, param_name)
        2. (mechanism.type, part_id, param_name)
        3. (None, part_id, param_name)
        4. (mechanism.name, None, param_name)
        5. (mechanism.type, None, param_name)
        6. (None, None, param_name)

        This allows setting default parameters at various levels of
        specificity. For example, a general 'kb' parameter can be
        overridden for specific mechanisms or parts.

        Examples
        --------
        >>> db = bcp.ParameterDatabase()
        >>> db.add_parameter('kb', 100.0)
        >>> db.add_parameter('kb', 200.0,
        ...     parameter_key={'mechanism': 'binding'})

        General lookup finds the general parameter

        >>> param = db.find_parameter(None, None, 'kb')
        >>> param.value
        100.0

        Mechanism-specific lookup finds the specific parameter

        >>> param = db.find_parameter('binding', None, 'kb')
        >>> param.value
        200.0

        """
        # Parameter defaulting hierarchy:
        # (mechanism_name, part_id, param_name) --> param_val.
        # If that particular parameter key cannot be found,
        # the software will default to the following keys:
        # (mechanism_type, part_id, param_name) >> (part_id, param_name) >>
        # (mechanism_name, param_name) >> (mechanism_type, param_name) >>
        # (param_name) and give a warning.
        #
        # As a note, mechanism_name refers to the .name variable of a
        # Mechanism.  mechanism_type refers to the .type variable of a
        # Mechanism.  Either of these can be used as a mechanism_id.  This
        # allows for models to be constructed easily using default parameter
        # values and for parameters to be shared between different Mechanisms
        # and/or Components.
        # this is imported here because otherwise there are import loops
        from .mechanism import Mechanism

        found_entry = None

        if isinstance(mechanism, str):
            mech_name = mechanism
            mech_type = mechanism
        elif isinstance(mechanism, Mechanism):
            mech_name = mechanism.name
            mech_type = mechanism.mechanism_type
        elif mechanism is not None:
            raise ValueError(
                f"mechanism keyword must be or string or have name and "
                f"mechanism_type attributes: recievied {mechanism}."
            )
        else:
            mech_name = None
            mech_type = None

        if not isinstance(part_id, list):
            part_id = [part_id]
        part_id = part_id + [None]
        parameter_key_list = []

        # Create a parameter key for each part_id
        for id in part_id:
            parameter_key_list += [
                ParameterKey(
                    mechanism=mech_name, part_id=id, name=param_name
                ),
                ParameterKey(
                    mechanism=mech_type, part_id=id, name=param_name
                ),
                ParameterKey(mechanism=None, part_id=id, name=param_name),
            ]

        for key in parameter_key_list:
            if key in self.parameters and found_entry is None:
                found_entry = self.parameters[key]
                found_key = key
                break

        if found_entry is None:
            return None
        else:
            return_param = ModelParameter(
                found_entry.parameter_name,
                found_entry.value,
                (mech_name, part_id, param_name),
                found_key,
                parameter_key=found_entry.parameter_key,
                parameter_info=found_entry.parameter_info,
                unit=found_entry.unit,
            )
            return return_param
