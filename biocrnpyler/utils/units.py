# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import warnings

import libsbml  # type: ignore


_SIMPLE_UNIT_ALIASES = {
    's': 'sec',
    'second': 'sec',
    'secs': 'sec',
    'seconds': 'sec',
    'minute': 'min',
    'mins': 'min',
    'minutes': 'min',
    'hour': 'hrs',
    'hr': 'hrs',
    'hours': 'hrs',
}

def normalize_unit_id(unit_id: str) -> str:
    if unit_id is None:
        return ''
    unit_id = unit_id.strip()
    if unit_id == '':
        return ''

    unit_id = unit_id.replace(' per_', '/').replace(' per ', '/')
    unit_id = unit_id.replace('(', '').replace(')', '').replace('*', '/')
    unit_id = unit_id.replace(' ', '')

    def canon(token):
        if token == '1':
            return '1'
        return _SIMPLE_UNIT_ALIASES.get(token, token)

    if '/' not in unit_id:
        return canon(unit_id)

    pieces = [canon(p) for p in unit_id.split('/') if p]

    if unit_id.startswith('/') or (pieces and pieces[0] == '1'):
        denominators = pieces[1:] if pieces and pieces[0] == '1' else pieces
        return 'per_' + '_per_'.join(denominators)

    numerator = pieces[0]
    denominators = pieces[1:]
    return numerator + ''.join(f'_per_{d}' for d in denominators)


def biocrnpyler_supported_units():
    supported_units = {
        # Volume units
        'nL': {
            'unit_kind': [libsbml.UNIT_KIND_LITRE],
            'unit_exponents': 1,
            'unit_scale': -9,
            'unit_multiplier': 1,
        },
        'uL': {
            'unit_kind': libsbml.UNIT_KIND_LITRE,
            'unit_exponents': [1],
            'unit_scale': -6,
            'unit_multiplier': 1,
        },
        'mL': {
            'unit_kind': libsbml.UNIT_KIND_LITRE,
            'unit_exponents': 1,
            'unit_scale': [-3],
            'unit_multiplier': 1,
        },
        'L': {
            'unit_kind': libsbml.UNIT_KIND_LITRE,
            'unit_exponents': 1,
            'unit_scale': 0,
            'unit_multiplier': [1],
        },
        # Concentration units
        'M': {
            'unit_kind': [libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
            'unit_exponents': [1, -1],
            'unit_scale': [0, 0],
            'unit_multiplier': [1, 1],
        },
        'mM': {
            'unit_kind': [libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
            'unit_exponents': [1, -1],
            'unit_scale': [-3, 0],
            'unit_multiplier': [1, 1],
        },
        'uM': {
            'unit_kind': [libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
            'unit_exponents': [1, -1],
            'unit_scale': [-6, 0],
            'unit_multiplier': [1, 1],
        },
        'nM': {
            'unit_kind': [libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
            'unit_exponents': [1, -1],
            'unit_scale': [-9, 0],
            'unit_multiplier': [1, 1],
        },
        # Time units
        'hour': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [3600],
        },
        'minute': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [60],
        },
        'second': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [1],
        },
        # Common parameter units
        '1': {
            'unit_kind': [libsbml.UNIT_KIND_DIMENSIONLESS],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [1],
        },
        'per_second': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [-1],
            'unit_scale': [0],
            'unit_multiplier': [1],
        },
        'per_minute': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [-1],
            'unit_scale': [0],
            'unit_multiplier': [1 / 60],
        },
        'per_hour': {
            'unit_kind': [libsbml.UNIT_KIND_SECOND],
            'unit_exponents': [-1],
            'unit_scale': [0],
            'unit_multiplier': [1 / 3600],
        },
        'mole_per_litre': {
            'unit_kind': [libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE],
            'unit_exponents': [-1, 1],
            'unit_scale': [0, 0],
            'unit_multiplier': [1, 1],
        },
        'litre_per_mole_per_second': {
            'unit_kind': [
                libsbml.UNIT_KIND_SECOND,
                libsbml.UNIT_KIND_LITRE,
                libsbml.UNIT_KIND_MOLE,
            ],
            'unit_exponents': [-1, 1, -1],
            'unit_scale': [0, 0, 0],
            'unit_multiplier': [1, 1, 1],
        },
        'litre_per_mole_per_hour': {
            'unit_kind': [
                libsbml.UNIT_KIND_SECOND,
                libsbml.UNIT_KIND_LITRE,
                libsbml.UNIT_KIND_MOLE,
            ],
            'unit_exponents': [-1, 1, -1],
            'unit_scale': [0, 0, 0],
            'unit_multiplier': [1 / 3600, 1, 1],
        },
        'amol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -18,
            'unit_multiplier': 1,
        },
        'fmol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -15,
            'unit_multiplier': 1,
        },
        'pmol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -12,
            'unit_multiplier': 1,
        },
        'nmol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -9,
            'unit_multiplier': 1,
        },
        'umol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -6,
            'unit_multiplier': 1,
        },
        'mmol': {
            'unit_kind': libsbml.UNIT_KIND_MOLE,
            'unit_exponents': 1,
            'unit_scale': -3,
            'unit_multiplier': 1,
        },
        'aa': {
            'unit_kind': [libsbml.UNIT_KIND_DIMENSIONLESS],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [1],
        },
        'nt': {
            'unit_kind': [libsbml.UNIT_KIND_DIMENSIONLESS],
            'unit_exponents': [1],
            'unit_scale': [0],
            'unit_multiplier': [1],
        },

    }
    ### Add your own units to this dictionary ###

    # Aliases
    for hour_alias in ['hrs', 'hr', 'hours']:
        supported_units[hour_alias] = supported_units['hour']
    for minute_alias in ['min', 'mins', 'minutes']:
        supported_units[minute_alias] = supported_units['minute']
    for second_alias in ['s', 'sec', 'secs', 'seconds']:
        supported_units[second_alias] = supported_units['second']
    # since aliases are allowed for time, 
    # need these for the per_ version too
    supported_units['per_sec'] = supported_units['per_second']
    supported_units['per_min'] = supported_units['per_minute']
    supported_units['per_hrs'] = supported_units['per_hour']
    
    # adding concentration per time units:
    conc_scales = {'nM': -9, 'uM': -6, 'mM': -3, 'M': 0}
    for conc_name, conc_scale in conc_scales.items():
        supported_units[f'{conc_name}_per_sec'] = {
            'unit_kind': [
                libsbml.UNIT_KIND_MOLE,
                libsbml.UNIT_KIND_LITRE,
                libsbml.UNIT_KIND_SECOND,
            ],
            'unit_exponents': [1, -1, -1],
            'unit_scale': [conc_scale, 0, 0],
            'unit_multiplier': [1, 1, 1],
        }
        supported_units[f'per_{conc_name}_per_sec'] = {
            'unit_kind': [
                libsbml.UNIT_KIND_SECOND,
                libsbml.UNIT_KIND_LITRE,
                libsbml.UNIT_KIND_MOLE,
            ],
            'unit_exponents': [-1, 1, -1],
            'unit_scale': [0, 0, conc_scale],
            'unit_multiplier': [1, 1, 1],
        }
        supported_units[f'per_sec_per_{conc_name}'] = {
            'unit_kind': [
                libsbml.UNIT_KIND_SECOND,
                libsbml.UNIT_KIND_LITRE,
                libsbml.UNIT_KIND_MOLE,
            ],
            'unit_exponents': [-1, 1, -1],
            'unit_scale': [0, 0, conc_scale],
            'unit_multiplier': [1, 1, 1],
        }



    return supported_units


def create_new_unit_definition(model, unit_id):
    """
    Creates UnitDefinition inside SBML Model object.

    Parameters
    ----------
    model : libsbml.Model
    unit_id : str

    Returns
    -------
    libsbml.UnitDefinition
        A pointer to the new libSBML object created for the unit type.

    """
    supported_units = biocrnpyler_supported_units()
    if not isinstance(unit_id, str):
        raise ValueError(
            "Arguments are not of expected type. unit_id must be a string."
        )
    if unit_id not in supported_units.keys():
        warnings.warn(
            f"The string identifier for the unit {unit_id} is not supported "
            "by BioCRNpyler. Add this to the dictionary in "
            "biocrnpyler/units.py if you want this unit."
        )
        return None
    unit_kind = supported_units[unit_id]['unit_kind']
    unit_exponents = supported_units[unit_id]['unit_exponents']
    unit_scale = supported_units[unit_id]['unit_scale']
    unit_multiplier = supported_units[unit_id]['unit_multiplier']

    unitdef = model.createUnitDefinition()
    unitdef.setId(unit_id)
    # Scale list
    if isinstance(unit_kind, list) and not isinstance(unit_scale, list):
        unit_scaleList = []
        for kind in unit_kind:
            unit_scaleList.append(unit_scale)

    elif isinstance(unit_kind, list) and isinstance(unit_scale, list):
        if len(unit_kind) != len(unit_scale):
            raise ValueError(
                "Lengths of unit_scale and unit kind lists are not equal"
            )
        unit_scaleList = unit_scale[:]

    # Multiplier list
    if isinstance(unit_kind, list) and not isinstance(unit_multiplier, list):
        unit_multiplierList = []
        for kind in unit_kind:
            unit_multiplierList.append(unit_multiplier)

    elif isinstance(unit_kind, list) and isinstance(unit_multiplier, list):
        if len(unit_kind) != len(unit_multiplier):
            raise ValueError(
                "Lengths of unit kind and unit_multiplier lists are not equal"
            )
        unit_multiplierList = unit_multiplier[:]

    if not isinstance(unit_kind, list):
        unit_kind = [unit_kind]
        if isinstance(unit_scale, list):
            unit_scale = unit_scale[0]
        if isinstance(unit_multiplier, list):
            unit_multiplier = unit_multiplier[0]
        if isinstance(unit_exponents, list):
            unit_exponents = unit_exponents[0]
        if not isinstance(unit_scale, int) or not isinstance(
            unit_multiplier, int
        ):
            raise ValueError(
                "Scale and unit_multiplier must be integers when there is "
                "only one unit kind"
            )
        unit_scaleList = [unit_scale]
        unit_multiplierList = [unit_multiplier]
    if not isinstance(unit_exponents, list):
        if not isinstance(unit_exponents, int):
            raise ValueError("All unit_exponentss should be integers")
        unit_exponents = [unit_exponents]
    if len(unit_kind) != len(unit_exponents):
        raise ValueError(
            "Lengths of unit kind and unit unit_exponents lists must be equal"
        )

    for kind, expo, unit_scale, unit_multiplier in zip(
        unit_kind, unit_exponents, unit_scaleList, unit_multiplierList
    ):
        unit = unitdef.createUnit()
        unit.setKind(kind)
        unit.setExponent(expo)
        unit.setScale(unit_scale)
        unit.setMultiplier(unit_multiplier)
    return unitdef


#
# Unit conversions
#
# To allow easy conversion of units, we define a number of variables
# that convert to the default units of uM, uL, sec.

# Concentration
nM = 1e-3
uM = 1
mM = 1e3
M = 1e6

# Volumes
nL = 1e-3
uL = 1
mL = 1e3
L = 1e6

# Time
sec = secs = 1
min = mins = 60
hr = hrs = 3600
