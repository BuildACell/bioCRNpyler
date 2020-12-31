# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import libsbml
import warnings
def biocrnpyler_supported_units():
    supported_units = {
        # Volume units
        'nL': 
        {'unit_kind':[libsbml.UNIT_KIND_LITRE],
        'unit_exponents':1,
        'unit_scale':-9,
        'unit_multiplier':1},

        'uL': 
        {'unit_kind':libsbml.UNIT_KIND_LITRE,
        'unit_exponents':[1],
        'unit_scale':-6,
        'unit_multiplier':1},

        'mL': 
        {'unit_kind':libsbml.UNIT_KIND_LITRE,
        'unit_exponents':1,
        'unit_scale':[-3],
        'unit_multiplier':1},

        'L': 
        {'unit_kind':libsbml.UNIT_KIND_LITRE,
        'unit_exponents':1,
        'unit_scale':0,
        'unit_multiplier':[1]},
        
        # Concentration units
        'M': 
        {'unit_kind':[libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
        'unit_exponents':[1, -1],
        'unit_scale':[0, 0],
        'unit_multiplier':[1,1]},

        'mM': 
        {'unit_kind':[libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
        'unit_exponents':[1,-1],
        'unit_scale':[-3,0],
        'unit_multiplier':[1,1]},

        'uM': 
        {'unit_kind':[libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
        'unit_exponents':[1, -1],
        'unit_scale':[-6, 0],
        'unit_multiplier':[1,1]},

        'nM': 
        {'unit_kind':[libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
        'unit_exponents':[1, -1],
        'unit_scale':[-9, 0],
        'unit_multiplier':[1,1]},

        'mM': 
        {'unit_kind':[libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITRE],
        'unit_exponents':[1, -1],
        'unit_scale':[-3, 0],
        'unit_multiplier':[1,1]},
        
        # Time units
        'hour': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[1],
        'unit_scale':[0],
        'unit_multiplier':[3600]},

        'minute': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[1],
        'unit_scale':[0],
        'unit_multiplier':[60]},

        'second': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[1],
        'unit_scale':[0],
        'unit_multiplier':[1]},
        
        # Common parameter units
        'per_second': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[-1],
        'unit_scale':[0],
        'unit_multiplier':[1]},

        'per_minute': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[-1],
        'unit_scale':[0],
        'unit_multiplier':[1/60]},

        'per_hour': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND],
        'unit_exponents':[-1],
        'unit_scale':[0],
        'unit_multiplier':[1/3600]},

        'mole_per_litre': 
        {'unit_kind':[libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE],
        'unit_exponents':[-1,1],
        'unit_scale':[0,0],
        'unit_multiplier':[1,1]},

        'litre_per_mole_per_second': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND, libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE],
        'unit_exponents':[-1,1,-1],
        'unit_scale':[0,0,0],
        'unit_multiplier':[1,1,1]},

        'litre_per_mole_per_hour': 
        {'unit_kind':[libsbml.UNIT_KIND_SECOND, libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE],
        'unit_exponents':[-1,1,-1],
        'unit_scale':[0,0,0],
        'unit_multiplier':[1/3600,1,1]},
    }
    ### Add your own units to this dictionary ###
    return supported_units

def create_new_unit_definition(model, unit_id):
    """
    Creates a new UnitDefinition inside the 
    SBML Model object passed as model argument and other attributes. 
    Returns a pointer to the new libSBML object created for the unit type.
    """
    supported_units = biocrnpyler_supported_units()
    if type(unit_id) is not str:
        raise ValueError(
            'The arguments are not of expected type. unit_id must be a string.')
    if unit_id not in supported_units.keys():
        warnings.warn('The string identifier for the unit {0} is not supported by BioCRNpyler. Add this to the dictionary in biocrnpyler/units.py if you want this unit.'.format(unit_id))
        return None
    unit_kind = supported_units[unit_id]['unit_kind']
    unit_exponents = supported_units[unit_id]['unit_exponents']
    unit_scale = supported_units[unit_id]['unit_scale']
    unit_multiplier = supported_units[unit_id]['unit_multiplier']

    unitdef = model.createUnitDefinition()
    unitdef.setId(unit_id)
    # Scale list
    if type(unit_kind) is list and type(unit_scale) is not list:
        unit_scaleList = []
        for kind in unit_kind:
            unit_scaleList.append(unit_scale)

    elif type(unit_kind) is list and type(unit_scale) is list:
        if len(unit_kind) != len(unit_scale):
            raise ValueError(
                'Lengths of unit_scale and unit kind lists are not equal')
        unit_scaleList = unit_scale[:]

    # Multiplier list
    if type(unit_kind) is list and type(unit_multiplier) is not list:
        unit_multiplierList = []
        for kind in unit_kind:
            unit_multiplierList.append(unit_multiplier)

    elif type(unit_kind) is list and type(unit_multiplier) is list:
        if len(unit_kind) != len(unit_multiplier):
            raise ValueError(
                'Lengths of unit kind and unit_multiplier lists are not equal')
        unit_multiplierList = unit_multiplier[:]

    if type(unit_kind) is not list:
        unit_kind = [unit_kind]
        if type(unit_scale) is list:
            unit_scale = unit_scale[0]
        if type(unit_multiplier) is list:
            unit_multiplier = unit_multiplier[0]
        if type(unit_exponents) is list:
            unit_exponents = unit_exponents[0]
        if type(unit_scale) is not int or type(unit_multiplier) is not int:
            raise ValueError(
                'Scale and unit_multiplier must be integers when there is only one unit kind')
        unit_scaleList = [unit_scale]
        unit_multiplierList = [unit_multiplier]
    if type(unit_exponents) is not list:
        if type(unit_exponents) is not int:
            raise ValueError('All unit_exponentss should be integers')
        unit_exponents = [unit_exponents]
    if len(unit_kind) != len(unit_exponents):
        raise ValueError(
            'Lengths of unit kind and unit unit_exponents lists must be equal')

    for kind, expo, unit_scale, unit_multiplier in zip(unit_kind, unit_exponents, unit_scaleList, unit_multiplierList):
        unit = unitdef.createUnit()
        unit.setKind(kind)
        unit.setExponent(expo)
        unit.setScale(unit_scale)
        unit.setMultiplier(unit_multiplier)
    return unitdef
