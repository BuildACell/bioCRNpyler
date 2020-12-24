import libsbml

# Volume units


def create_unit_nL(model):
    nL = create_new_unit_definition(
        model, unit_id='nL', unit_kind=libsbml.UNIT_KIND_LITER, unit_exponents=1, unit_scale=-9)
    return nL


def create_unit_uL(model):
    uL = create_new_unit_definition(
        model, unit_id='uL', unit_kind=libsbml.UNIT_KIND_LITER, unit_exponents=1, unit_scale=-6)
    return uL


def create_unit_mL(model):
    mL = create_new_unit_definition(
        model, unit_id='mL', unit_kind=libsbml.UNIT_KIND_LITER, unit_exponents=1, unit_scale=-3)
    return mL


def create_unit_L(model):
    L = create_new_unit_definition(
        model, unit_id='L', unit_kind=libsbml.UNIT_KIND_LITER, unit_exponents=1, unit_scale=0)
    return L

# Concentration units


def create_unit_M(model):
    M = create_new_unit_definition(model, unit_id='M', unit_kind=[
        libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITER], unit_exponents=[1, -1], unit_scale=[0, 0])
    return M


def create_unit_mM(model):
    mM = create_new_unit_definition(model, unit_id='mM', unit_kind=[
                                    libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITER], unit_exponents=[1, -1], unit_scale=[-3, 0])
    return mM


def create_unit_uM(model):
    uM = create_new_unit_definition(model, unit_id='uM', unit_kind=[
                                    libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITER], unit_exponents=[1, -1], unit_scale=[-6, 0])
    return uM


def create_unit_nM(model):
    nM = create_new_unit_definition(model, unit_id='nM', unit_kind=[
                                    libsbml.UNIT_KIND_MOLE, libsbml.UNIT_KIND_LITER], unit_exponents=[1, -1], unit_scale=[-9, 0])
    return nM

# Time units


def create_unit_hour(model):
    hour = create_new_unit_definition(
        model, unit_id='hour', unit_kind=libsbml.UNIT_KIND_SECOND, unit_exponents=1, unit_scale=0, unit_multipler=3600)
    return hour


def create_unit_minute(model):
    minute = create_new_unit_definition(
        model, unit_id='minute', unit_kind=libsbml.UNIT_KIND_SECOND, unit_exponents=1, unit_scale=0, unit_multipler=60)
    return minute


def create_unit_second(model):
    second = create_new_unit_definition(
        model, unit_id='second', unit_kind=libsbml.UNIT_KIND_SECOND, unit_exponents=1, unit_scale=0, unit_multipler=1)
    return second

# Common parameter units

def create_unit_per_second(model):
    per_second = create_new_unit_definition(
        model, unit_id='per_second', unit_kind=libsbml.UNIT_KIND_SECOND, unit_exponents=-1, unit_scale=0, unit_multipler=1)
    return per_second

def create_unit_per_hour(model):
    per_hour = create_new_unit_definition(
        model, unit_id='per_hour', unit_kind=libsbml.UNIT_KIND_SECOND, unit_exponents=-1, unit_scale=0, unit_multipler=3600)
    return per_hour

def create_unit_mole_per_litre(model):
    mole_per_litre = create_new_unit_definition(
        model, unit_id='mole_per_litre', unit_kind=[libsbml.UNIT_KIND_SECOND, libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE], 
        unit_exponents=[-1,1,-1], unit_scale=[0,0,0], unit_multipler=[1,1,1])
    return mole_per_litre

def create_unit_litre_per_mole_per_second(model):
    litre_per_mole_per_second = create_new_unit_definition(
        model, unit_id='litre_per_mole_per_second', unit_kind=[libsbml.UNIT_KIND_SECOND, libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE], 
        unit_exponents=[-1,1,-1], unit_scale=[0,0,0], unit_multipler=[1,1,1])
    return litre_per_mole_per_second

def create_unit_litre_per_mole_per_hour(model):
    litre_per_mole_per_hour = create_new_unit_definition(
        model, unit_id='litre_per_mole_per_hour', unit_kind=[libsbml.UNIT_KIND_SECOND, libsbml.UNIT_KIND_LITRE, libsbml.UNIT_KIND_MOLE], 
        unit_exponents=[-1,1,-1], unit_scale=[0,0,0], unit_multipler=[3600,1,1])
    return litre_per_mole_per_hour
#
### Add your own units here ###


def create_new_unit_definition(model, unit_id, unit_kind, unit_exponents, unit_scale=0, unit_multipler=1):
    """
    Creates a new UnitDefinition inside the 
    SBML Model object passed as model argument and other attributes. 
    Returns a pointer to the new libSBML object created for the unit type.
    """
    if type(unit_id) is not str:
        raise ValueError(
            'The arguments are not of expected type. unit_id must be a string of valid SId format')

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
    if type(unit_kind) is list and type(unit_multipler) is not list:
        unit_multiplerList = []
        for kind in unit_kind:
            unit_multiplerList.append(unit_multipler)

    elif type(unit_kind) is list and type(unit_multipler) is list:
        if len(unit_kind) != len(unit_multipler):
            raise ValueError(
                'Lengths of unit kind and unit_multipler lists are not equal')
        unit_multiplerList = unit_multipler[:]

    if type(unit_kind) is not list:
        unit_kind = [unit_kind]
        if type(unit_scale) is not int or type(unit_multipler) is not int:
            raise ValueError(
                'Scale and unit_multipler must be integers when there is only one unit kind')
        unit_scaleList = [unit_scale]
        unit_multiplerList = [unit_multipler]
    if type(unit_exponents) is not list:
        if type(unit_exponents) is not int:
            raise ValueError('All unit_exponentss should be integers')
        unit_exponents = [unit_exponents]
    if len(unit_kind) != len(unit_exponents):
        raise ValueError(
            'Lengths of unit kind and unit unit_exponents lists must be equal')

    for kind, expo, unit_scale, unit_multipler in zip(unit_kind, unit_exponents, unit_scaleList, unit_multiplerList):
        unit = unitdef.createUnit()
        unit.setKind(kind)
        unit.setExponent(expo)
        unit.setScale(unit_scale)
        unit.setMultiplier(unit_multipler)
    return unitdef
