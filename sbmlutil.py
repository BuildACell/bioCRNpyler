# sbmlutil.py - libsbml helper functions
# RMM, 14 Aug 2018
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import libsbml
import re
from warnings import warn

# Reaction ID number (global)
reaction_id = 0

# Create an SBML model
def create_sbml_model():
    document = libsbml.SBMLDocument(3, 1)
    model = document.createModel()

    # Define units for area (not used, but keeps COPASI from complaining)
    unitdef = model.createUnitDefinition()
    unitdef.setId('square_metre')
    unit = unitdef.createUnit()
    unit.setKind(libsbml.UNIT_KIND_METRE)
    unit.setExponent(2)
    unit.setScale(0)
    unit.setMultiplier(1)

    # Set up required units and containers
    model.setTimeUnits('second')            # set model-wide time units
    model.setExtentUnits('mole')            # set model units of extent
    model.setSubstanceUnits('mole')         # set model substance units
    model.setLengthUnits('metre')           # area units (never used?)
    model.setAreaUnits('square_metre')      # area units (never used?)
    model.setVolumeUnits('litre')           # default volume unit

    # Define the default compartment
    compartment = model.createCompartment()
    compartment.setId('txtl')
    compartment.setConstant(True)           # keep compartment size constant
    compartment.setSpatialDimensions(3)     # 3 dimensional compartment
    compartment.setVolume(1e-6)             # 1 microliter

    return document, model, compartment

# Helper function to add a species to the model
def add_species(mixture, type, name, ic=None, debug=False):
    model = mixture.model   # Get the model where we will store results
    
    # Construct the species name
    prefix = type + " " if type is not None else ""
    species_name = prefix + name
    
    # Construct the species ID
    species_id = _id_from_name(species_name)
    
    # Check to see if this species is already present
    species = find_species(mixture, species_id)
    if species is None:
        if debug: print("Adding species %s" % species_name)
        species = model.createSpecies()
        species.setName(species_name)
        species.setId(species_id)
        species.setCompartment(mixture.compartment.getId())
        species.setConstant(False)
        species.setBoundaryCondition(False)
        species.setHasOnlySubstanceUnits(False)

    else:
        if debug: print("add_species: %s already exists", species.getId())
        
    # Set the initial concentration (if specified)
    #! TODO: Decide whether to warn if species is already present
    #! TODO: add initial concentrations if species is already present
    if ic != None:
        if debug: print("    %s IC = %s" % (species_name, ic))
        species.setInitialConcentration(float(ic))

    return species

# Look for a species in the current mixture
def find_species(mixture, species_name):
    model = mixture.model   # Get the model where we will store results
    
    # Construct the species ID (no-op if already a species ID)
    species_id = _id_from_name(species_name)

    #! TODO: Add error checking
    return model.getSpecies(species_id)

# Helper function to add a parameter to the model
def add_parameter(mixture, name, value=0, debug=False):
    model = mixture.model   # Get the model where we will store results

    # Check to see if this parameter is already present
    parameter = find_parameter(mixture, name)   #! TODO: add error checking
    if parameter is None:
        if debug: print("Adding parameter %s" % name)
        parameter = model.createParameter()
        parameter.setId(name)                   #! TODO: add error checking

    else:
        if debug: print("add_parameter: %s already exists", parameter.getId())

    # Set the value of the parameter
    parameter.setValue(float(value))
    parameter.setConstant(True)

    return parameter

# Look for a parameter in the current model
def find_parameter(mixture, id):
    model = mixture.model               # Get model where parameters are stored
    return model.getParameter(id)       #! TODO: add error checking

# Helper function to add a reaction to a model
#! Add stochiometry argument to allow non-unitary stochiometries
def add_reaction(mixture, reactants, products, kf, kr=None, id=None,
                 parameters={}, prefix="r", debug=False):
    """Add a reaction to a model

    The `add_reaction` function is used to add a reaction to a model.
    It allows specifation for reaction rates as either symbolic or
    numeric entries, including the use of global symbolic names.
    Reactions can be unidirectional or reversible.

    Parameters
    ----------
    model       SBML model
    reactants   List of SBML species that are reactants in the reaction
    projects    List of SBML species that are products of the reaction
    kf          Forward rate constant (parameter, string, number, or list)
    kf          Reverse rate constant (None if non-reversible)
    id          Optional parameter to specify reaction id (otherwise numbered)

    Note: the current implementation requires that non-unitary
    stochiometries be represented by repeated entries in the reactants
    and/or products list.

    """
    model = mixture.model  # Get the model where we will store results

    # Create the reaction
    reaction = model.createReaction()
    reaction.setReversible(False)
    reaction.setFast(False)

    # Store the reaction id
    global reaction_id
    reaction.setId("%s%d" % (prefix, reaction_id));
    reaction_id += 1

    if debug: print("Creating reaction: ",
                    reactants, "-->[", kf, "] ", products)

    #! Sort out the reaction rates and parameter names
    #
    # Reactions can be specified in multiple forms
    #
    #   "string"                Reaction is a global paramater
    #   ["string", float]       Local reaction rate, with initial value
    #   float                   Local reaction rate

    #
    # Create forward and reverse rate strings
    #
    # The rate expression for the kinetic law in the expression is
    # created by building up the string that specifies the kinetic
    # law.  If the parameter `kf` is a string, we assume it is a
    # global parameter value that will be set later.  if it is a
    # number or a Parameter, then we create a local parameter within
    # this reaction.
    #
    if isinstance(kf, (float, int)):
        kfname = "k"
    elif isinstance(kf, Parameter):
        kfname = kf.name
    elif isinstance(kf, str):
        kfname = kf
    else:
        raise TypeError("add reaction: unknown parameter type", kf)

    # Create the reactants
    ratestring = kfname
    for species in reactants:
        reactant = reaction.createReactant()
        reactant.setSpecies(species.getId())    #! TODO: add error checking
        reactant.setConstant(True)
        ratestring += " * " + species.getId()

    # Create the products
    for species in products:
        product = reaction.createProduct()
        product.setSpecies(species.getId())     #! TODO: add error checking
        product.setConstant(True)

    # Create a kinetic law for the reaction
    if debug: print("    Creating kinetic law (%s): %s" %
                    (reaction.getId(), ratestring))
    ratelaw = reaction.createKineticLaw();
    if isinstance(kf, Parameter):
        param = ratelaw.createParameter();
        param.setId(kf.name)
        
        # Set the parameter value
        if kf.type == 'Numeric':
            param.setValue(kf.value)
        elif kf.type == 'Expression':
            #! TODO: handle more general expressions
            param.setValue(float(eval(kf.value)))
        else:
            warn("add_reaction: parameter type %s not supported" % kf.type)
            
    elif isinstance(kf, (float, int)):
        param = ratelaw.createParameter();
        param.setId(kfname)
        param.setValue(float(kf))
    ratelaw.setFormula(ratestring);

    # If the reverse rate is given, switch things around create reverse reaction
    if kr is not None:
        revreaction = add_reaction(mixture, products, reactants, kr, None,
                                   prefix=prefix)
        return reaction, revreaction

    return reaction

# Utility function to convert name to id
def _id_from_name(name):
    "Convert name to a id (remove spaces and other characters)"
    id = re.sub(" ", "_", name)
    id = re.sub("--", "_", id)
    id = re.sub(":", "_", id)
    return id
