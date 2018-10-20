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
def create_sbml_model(compartment_id = "default", time_units = 'second', extent_units = 'mole', substance_units = 'mole',
                      length_units = 'metre', area_units = 'square_metre', volume_units = 'litre'):
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
    model.setTimeUnits(time_units)            # set model-wide time units
    model.setExtentUnits(extent_units)            # set model units of extent
    model.setSubstanceUnits(substance_units)         # set model substance units
    model.setLengthUnits(length_units)           # area units (never used?)
    model.setAreaUnits(area_units)      # area units (never used?)
    model.setVolumeUnits(volume_units)           # default volume unit

    # Define the default compartment
    compartment = model.createCompartment()
    compartment.setId(compartment_id)
    compartment.setConstant(True)           # keep compartment size constant
    compartment.setSpatialDimensions(3)     # 3 dimensional compartment
    compartment.setVolume(1e-6)             # 1 microliter

    return document, model, compartment


#Creates an SBML id from a chemical_reaction_network.specie object
def specie_sbml_id(specie):
    # Construct the species ID
    specie_id = repr(specie).replace(" ", "_").replace(":", "_").replace("--", "_").replace("-", "_").replace("'", "")
    return specie_id

# Helper function to add a species to the model
# species must be chemical_reaction_network.species objects
def add_species(model, compartment, specie, debug=False):
    model = model   # Get the model where we will store results
    
    # Construct the species name
    specie_name = repr(specie)
    
    # Construct the species ID
    specie_id = specie_sbml_id(specie)

    if debug: print("Adding species", specie_name, specie_id)
    sbml_species = model.createSpecies()
    sbml_species.setName(specie_name)
    sbml_species.setId(specie_id)
    sbml_species.setCompartment(compartment.getId())
    sbml_species.setConstant(False)
    sbml_species.setBoundaryCondition(False)
    sbml_species.setHasOnlySubstanceUnits(False)

    return sbml_species

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
# reaction must be a chemical_reaction_network.reaction object
def add_reaction(model, inputs, input_coefs, outputs, output_coefs, k, reaction_id, kname = None,
                stochastic = False, mass_action = True):

    # Create the reaction
    reaction = model.createReaction()
    reaction.setReversible(False)
    reaction.setFast(False)
    reaction.setId(reaction_id)

    if kname == None:
        kname = "k"
        ratestring = kname

    if mass_action:
        # Create a kinetic law for the reaction
            ratelaw = reaction.createKineticLaw()
            param = ratelaw.createParameter()
            param.setId(kname)
            param.setConstant(True)
            param.setValue(k)
    else:
        raise NotImplementedError("SBML Writing of non-massaction ratelaws not implemented yet")


    # Create the reactants
    for i in range(len(inputs)):
        specie = str(inputs[i]).replace("'", "")
        stoichiometry = input_coefs[i]
        specie_id = specie_sbml_id(specie)
        reactant = reaction.createReactant()
        reactant.setSpecies(specie_id)    #! TODO: add error checking
        reactant.setConstant(True)
        reactant.setStoichiometry(stoichiometry)

        if mass_action and stochastic:
            for i in range(stoichiometry):
                if i > 0:
                    ratestring += " * "+ "( "+specie_id+" - "+str(i)+" )"
                else:
                    ratestring += " * " + specie_id

        elif mass_action and not stochastic:
            ratestring += " * " + specie_id

    # Create the products
    for i in range(len(outputs)):
        specie = str(outputs[i]).replace("'", "")
        stoichiometry = output_coefs[i]
        product = reaction.createProduct()
        specie_id = specie_sbml_id(specie)
        product.setSpecies(specie_id)
        reactant.setStoichiometry(stoichiometry)
        product.setConstant(True)

    math_ast = libsbml.parseL3Formula(ratestring)

    #Set the ratelaw to the ratestring
    math_ast =  libsbml.parseL3Formula(ratestring)
    ratelaw.setMath(math_ast)
    return reaction

