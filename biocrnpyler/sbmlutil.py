# sbmlutil.py - libsbml helper functions
# RMM, 14 Aug 2018
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import libsbml
import numpy as np
from warnings import warn

# Reaction ID number (global)
reaction_id = 0


# Create an SBML model
def create_sbml_model(compartment_id="default", time_units='second', extent_units='mole', substance_units='mole',
                      length_units='metre', area_units='square_metre', volume_units='litre', volume = 1e-6, model_id = None):
    '''
    Creates an SBML Level 3 Version 2 model with some fixed standard settings.
    Returns the SBMLDocument and the Model object as a tuple.
    Refer to python-libsbml for more information on SBML API.
    '''
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    if model_id is None:
        model_id = 'biocrnpyler_'+str(np.random.randint(1e6))
    model.setId(model_id)
    model.setName(model_id)
    # Define units for area (not used, but keeps COPASI from complaining)
    unitdef = model.createUnitDefinition()
    unitdef.setId('square_metre')
    unit = unitdef.createUnit()
    unit.setKind(libsbml.UNIT_KIND_METRE)
    unit.setExponent(2)
    unit.setScale(0)
    unit.setMultiplier(1)

    # Set up required units and containers
    model.setTimeUnits(time_units)  # set model-wide time units
    model.setExtentUnits(extent_units)  # set model units of extent
    model.setSubstanceUnits(substance_units)  # set model substance units
    model.setLengthUnits(length_units)  # area units (never used?)
    model.setAreaUnits(area_units)  # area units (never used?)
    model.setVolumeUnits(volume_units)  # default volume unit

    # Define the default compartment
    compartment = model.createCompartment()
    compartment.setId(compartment_id)
    compartment.setName(compartment_id)
    compartment.setConstant(True)  # keep compartment size constant
    compartment.setSpatialDimensions(3)  # 3 dimensional compartment
    compartment.setVolume(volume)  # 1 microliter

    # Returning document is enough. document.getModel() gives the model, and model.getCompartment(0) gives the compartment.
    return document, model


# Creates an SBML id from a chemical_reaction_network.species object
def species_sbml_id(species, document=None):
    # Construct the species ID
    all_ids = []
    if document:
        all_ids = getAllIds(document.getListOfAllElements())
    trans = SetIdFromNames(all_ids)
    species_id = trans.getValidIdForName(repr(species))
    return species_id


# Helper function to add a species to the model
# species must be chemical_reaction_network.species objects
def add_species(model, compartment, species, debug=False, initial_concentration=None):
    model = model  # Get the model where we will store results

    # Construct the species name
    species_name = repr(species)

    # Construct the species ID
    species_id = species_sbml_id(species, model.getSBMLDocument())

    if debug: print("Adding species", species_name, species_id)
    sbml_species = model.createSpecies()
    sbml_species.setName(species_name)
    sbml_species.setId(species_id)
    sbml_species.setCompartment(compartment.getId())
    sbml_species.setConstant(False)
    sbml_species.setBoundaryCondition(False)
    sbml_species.setHasOnlySubstanceUnits(False)
    if initial_concentration is None:
        initial_concentration = 0
    sbml_species.setInitialConcentration(initial_concentration)

    return sbml_species


# Helper function to add a parameter to the model
def add_parameter(mixture, name, value=0, debug=False):
    model = mixture.model  # Get the model where we will store results

    # Check to see if this parameter is already present
    parameter = find_parameter(mixture, name)  # ! TODO: add error checking
    if parameter is None:
        if debug: print("Adding parameter %s" % name)
        parameter = model.createParameter()
        all_ids = getAllIds(model.getSBMLDocument().getListOfAllElements())
        trans = SetIdFromNames(all_ids)
        parameter.setId(trans.getValidIdForName(name))  # ! TODO: add error checking

    else:
        if debug: print("add_parameter: %s already exists", parameter.getId())

    # Set the value of the parameter
    parameter.setValue(float(value))
    parameter.setConstant(True)

    return parameter


# Look for a parameter in the current model
def find_parameter(mixture, id):
    model = mixture.model  # Get model where parameters are stored
    return model.getParameter(id)  # ! TODO: add error checking


# Helper function to add a reaction to a model
# reaction must be a chemical_reaction_network.reaction object
#propensity params is a dictionary of the parameters for non-massaction propensities.
def add_reaction(model, inputs, input_coefs, outputs, output_coefs,
                 reaction_id, k = None, kname = None,
                 stochastic = False, propensity_type = "massaction", 
                 propensity_params = None, propensity_annotation = True):

    # Create the reaction
    reaction = model.createReaction()
    reaction.setReversible(False)
    # reaction.setFast(False) # Deprecated in SBML
    all_ids = getAllIds(model.getSBMLDocument().getListOfAllElements())
    trans = SetIdFromNames(all_ids)
    reaction.setId(trans.getValidIdForName(reaction_id))
    reaction.setName(reaction.getId())
    ratestring = "" #Stores the string representing the rate function
    annotation_dict = {"type":propensity_type}
    # Create a kinetic law for the reaction
    ratelaw = reaction.createKineticLaw()
    #Create Local Propensity Parameters
    if propensity_type=="massaction":
        if kname is None:
            kname = "k"
        
        param = ratelaw.createParameter()
        param.setId(kname)
        param.setConstant(True)
        if k is not None and propensity_params is None:
            param.setValue(k)
            annotation_dict["k"] = k
        elif 'k' in propensity_params:
            param.setValue(propensity_params['k'])
            annotation_dict["k"] = propensity_params['k']
        elif k is not None and "k" in propensity_params and propensity_params['k'] != k:
            raise ValueError("Keyword k and propensity_params['k'] have different values. Only one of these arguments is needed or they must match.")
        else:
            raise ValueError("Massaction propensities require a rate k which can be passed into add_reaction as a keyword k= or inside the propensity_params keyword dictionary.")
        
        ratestring = kname

    #Hill Function Propensities
    elif propensity_type in ["hillpositive", "hillnegative", "proportionalhillpositive", "proportionalhillnegative"]:
        if not ("k" in propensity_params and "K" in propensity_params and "n" in propensity_params):
            raise ValueError(propensity_type+" requires the following keys in the propensity_params dictionary: "
                    "'k':rate constant (float)"
                    "'n':cooperativity(float), "
                    "and 'K':dissociationc constant (float).")
        param_k = ratelaw.createParameter()
        param_k.setId("k")
        param_k.setConstant(True)
        param_k.setValue(propensity_params['k'])

        param_n = ratelaw.createParameter()
        param_n.setId("n")
        param_n.setConstant(True)
        param_n.setValue(propensity_params['n'])

        param_K = ratelaw.createParameter()
        param_K.setId("K")
        param_K.setConstant(True)
        param_K.setValue(propensity_params['K'])

        ratestring = "k"

        annotation_dict["k"] = propensity_params['k']
        annotation_dict["K"] = propensity_params['K']
        annotation_dict["n"] = propensity_params['n']


    elif propensity_type == "general":
        raise NotImplementedError("SBML writing of general propensities not implemented")
    else:
        raise ValueError(propensity_type+" is not a supported propensity_type")

    # Create the reactants
    for i in range(len(inputs)):
        species = str(inputs[i]).replace("'", "")
        stoichiometry = input_coefs[i]
        # species_id = species_sbml_id(species, model.getSBMLDocument())

        # What to do when there are multiple species with same name?
        species_id = getSpeciesByName(model,species).getId()
        reactant = reaction.createReactant()
        reactant.setSpecies(species_id)  # ! TODO: add error checking
        reactant.setConstant(False)
        if stoichiometry is None or stoichiometry is np.nan:
            stoichiometry = 1.0
        reactant.setStoichiometry(stoichiometry)
        #Create Rate-strings for massaction propensities
        if propensity_type=="massaction" and stochastic:
            for i in range(stoichiometry):
                if i > 0:
                    ratestring += f" * ( {species_id} - {i} )"
                else:
                    ratestring += f" * {species_id}"

        elif propensity_type=="massaction" and not stochastic:
            if stoichiometry > 1:
                ratestring += f" * {species_id}^{stoichiometry}"
            else:
                ratestring += f" * {species_id}"

    #Create ratestring for non-massaction propensities
    if propensity_type == "hillpositive":
        if not ("s1" in propensity_params):
            raise ValueError("hillpositive propensities, p(s1; k, K, n) "
                    "= k*s1^n/(s1^n + K), require the following key in the propensity_params dictionary:"
                    "'s1':species (chemical_reaction_network.species)")

        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()

        ratestring+=f"*{s_species_id}^n/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id

    elif propensity_type == "hillnegative":
        if not ("s1" in propensity_params):
            raise ValueError("hillnegative propensities, "
                    "p(s1; k, K, n) = k*1/(s1^n + K), require the following key in the propensity_params dictionary:"
                    "'s1':species (chemical_reaction_network.species)")
        s = str(propensity_params['s1']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()

        ratestring+=f"/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id

    elif propensity_type == "proportionalhillpositive":
        if not ("s1" in propensity_params and "d" in propensity_params):
            raise ValueError("proportionalhillpositive propensities, "
                "p(s1, d; k, K, n) = k*d*s1^n/(s1^n + K), require the following key in the propensity_params dictionary:"
                "'s1':species (chemical_reaction_network.species)"
                "'d':species (chemical_reaction_network.species), ")

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        d_species_id = getSpeciesByName(model,d).getId()

        ratestring+=f"*{d_species_id}*{s_species_id}^n/({s_species_id}^n + K)"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id

    elif propensity_type == "proportionalhillnegative":
        if not ("s1" in propensity_params and "d" in propensity_params):
            raise ValueError("proportionalhillnegative propensities, "
                "p(s1, d; k, K, n) = k*d/(s1^n + K), require the following key in the propensity_params dictionary:"
                "'s1':species (chemical_reaction_network.species)"
                "'d':species (chemical_reaction_network.species), ")

        s = str(propensity_params['s1']).replace("'", "")
        d = str(propensity_params['d']).replace("'", "")
        s_species_id = getSpeciesByName(model,s).getId()
        d_species_id = getSpeciesByName(model,d).getId()

        ratestring+=f"*{d_species_id}/({s_species_id}^n+K)"

        annotation_dict["s1"] = s_species_id
        annotation_dict["d"] = d_species_id

    elif propensity_type == "general":
        raise NotImplementedError("General propensity SBML Writing Not Implemented")


    # Create the products
    for i in range(len(outputs)):
        species = str(outputs[i]).replace("'", "")
        stoichiometry = output_coefs[i]
        product = reaction.createProduct()
        species_id = getSpeciesByName(model, species).getId()
        product.setSpecies(species_id)
        if stoichiometry is None or stoichiometry is np.nan:
            stoichiometry = 1.0
        product.setStoichiometry(stoichiometry)
        product.setConstant(False)

    # Set the ratelaw to the ratestring
    math_ast = libsbml.parseL3Formula(ratestring)
    ratelaw.setMath(math_ast)
    if propensity_annotation:
        annotation_string = "<PropensityType>"
        for k in annotation_dict:
            annotation_string += " "+k + "=" + str(annotation_dict[k])
        annotation_string += "</PropensityType>"
        reaction.appendAnnotation(annotation_string)

    return reaction


# !/usr/bin/env python
##
## @file    setIdFromNames.py
## @brief   Utility program, renaming all SIds that also has
##          names specified. The new id will be derived from
##          the name, with all invalid characters removed.
##
## @author  Frank T. Bergmann
##
##
## <!--------------------------------------------------------------------------
## This sample program is distributed under a different license than the rest
## of libSBML.  This program uses the open-source MIT license, as follows:
##
## Copyright (c) 2013-2017 by the California Institute of Technology
## (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
## and the University of Heidelberg (Germany), with support from the National
## Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
## Permission is hereby granted, free of charge, to any person obtaining a
## copy of this software and associated documentation files (the "Software"),
## to deal in the Software without restriction, including without limitation
## the rights to use, copy, modify, merge, publish, distribute, sublicense,
## and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
##
## The above copyright notice and this permission notice shall be included in
## all copies or substantial portions of the Software.
##
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
## THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
## DEALINGS IN THE SOFTWARE.
##
## Neither the name of the California Institute of Technology (Caltech), nor
## of the European Bioinformatics Institute (EMBL-EBI), nor of the University
## of Heidelberg, nor the names of any contributors, may be used to endorse
## or promote products derived from this software without specific prior
## written permission.
## ------------------------------------------------------------------------ -->
##
##

import sys
import os.path
import time


# This class implements an identifier transformer, that means it can be used
# to rename all sbase elements.
class SetIdFromNames(libsbml.IdentifierTransformer):
    def __init__(self, ids):
        # call the constructor of the base class
        libsbml.IdentifierTransformer.__init__(self)
        # remember existing ids ...
        self.existingIds = ids

        # The function actually doing the transforming. This function is called

    # once for each SBase element in the model.
    def transform(self, element):
        # return in case we don't have a valid element
        if (element is None \
           or element.getTypeCode() == libsbml.SBML_LOCAL_PARAMETER):
            return libsbml.LIBSBML_OPERATION_SUCCESS

            # or if there is nothing to do
        if (element.isSetName() == False \
           or element.getId() == element.getName()):
            return libsbml.LIBSBML_OPERATION_SUCCESS

            # find the new id
        newId = self.getValidIdForName(element.getName())

        # set it
        element.setId(newId)

        # remember it
        self.existingIds.append(newId)

        return libsbml.LIBSBML_OPERATION_SUCCESS

    def nameToSbmlId(self, name):
        IdStream = []
        count = 0
        end = len(name)

        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('x_')
        if '*' in name:
            IdStream.append('xx')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if (Id[len(Id) - 1] != '_'):
            return Id

        return Id[:-1]

    #
    # Generates the id out of the name, and ensures it is unique.
    # It does so by appending numbers to the original name.
    #
    def getValidIdForName(self, name):
        baseString = self.nameToSbmlId(name)
        id = baseString
        count = 1
        while (self.existingIds.count(id) != 0):
            id = "{0}_{1}".format(baseString, count)
            count = count + 1
        return id

    #  #


#  # Returns a list of all ids from the given list of elements
#  #
def getAllIds(allElements):
    result = []
    if allElements is None or allElements.getSize() == 0:
        return result

    for i in range(0, allElements.getSize()):
        current = allElements.get(i)
        if (current.isSetId() \
           and current.getTypeCode() != libsbml.SBML_LOCAL_PARAMETER):
            result.append(current.getId())
    return result


def getSpeciesByName(model, name, compartment=''):
    '''
    Returns a list of species in the Model with the given name
    compartment : (Optional) argument to specify the compartment name in which
    to look for the species.
    '''
    if type(name) is not str:
        raise ValueError('"name" must be a string.')
    species_found = []
    for species in model.getListOfSpecies():
        if species.getName() == name:
            if compartment != '':
                comp_elem = species.getCompartment()
                comp_name = model.getElementBySId(comp_elem).getName()
                if comp_name == compartment:
                    species_found.append(species)
                else:
                    continue
            else:
                species_found.append(species)

    if len(species_found) == 1:
        return species_found[0]
    elif not species_found:
        raise ValueError('The species ' + name + ' not found.')
    else:
        warn('Multiple species with name ' + name + ' found. Returning a list')
        return species_found
