# sbmlutil.py - libsbml helper functions
# RMM, 14 Aug 2018
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import logging
import warnings
from random import randint
from typing import List
from warnings import warn
from .general import parameter_to_value
from ..utils.units import create_new_unit_definition
import libsbml  # type: ignore

# Reaction ID number (global)
reaction_id = 0

logger = logging.getLogger(__name__)


def create_sbml_model(
    compartment_id="default",
    time_units="second",
    extent_units="mole",
    substance_units="mole",
    length_units="metre",
    area_units="square_metre",
    volume_units="litre",
    volume=1e-6,
    model_id=None,
    **kwargs,
):
    """Creates an SBML Level 3 Version 2 model with some fixed standard settings.
    Refer to python-libsbml for more information on SBML API.
    :param compartment_id:
    :param time_units:
    :param extent_units:
    :param substance_units:
    :param length_units:
    :param area_units:
    :param volume_units:
    :param volume:
    :param model_id:
    :return:  the SBMLDocument and the Model object as a tuple
    """
    document = libsbml.SBMLDocument(3, 2)
    model = document.createModel()
    if model_id is None:
        model_id = "biocrnpyler_" + str(randint(1, int(1e6)))
    model.setId(model_id)
    model.setName(model_id)
    # Define units for area (not used, but keeps COPASI from complaining)
    unitdef = model.createUnitDefinition()
    unitdef.setId("square_metre")
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

    return document, model


# Creates an SBML id from a chemical_reaction_network.species object
def valid_sbml_id(given_id, document=None):
    # Construct the species ID
    all_ids = []
    if document:
        all_ids = getAllIds(document.getListOfAllElements())
    trans = SetIdFromNames(all_ids)
    valid_id = trans.getValidIdForName(repr(given_id))
    return valid_id


def add_all_species(
    model, species: List, initial_condition_dictionary: dict, compartment=None, **kwargs
):
    """adds a list of Species to the SBML model.
    :param model: valid SBML model
    :param species: list of species to be added to the SBML model
    :param compartment: compartment id, if empty species go to the first compartment
    :param initial_concentration_dict: a dictionary s --> initial_concentration
    :return: None
    """

    elementlist = model.getSBMLDocument().getListOfAllElements()

    all_ids = getAllIds(elementlist)

    trans = SetIdFromNames(all_ids)

    all_ids = [str(s) for s in species]
    all_sbml_ids = [trans.getValidIdForName(species_id) for species_id in all_ids]

    for s_ind, s in enumerate(species):
        s_id = all_sbml_ids[s_ind]
        if compartment is None or s.compartment is not None:
            # If no compartment was passed in or if species (s) has its own compartment set:
            compartment = get_compartment_by_name(model, s.compartment.name)
            if compartment is None:
                compartment = add_compartment(model, s.compartment)
        if s in initial_condition_dictionary:
            initial_concentration = parameter_to_value(initial_condition_dictionary[s])
        else:
            initial_concentration = 0
        add_species(
            model=model,
            compartment=compartment,
            species_name=str(s),
            species_id=s_id,
            initial_concentration=initial_concentration,
        )


def add_species(
    model, compartment, species_name, species_id, initial_concentration=None, **kwargs
):
    """Helper function to add a species to the sbml model.
    :param model:
    :param compartment: a compartment in the SBML model
    :param species: must be chemical_reaction_network.species objects
    :param initial_concentration: initial concentration of the species in the SBML model
    :return: SBML species object
    """

    # Construct the species ID

    logger.debug(f"Adding species: {species_name}, id: {species_id}")
    sbml_species = model.createSpecies()
    sbml_species.setName(species_name)
    sbml_species.setId(species_id)
    sbml_species.setCompartment(compartment.getId())
    sbml_species.setConstant(False)
    sbml_species.setBoundaryCondition(False)
    sbml_species.setHasOnlySubstanceUnits(False)
    sbml_species.setSubstanceUnits("mole")
    if initial_concentration is None:
        initial_concentration = 0
    sbml_species.setInitialConcentration(initial_concentration)

    return sbml_species


def add_all_compartments(model, compartments: List, **keywords):
    """Adds the list of Compartment objects to the SBML model
    :param model: valid SBML model
    :param compartments: list of compartments to be added to the SBML model
    :return: None
    """
    for compartment in compartments:
        add_compartment(model=model, compartment=compartment, **keywords)


def add_compartment(model, compartment, **keywords):
    """Helper function to add a compartment to the SBML model.
    :param model: a valid SBML model
    :param compartment: a Compartment object
    :return: SBML compartment object
    """
    sbml_compartment = model.createCompartment()
    compartment_id = compartment.name
    sbml_compartment.setId(compartment_id)
    sbml_compartment.setName(compartment.name)
    sbml_compartment.setConstant(True)  # keep compartment size constant
    # For example, 3 dimensional compartment
    sbml_compartment.setSpatialDimensions(compartment.spatial_dimensions)
    sbml_compartment.setSize(compartment.size)  # For example, 1e-6 liter
    if compartment.unit is not None:
        sbml_compartment.setUnits(compartment.unit)
    return sbml_compartment


def get_compartment_by_name(model, compartment_name):
    """Helper function to find the SBML compartment object
    given the compartment name in the SBML file
    """
    for compartment in model.getListOfCompartments():
        if compartment.getName() == compartment_name:
            return compartment


# Helper function to add a parameter to the model


def add_parameter(mixture, name, value=0, debug=False):
    model = mixture.model  # Get the model where we will store results

    # Check to see if this parameter is already present
    parameter = find_parameter(mixture, name)  # ! TODO: add error checking
    if parameter is None:
        parameter = model.createParameter()
        all_ids = getAllIds(model.getSBMLDocument().getListOfAllElements())
        trans = SetIdFromNames(all_ids)
        # ! TODO: add error checking
        parameter.setId(trans.getValidIdForName(name))
    # Set the value of the parameter
    parameter.setValue(float(value))
    parameter.setConstant(True)

    return parameter


# Look for a parameter in the current model
def find_parameter(mixture, id):
    model = mixture.model  # Get model where parameters are stored
    return model.getParameter(id)  # ! TODO: add error checking


def add_all_reactions(model, reactions: List, stochastic=False, **kwargs):
    """adds a list of reactions to the SBML model.
    :param model: an sbml model created by create_sbml_model()
    :param reactions: list of Reactions
    :param stochastic: binary flag for stochastic models
    :return: None
    """

    elementlist = model.getSBMLDocument().getListOfAllElements()

    all_ids = getAllIds(elementlist)

    trans = SetIdFromNames(all_ids)

    all_ids = [f"r{i}" for i in range(len(reactions))]
    all_ids_rev = [f"r{i}rev" for i in range(len(reactions))]
    all_sbml_ids = [trans.getValidIdForName(reaction_id) for reaction_id in all_ids]
    all_sbml_ids_rev = [
        trans.getValidIdForName(reaction_id) for reaction_id in all_ids_rev
    ]

    for rxn_count, r in enumerate(reactions):

        add_reaction(
            model=model,
            crn_reaction=r,
            reaction_id=all_sbml_ids[rxn_count],
            stochastic=stochastic,
            **kwargs,
        )
        # Reversible reactions are always seperated into two seperate reactions
        if r.is_reversible:
            add_reaction(
                model=model,
                crn_reaction=r,
                reaction_id=all_sbml_ids_rev[rxn_count],
                stochastic=stochastic,
                reverse_reaction=True,
                **kwargs,
            )


def add_reaction(
    model,
    crn_reaction,
    reaction_id: str,
    stochastic: bool = False,
    reverse_reaction: bool = False,
    **kwargs,
):
    """adds a sbml_reaction to an sbml model.
    :param model: an sbml model created by create_sbml_model()
    :param crn_reaction: must be a chemical_reaction_network.reaction object
    :param reaction_id: unique id of the reaction
    :param stochastic: stochastic model flag
    :param reverse_reaction:
    :return: SBML Reaction object
    """
    # Create the sbml_reaction in SBML
    sbml_reaction = model.createReaction()

    sbml_reaction.setId(reaction_id)
    sbml_reaction.setName(sbml_reaction.getId())
    # all reactions are set to be non-reversible in BioCRNpyler because this is correct in deterministic and stochastic simulation.
    sbml_reaction.setReversible(False)
    # Create the reactants and products for the sbml_reaction
    if not reverse_reaction:
        _create_reactants(
            reactant_list=crn_reaction.inputs, sbml_reaction=sbml_reaction, model=model
        )
        _create_products(
            product_list=crn_reaction.outputs, sbml_reaction=sbml_reaction, model=model
        )
    else:
        _create_reactants(
            reactant_list=crn_reaction.outputs, sbml_reaction=sbml_reaction, model=model
        )
        _create_products(
            product_list=crn_reaction.inputs, sbml_reaction=sbml_reaction, model=model
        )

    # Create the kinetic law and corresponding local propensity parameters
    crn_reaction.propensity_type.create_kinetic_law(
        model=model,
        sbml_reaction=sbml_reaction,
        stochastic=stochastic,
        crn_reaction=crn_reaction,
        reverse_reaction=reverse_reaction,
        **kwargs,
    )

    # Create SpeciesModifierReference in SBML for species that are referred by the
    # KineticLaw but not in reactants or products
    _create_modifiers(
        crn_reaction=crn_reaction, sbml_reaction=sbml_reaction, model=model
    )
    return sbml_reaction


def _create_reactants(reactant_list, sbml_reaction, model):
    for input in reactant_list:
        # What to do when there are multiple species with same name?
        species_id = str(input.species)
        reactant = sbml_reaction.createReactant()
        reactant.setSpecies(species_id)
        reactant.setConstant(False)
        reactant.setStoichiometry(input.stoichiometry)


def _create_products(product_list, sbml_reaction, model):
    for output in product_list:
        species_id = str(output.species)
        product = sbml_reaction.createProduct()
        product.setSpecies(species_id)
        product.setStoichiometry(output.stoichiometry)
        product.setConstant(False)


def _create_modifiers(crn_reaction, sbml_reaction, model):
    reactants_list = [str(i.species) for i in crn_reaction.inputs]
    products_list = [str(i.species) for i in crn_reaction.outputs]
    modifier_species = [
        str(i) for i in crn_reaction.propensity_type.propensity_dict["species"].values()
    ]
    for modifier_id in modifier_species:
        modifier_id = str(modifier_id)
        if modifier_id not in reactants_list and modifier_id not in products_list:
            modifier = sbml_reaction.createModifier()
            modifier.setSpecies(modifier_id)


# Creates a local parameter SBML kinetic rate law


def _create_local_parameter(ratelaw, name, value, constant=True):
    param = ratelaw.createParameter()
    param.setId(name)
    param.setConstant(constant)
    param.setValue(value)
    return param


# Creates a global parameter SBML model


def _create_global_parameter(model, name, value, p_unit=None, constant=True):
    if p_unit is not None:
        if not isinstance(p_unit, str):
            raise ValueError("Units for a parameter must be passed as strings.")
        unit_added = False
        for unit_definition in model.getListOfUnitDefinitions():
            if unit_definition.getId() == p_unit:
                unit_added = True
        if not unit_added:
            unit_created = create_new_unit_definition(model, p_unit)
            model.addUnitDefinition(unit_created)

    if model.getParameter(name) is None:
        param = model.createParameter()
        param.setId(name)
        param.setConstant(constant)
        param.setValue(value)
        if p_unit is not None:
            param.setUnits(str(p_unit))
    else:
        param = model.getParameter(name)
        if p_unit is not None:
            param.setUnits(str(p_unit))
    return param


##
# @file    setIdFromNames.py
# @brief   Utility program, renaming all SIds that also has
# names specified. The new id will be derived from
# the name, with all invalid characters removed.
##
# @author  Frank T. Bergmann
##
##
# <!--------------------------------------------------------------------------
# This sample program is distributed under a different license than the rest
# of libSBML.  This program uses the open-source MIT license, as follows:
##
# Copyright (c) 2013-2017 by the California Institute of Technology
# (California, USA), the European Bioinformatics Institute (EMBL-EBI, UK)
# and the University of Heidelberg (Germany), with support from the National
# Institutes of Health (USA) under grant R01GM070923.  All rights reserved.
##
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
##
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
##
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
##
# Neither the name of the California Institute of Technology (Caltech), nor
# of the European Bioinformatics Institute (EMBL-EBI), nor of the University
# of Heidelberg, nor the names of any contributors, may be used to endorse
# or promote products derived from this software without specific prior
# written permission.
# ------------------------------------------------------------------------ -->
##
##

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
        if element is None or element.getTypeCode() == libsbml.SBML_LOCAL_PARAMETER:
            return libsbml.LIBSBML_OPERATION_SUCCESS

            # or if there is nothing to do
        if element.isSetName() == False or element.getId() == element.getName():
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

        if "0" <= name[count] and name[count] <= "9":
            IdStream.append("x_")
        if "*" in name:
            IdStream.append("xx")
        for count in range(0, end):
            if (
                ("0" <= name[count] and name[count] <= "9")
                or ("a" <= name[count] and name[count] <= "z")
                or ("A" <= name[count] and name[count] <= "Z")
            ):
                IdStream.append(name[count])
            else:
                IdStream.append("_")
        Id = "".join(IdStream)
        if Id[len(Id) - 1] != "_":
            return Id

        return Id
        # return Id[:-1] #this code was removing underscores at the end of ComplexSpecies needed for imbedded ComplexSpecies.

    #
    # Generates the id out of the name, and ensures it is unique.
    # It does so by appending numbers to the original name.
    #
    def getValidIdForName(self, name):
        baseString = self.nameToSbmlId(name)
        id = baseString
        count = 1
        while self.existingIds.count(id) != 0:
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
        if current.isSetId() and current.getTypeCode() != libsbml.SBML_LOCAL_PARAMETER:
            result.append(current.getId())
    return result


def getSpeciesByName(model, name, compartment=""):
    """
    Returns a list of species in the Model with the given name
    compartment : (Optional) argument to specify the compartment name in which
    to look for the species.
    """
    if type(name) is not str:
        raise ValueError(f'"name" must be a string. Received {name}.')
    species_found = []
    for species in model.getListOfSpecies():
        if species.getName() == name:
            if compartment != "":
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
        raise ValueError("The species " + name + " not found.")
    else:
        warn("Multiple species with name " + name + " found. Returning a list")
        return species_found


# Validate SBML


class validateSBML(object):
    """
    libSBML class to validate the generated SBML models
    ## @brief   Validates SBMLDocument
    ## @author  Akiya Jouraku (translated from libSBML C++ examples)
    ## @author  Ben Bornstein
    ## @author  Michael Hucka
    """

    def __init__(self, ucheck):
        self.reader = libsbml.SBMLReader()
        self.ucheck = ucheck

    def validate(self, sbml_document, print_results=False):
        """sbml_document: libSBML SBMLDocument object.
        print_results: Print toggle for validation warnings.
        """
        sbmlDoc = sbml_document
        errors = sbmlDoc.getNumErrors()
        if print_results:
            print(
                "Validating SBML model with ID: {0}...".format(
                    sbmlDoc.getModel().getId()
                )
            )
        seriousErrors = False

        numReadErr = 0
        numReadWarn = 0
        errMsgRead = ""

        if errors > 0:
            for i in range(errors):
                severity = sbmlDoc.getError(i).getSeverity()
                if (severity == libsbml.LIBSBML_SEV_ERROR) or (
                    severity == libsbml.LIBSBML_SEV_FATAL
                ):
                    seriousErrors = True
                    numReadErr += 1
                else:
                    numReadWarn += 1
            errMsgRead = sbmlDoc.getErrorLog().toString()

        # If serious errors are encountered while reading an SBML document, it
        # does not make sense to go on and do full consistency checking because
        # the model may be nonsense in the first place.

        numCCErr = 0
        numCCWarn = 0
        errMsgCC = ""

        if seriousErrors:
            errMsgRead += "Further consistency checking and validation aborted."
        else:
            sbmlDoc.setConsistencyChecks(
                libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, self.ucheck
            )
            failures = sbmlDoc.checkConsistency()
            if failures > 0:
                isinvalid = False
                for i in range(failures):
                    severity = sbmlDoc.getError(i).getSeverity()
                    if (severity == libsbml.LIBSBML_SEV_ERROR) or (
                        severity == libsbml.LIBSBML_SEV_FATAL
                    ):
                        numCCErr += 1
                        isinvalid = True
                    else:
                        numCCWarn += 1
                if isinvalid:
                    errMsgCC = sbmlDoc.getErrorLog().toString()
        if errMsgRead or errMsgCC:
            if print_results:
                print()
                print("===== validation error/warning messages =====\n")
            if errMsgRead:
                if print_results:
                    print(errMsgRead)
            if errMsgCC:
                if print_results:
                    print("*** consistency check ***\n")
                    print(errMsgCC)
        if not (numReadErr + numCCErr):
            print("Successful!")
        return numReadErr + numCCErr


def validate_sbml(sbml_document, enable_unit_check=False, print_results=True):
    """Validates the generated SBML model by using libSBML SBML validation code."""
    validator = validateSBML(enable_unit_check)
    validation_result = validator.validate(sbml_document, print_results=print_results)
    if validation_result > 0:
        warn(
            "SBML model invalid. Run with print_results = False to hide print statements"
        )
    return validation_result
