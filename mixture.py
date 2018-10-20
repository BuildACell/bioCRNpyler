# mixture.py - Mixture class and related functions
# RMM, 11 Aug 2018
#
# This file default the mixture class, which is used to hold a
# collection of components.  The docstring from the Mixture class
# describes the use of this function in more detail.
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.
from warnings import warn


import component as comp
import chemical_reaction_network as crn

"""Container for components (extract, genes, etc)

    The Mixture class is used as a container for a set of components
    that define the species and reactions to implement a TX-TL system.
    Mixtures can be added and scaled to create new mixtures (not yet
    implemented).

    Data attributes
    ---------------
    name                Name of the mixture
    components          List of components in the mixture (list of Components)
    concentrations      Concentration of each component (list of floats)
    default_mechanisms  Default mechanisms for this mixture (dict)
    custom_mechanisms   User-specified mechanisms for this mixture (dict)
    parameters          Global parameters for the mixture (dict)
"""
class Mixture():
    def __init__(self, name="", mechanisms={}, components = [], parameters = {}, default_mechanisms = {}, global_mechanisms = {}, **kwargs):
        "Create a new mixture"

        # Initialize instance variables
        self.name = name  # Save the name of the mixture
        self.parameters = parameters

        # Override the default mechanisms with anything we were passed
        self.default_mechanisms = default_mechanisms #default parameters are used by mixture subclasses
        self.custom_mechanisms = mechanisms

        #Mechaninisms stores the mechanisms used for compilation where defaults are overwritten by custom mechanisms
        self.mechanisms = dict(self.default_mechanisms)

        if isinstance(self.custom_mechanisms, dict):
            for mech_type in self.custom_mechanisms:
                self.mechanisms[mech_type] = self.custom_mechanisms[mech_type]
        elif isinstance(self.custom_mechanisms, list):
            for mech in self.custom_mechanisms:
                self.mechanisms[mech.type] = mech
        else:
            raise ValueError("Mechanisms must be passed as a list of instantiated objects or a dictionary {type:mechanism}")

        #Global mechanisms are applied just once ALL species generated from components inside a mixture
        #Global mechanisms should be used rarely, and with care. An example usecase is degredation via dilution.
        self.global_mechanisms = global_mechanisms

        self.components = []  # components contained in mixture
        self.added_species = [] # if chemcial_reaction_network.specie objects are passed in as components they are stored here
        for component in components:
            if isinstance(component, comp.Component):
                self.add_components(component)
            elif isinstance(component, crn.specie):
                self.added_species += component
            else:
                raise ValueError("Objects passed into mixture as Components must be of the class Component or chemical_reaction_network.specie")

    def add_components(self, components):
        if isinstance(components, comp.Component):
            components = [components]

        for component in components:
            if isinstance(component, comp.Component):
                self.components.append(component)
                component.update_mechanisms(mixture_mechanisms = self.mechanisms)
                component.update_parameters(mixture_parameters=self.parameters)
            else:
                warn("Non-component added to mixture "+self.name, RuntimeWarning)

    def update_species(self):
        self.crn_species = self.added_species
        for component in self.components:
            new_species = component.update_species()
            self.crn_species += component.update_species()

        #Update Global Mechanisms
        for mech in self.global_mechanisms:
            self.crn_species += self.global_mechanisms[mech].update_species_global(self.crn_species, self.parameters)

        return self.crn_species

    def update_reactions(self):
        if self.crn_species == None:
            raise AttributeError("Mixture.crn_species not defined. mixture.update_species() must be called before mixture.update_reactions()")

        self.crn_reactions = []
        for component in self.components:
            self.crn_reactions += component.update_reactions()

        #update with global mechanisms
        for mech in self.global_mechanisms:
            self.crn_reactions += self.global_mechanisms[mech].update_reactions_global(self.crn_species, self.parameters)
        return self.crn_reactions

    def compile_crn(self):
        species = self.update_species()
        reactions = self.update_reactions()
        CRN = crn.chemical_reaction_network(species, reactions)
        return CRN

    def __str__(self):
        return type(self).__name__ +': '+self.name

    def __repr__(self):
        txt = str(self)+"\n"
        txt += "Components = ["
        for comp in self.components:
            txt+="\n\t"+str(comp)
        txt+=" ]\nMechanisms = {"
        for mech in self.mechanisms:
            txt+="\n\t"+mech+":"+self.mechanisms[mech].name
        txt+=" }\nGlobal Mechanisms = {"
        for mech in self.global_mechanisms:
            txt+="\n\t"+mech+":"+self.mechanisms[mech].name
        txt+=" }"
        return txt
