
# Copyright (c) 2019, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from warnings import resetwarnings

from .component import Component
from .chemical_reaction_network import ChemicalReactionNetwork, Species, Reaction
from .parameter import Parameter
from typing import List, Union


class Mixture(object):
    def __init__(self, name="", mechanisms={}, components = [], parameters = {},
                 parameter_file = None, default_mechanisms = {},
                 global_mechanisms = {}, default_components = [],
                 species = [],
                 parameter_warnings = None, **kwargs):
        """
        A Mixture object holds together all the components (DNA,Protein, etc), mechanisms (Transcription, Translation),
        and parameters related to the mixture itself (e.g. Transcription rate). Default components and mechanisms can be
        added as well as global mechanisms that impacts all species (e.g. cell growth).

        :param name: Name of the mixture
        :param mechanisms: Dictionary of mechanisms
        :param components: List of components in the mixture (list of Components)
        :param parameters: Dictionary of parameters (check parameters documentation for the keys)
        :param parameter_file: Parameters can be loaded from a parameter file
        :param default_mechanisms:
        :param global_mechanisms: dict of global mechanisms that impacts all species (e.g. cell growth)
        :param default_components:
        :param parameter_warnings: suppressing parameter related warnings
        """
        # Initialize instance variables
        self.name = name  # Save the name of the mixture

        self.parameters = Parameter.create_parameter_dictionary(parameters,
                                                      parameter_file)
        # Toggles whether parameter warnings are raised. if None (default) this
        # can be toggled component by component.
        self.parameter_warnings = parameter_warnings

        # Override the default mechanisms with anything we were passed
        # default parameters are used by mixture subclasses.
        self.default_mechanisms = default_mechanisms
        self.custom_mechanisms = mechanisms

        # Mechanisms stores the mechanisms used for compilation where defaults
        # are overwritten by custom mechanisms.
        self.mechanisms = dict(self.default_mechanisms)

        if isinstance(self.custom_mechanisms, dict):
            for mech_type in self.custom_mechanisms:
                self.mechanisms[mech_type] = self.custom_mechanisms[mech_type]
        elif isinstance(self.custom_mechanisms, list):
            for mech in self.custom_mechanisms:
                self.mechanisms[mech.type] = mech
        else:
            raise ValueError("Mechanisms must be passed as a list of "
                             "instantiated objects or a dictionary "
                             "{type:mechanism}")

        # Global mechanisms are applied just once ALL species generated from
        # components inside a mixture
        # Global mechanisms should be used rarely, and with care. An example
        # usecase is degradation via dilution.
        self.global_mechanisms = global_mechanisms

        self.components = []  # components contained in mixture
        # if chemical_reaction_network.species objects are passed in as
        # components they are stored here
        self.added_species = []
        # process the species
        self.add_species(species)
        # TODO find out why do we need default_components!
        # process the components
        self.add_components(components+default_components)

        # internal lists for the species and reactions
        self.crn_species = None
        self.crn_reactions = None

    def add_species(self, species: Union[List[Species], Species]):
        if not isinstance(species, list):
            species_list = [species]
        else:
            species_list = species

        assert all(isinstance(x, Species) for x in species_list), 'only Species type is accepted!'

        self.added_species += species_list

    def add_components(self, components):
        if not isinstance(components, list):
            components = [components]

        for component in components:
            assert isinstance(component, Component), \
                "the object: %s passed into mixture as component must be of the class Component" % str(component)
            self.components.append(component)
            component.update_mechanisms(mixture_mechanisms=self.mechanisms)
            component.update_parameters(mixture_parameters=self.parameters)
            if self.parameter_warnings is not None:
                component.set_parameter_warnings(self.parameter_warnings)

    def update_species(self) -> List[Species]:
        """ it generates the list of species based on all the mechanisms and global mechanisms

        :return: list of species generated by all the mechanisms and global mechanisms
        """
        # TODO check if we can merge the two variables
        self.crn_species = self.added_species
        for component in self.components:
            self.crn_species += component.update_species()

        # Update Global Mechanisms
        for mech in self.global_mechanisms:
            self.crn_species += \
                self.global_mechanisms[mech].update_species_global(
                                                            self.crn_species,
                                                            self.parameters)

        return self.crn_species

    def update_reactions(self) -> List[Reaction]:
        """ it generates the list of reactions based on all the mechanisms and global mechanisms
        it **must be** called after update_species() was called!

        :raise: AttributeError if it was called before update_species()
        :return: list of reactions generated by all the mechanisms and global mechanisms
        """
        if self.crn_species is None:
            raise AttributeError("Mixture.crn_species not defined. "
                                 "mixture.update_species() must be called "
                                 "before mixture.update_reactions()")

        self.crn_reactions = []
        for component in self.components:
            if self.parameter_warnings is not None:
                component.set_parameter_warnings(self.parameter_warnings)

            self.crn_reactions += component.update_reactions()

        # update with global mechanisms
        for mech in self.global_mechanisms:
            self.crn_reactions += \
                self.global_mechanisms[mech].update_reactions_global(
                                            self.crn_species, self.parameters)
        return self.crn_reactions

    def compile_crn(self) -> ChemicalReactionNetwork:
        """ Creates a chemical reaction network from the species and reactions associated with a mixture object
        :return: ChemicalReactionNetwork
        """
        resetwarnings()#Reset warnings - better to toggle them off manually.
        species = self.update_species()
        reactions = self.update_reactions()
        CRN = ChemicalReactionNetwork(species, reactions)
        return CRN

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

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
