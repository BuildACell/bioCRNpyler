
# Copyright (c) 2019, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from warnings import resetwarnings

from .component import Component
from .chemical_reaction_network import ChemicalReactionNetwork, Species, Reaction
from .parameter import Parameter
from typing import List, Union


class Mixture(object):
    def __init__(self, name="", mechanisms={}, components = [], parameters=None,
                 parameter_file = None, default_mechanisms = {},
                 global_mechanisms = {}, default_components = [],
                 species = [], custom_initial_condition = {},
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

        init = kwargs.get('init')
        parameter_warnings = kwargs.get('parameter_warnings')
        if parameter_warnings:
            warn('Parameter warnings have been set True. Verbose warnings regarding parameter files will be displayed.')
        else:
            parameter_warnings = False
            kwargs['parameter_warnings'] = parameter_warnings
        if not init and parameter_warnings:
            warn('Initial concentrations for extract species will all be set to zero.')


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

        #Initial conditions are searched for by defauled in the parameter file
        #see Mixture.set_initial_condition(self)
        #These can be overloaded with custom_initial_condition dictionary: component.name --> initial amount
        self.custom_initial_condition = custom_initial_condition

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


    

    #Used to set internal species froms strings, Species or Components
    def set_species(self, species, material_type = None, attributes = None):
        if isinstance(species, Species):
                return species
        elif isinstance(species, str):
            return Species(name = species, material_type = material_type, attributes = attributes)
        elif isinstance(species, Component) and species.get_species() is not None:
            return species.get_species()
        else:
            raise ValueError("Invalid Species: string, chemical_reaction_network.Species or Component with implemented .get_species() required as input.")


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

    #Sets the initial condition for all components with internal species
    #Does this for the species returned during compilation to prevent errors
    # First checks if (mixture.name, repr(species) is in the self.custom_initial_condition_dict
    # Then checks if (repr(species) is in the self.custom_initial_condition_dict
    # First checks if (mixture.name, component.name) is in the self.custom_initial_condition_dictionary
    # Then checks if (component.name) is in the self.custom_initial_condition_dictionary

    # First checks if (mixture.name, repr(species) is in the parameter dictionary
    # Then checks if repr(species) is in the parameter dictionary
    # Then checks if (mixture.name, component.name) is in the parameter dictionary
    # Then checks if component.name is in the parameter dictionary
    # Then defaults to 0
    def set_initial_condition(self, species):
        return_species = []
        for s in species:
            found = False
            if not found and (self.name, repr(s)) in self.custom_initial_condition:
                s.initial_concentration = self.custom_initial_condition[self.name, repr(s)]
                found = True
            elif not found and repr(s) in self.custom_initial_condition:
                s.initial_concentration = self.custom_initial_condition[repr(s)]
                found = True
            elif not found:
                for comp in self.components:

                    s_comp = comp.get_species()
                    if repr(s_comp) == repr(s):
                        if not found and (self.name, comp.name) in self.custom_initial_condition:
                            s.initial_concentration = self.custom_initial_condition[(self.name, comp.name)]
                            s_comp.initial_concentration = self.custom_initial_condition[(self.name, comp.name)]
                            found = True
                        elif not found and comp.name in self.custom_initial_condition:
                            s.initial_concentration = self.custom_initial_condition[comp.name]
                            s_comp.initial_concentration = self.custom_initial_condition[comp.name]
                            found = True

            if not found and (self.name, repr(s)) in self.parameters:
                s.initial_concentration = self.parameters[self.name, repr(s)]
                s_comp.initial_concentration = self.parameters[self.name, repr(s)]
                found = True
            elif not found and repr(s) in self.parameters:
                s.initial_concentration = self.parameters[repr(s)]
                s_comp.initial_concentration = self.parameters[repr(s)]
                found = True
            elif not found:
                for comp in self.components:
                    s_comp = comp.get_species()
                    if not found and repr(s_comp) == repr(s):
                        if not found and (self.name, comp.name) in self.parameters:
                            s.initial_concentration = self.parameters[(self.name, comp.name)]
                            s_comp.initial_concentration = self.parameters[(self.name, comp.name)]
                            found = True
                        elif not found and comp.name in self.parameters:
                            s.initial_concentration = self.parameters[comp.name]
                            s_comp.initial_concentration = self.parameters[comp.name]
                            found = True
        return species
    
    def append_species(self, new_species):
        self.crn_species += [s for s in new_species if s not in self.crn_species]

    def update_species(self) -> List[Species]:
        """ it generates the list of species based on all the mechanisms and global mechanisms

        :return: list of species generated by all the mechanisms and global mechanisms
        """
        # TODO check if we can merge the two variables
        self.crn_species = self.added_species
        for component in self.components:
            self.append_species(component.update_species())

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
                self.global_mechanisms[mech].update_reactions_global(self.crn_species, self.parameters)
        return self.crn_reactions

    def compile_crn(self) -> ChemicalReactionNetwork:
        """ Creates a chemical reaction network from the species and reactions associated with a mixture object
        :return: ChemicalReactionNetwork
        """
        resetwarnings()#Reset warnings - better to toggle them off manually.

        species = self.update_species()
        reactions = self.update_reactions()
        species = self.set_initial_condition(species)
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
            txt+="\n\t"+mech+":"+self.global_mechanisms[mech].name
        txt+=" }"
        return txt




