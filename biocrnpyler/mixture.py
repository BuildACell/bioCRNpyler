
# Copyright (c) 2019, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from warnings import warn
from warnings import resetwarnings

from .component import Component
from .chemical_reaction_network import ChemicalReactionNetwork, Species, Reaction
from .parameter import ParameterDatabase, ParameterEntry
from typing import List, Union


class Mixture(object):
    def __init__(self, name="", mechanisms=None, components=None, parameters=None, parameter_file=None,
                 default_mechanisms=None, global_mechanisms=None, species=None, initial_condition_dictionary=None,
                 parameter_warnings=None, overwrite_parameters = False,**kwargs):
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
        :param parameter_warnings: suppressing parameter related warnings
        """
        if components is None:
            components = []

        init = kwargs.get('init')
        parameter_warnings = kwargs.get('parameter_warnings')
        if parameter_warnings:
            warn('Parameter warnings have been set True. Verbose warnings regarding parameter files will be displayed.')
        else:
            parameter_warnings = False
            kwargs['parameter_warnings'] = parameter_warnings


        # Initialize instance variables
        self.name = name  # Save the name of the mixture

        self.parameter_database = ParameterDatabase(parameter_file = parameter_file, parameter_dictionary = parameters, overwrite_parameters = overwrite_parameters)
        
        # Toggles whether parameter warnings are raised. if None (default) this
        # can be toggled component by component.
        self.parameter_warnings = parameter_warnings

        # Override the default mechanisms with anything we were passed
        # default parameters are used by mixture subclasses.
        self.default_mechanisms = default_mechanisms
        self.custom_mechanisms = mechanisms

        # Initial conditions are searched for by defauled in the parameter file
        # see Mixture.set_initial_condition(self)
        # These can be overloaded with custom_initial_condition dictionary: component.name --> initial amount
        if initial_condition_dictionary is None:
            self.initial_condition_dictionary = {}
        else:
            self.initial_condition_dictionary = dict(initial_condition_dictionary)

        # Mechanisms stores the mechanisms used for compilation where defaults
        # are overwritten by custom mechanisms.
        self.mechanisms = self.default_mechanisms

        if self.custom_mechanisms:
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
        # process the components
        self.add_components(components)

        # internal lists for the species and reactions
        self.crn_species = None
        self.crn_reactions = None

    def add_species(self, species: Union[List[Species], Species]):
        if species is not None:
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


    def add_components(self, components: Union[List[Component],Component]):
        if not isinstance(components, list):
            components = [components]

        for component in components:
            assert isinstance(component, Component), \
                "the object: %s passed into mixture as component must be of the class Component" % str(component)
            self.components.append(component)

            #Reset components Mixtures
            component.set_mixture(self)

            component.update_mechanisms(mixture_mechanisms=self.mechanisms, overwrite_custom_mechanisms = False)
            #component.update_parameters(mixture_parameters=self.parameters)
            if self.parameter_warnings is not None:
                component.set_parameter_warnings(self.parameter_warnings)

    def update_parameters(self, parameter_file = None, parameters = None, overwrite_parameters = True):
        if parameter_file is not None:
            self.load_parameters_from_file.load_parameters_from_dictionary(parameter_file, overwrite_parameters = overwrite_parameters)

        if parameters is not None:
            self.parameter_database.load_parameters_from_dictionary(parameter_dictionary, overwrite_parameters = overwrite_parameters)
    
    def get_parameter(self, mechanism, part_id, param_name, parameter_warnings = False):
        param = self.parameter_database.find_parameter(mechanism, part_id, param_name, parameter_warnings = parameter_warnings)

        #TODO replace this with just returning the parameter, when reaction overhaul is complete
        if isinstance(param, ParameterEntry):
            return param.value
        else:
            return param

    #Sets the initial condition for all components with internal species
    #Does this for the species returned during compilation to prevent errors
    # First checks if (mixture.name, repr(species) is in the self.custom_initial_condition_dict
    # Then checks if (repr(species) is in the self.custom_initial_condition_dict
    # First checks if (mixture.name, component.name) is in the self.custom_initial_condition_dictionary
    # Then checks if (component.name) is in the self.custom_initial_condition_dictionary

    # First checks if (mixture.name, repr(species) is in the ParameterDatabase
    # Then checks if repr(species) is in the parameter dictionary
    # Then checks if (mixture.name, component.name) is in the ParameterDatabase
    # Then checks if component.name is in the parameter dictionary
    # Then defaults to 0
    def set_initial_condition_old(self, species):
        return_species = []
        for s in species:
            found = False
            if not found and (self.name, repr(s)) in self.initial_condition_dictionary:
                s.initial_concentration = self.custom_initial_condition[self.name, repr(s)]
                found = True
            elif not found and (self.name, s) in self.initial_condition_dictionary:
                s.initial_concentration = self.custom_initial_condition[self.name, s]
                found = True
            elif not found and repr(s) in self.initial_condition_dictionary:
                s.initial_concentration = self.initial_condition_dictionary[repr(s)]
                found = True
            elif not found and s in self.initial_condition_dictionary:
                s.initial_concentration = self.initial_condition_dictionary[s]
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
    
    #Tries to find an initial condition of species s using the parameter heirarchy
    # 1. Tries to find the initial concentration in the Component initial_Concentration_dictionary and ParameterDatabase
    # 2. Tries to find self.name, repr(s) in self.initial_condition_dictionary
    # 3. Tries to find repr(s) in self.initial_condition_dictionary
    # 4. tries to find self.name, repr(s) in self.parameter_database
    # 5. tries to find repr(s) in self.parameter_database
    # 6. defaults to 0
    def set_initial_condition(self, s, component = None):
        if not isinstance(s, Species):
            raise ValueError(f"{s} is not a Species! Can only set initial concentration of a Species.")
        init_conc = None
        if component is not None:
            init_conc = component.find_initial_condition(s)

        if init_conc is None:
            if (self.name, repr(s)) in self.initial_condition_dictionary:
                init_conc = self.initial_condition_dictionary[(self.name, repr(s))]
            elif repr(s) in self.initial_condition_dictionary:
                init_conc = self.initial_condition_dictionary[repr(s)]
            elif self.parameter_database.find_parameter(None, self.name, repr(s)) is not None:
                init_conc = self.parameter_database.find_parameter(None, self.name, repr(s))
            elif self.parameter_database.find_parameter(None, None, repr(s)) is not None:
                init_conc = self.parameter_database.find_parameter(None, None, repr(s))
            else:
                init_conc = 0

        s.initial_concentration = init_conc


    #Allows mechanisms to return nested lists of species. These lists are flattened.
    def append_species(self, new_species, component):
        for s in new_species:
            if isinstance(s, Species):
                self.set_initial_condition(s, component)
                self.crn_species.append(s)
            elif isinstance(s, list) and(all(isinstance(ss, Species) for ss in s) or len(s) == 0):
                for ss in s: 
                    self.set_initial_condition(ss, component)
                self.crn_species+=s
            elif s is not None:
                raise ValueError(f"Invalid Species Returned in {component}.update_species(): {s}.")
        #Old Version
        #self.crn_species += [s for s in new_species if s not in self.crn_species]


    def update_species(self) -> List[Species]:
        """ it generates the list of species based on all the mechanisms in each Component
        :return: list of species generated by all the mechanisms and global mechanisms
        """
        self.crn_species = []
        #Append Species added manually
        self.append_species(self.added_species, None)

        #Appendy species from each Component
        for component in self.components:
            self.append_species(component.update_species(), component)

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

        return self.crn_reactions

        

    def apply_global_mechanisms(self) -> (List[Species], List[Reaction]):
        # update with global mechanisms

        if self.crn_species is None:
            raise AttributeError("Mixture.crn_species not defined. "
                                 "mixture.update_species() must be called "
                                 "before mixture.apply_global_mechanisms()")

        global_mech_species = []
        global_mech_reactions = []
        if self.global_mechanisms:
            for mech in self.global_mechanisms:
                # Update Global Mechanisms
                global_mech_species += self.global_mechanisms[mech].update_species_global(self.crn_species, self)
                global_mech_reactions += self.global_mechanisms[mech].update_reactions_global(self.crn_species, self)

        return global_mech_species, global_mech_reactions


    def compile_crn(self) -> ChemicalReactionNetwork:
        """ Creates a chemical reaction network from the species and reactions associated with a mixture object
        :return: ChemicalReactionNetwork
        """
        resetwarnings()#Reset warnings - better to toggle them off manually.

        #reset the Components' mixture to self - in case they have been added to other Mixtures
        for c in self.components:
            c.set_mixture(self)

        self.update_species() #updates species to self.crn_species and sets initial concentrations
        self.update_reactions() #updates reactions to self.crn_reactions

        #global mechanisms are applied last and only to all the species 
        global_mech_species, global_mech_reactions = self.apply_global_mechanisms()

        #append global species to self.crn_species and update initial concentrations
        if isinstance(global_mech_species, list) and len(global_mech_species) > 0:
            self.append_species(global_mech_species, component = None)

        #append global reactions
        if isinstance(global_mech_reactions, list) and len(global_mech_reactions)>0:
            self.crn_reactions += global_mech_reactions 

        CRN = ChemicalReactionNetwork(list(self.crn_species), list(self.crn_reactions))
        return CRN

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self)+"\n"
        if self.components:
            txt += "Components = ["
            for comp in self.components:
                txt+="\n\t"+str(comp)
        if self.mechanisms:
            txt+=" ]\nMechanisms = {"
            for mech in self.mechanisms:
                txt+="\n\t"+mech+":"+self.mechanisms[mech].name
        if self.global_mechanisms:
            txt+=" }\nGlobal Mechanisms = {"
            for mech in self.global_mechanisms:
                txt+="\n\t"+mech+":"+self.global_mechanisms[mech].name
        txt+=" }"
        return txt




