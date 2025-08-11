#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import warnings
from typing import Dict, List, Tuple, Union
from warnings import warn
import numbers

import libsbml  # type: ignore

from .reaction import Reaction
from ..utils.sbmlutil import (
    add_all_reactions,
    add_all_species,
    add_all_compartments,
    create_sbml_model,
)
from .species import Species
from ..utils import (
    process_initial_concentration_dict,
    parameter_to_value,
    remove_bindloc,
)
from .parameter import ModelParameter, Parameter


class ChemicalReactionNetwork(object):
    r"""A chemical reaction network is a container of species and reactions
    chemical reaction networks can be compiled into SBML.

    reaction types:
    mass action: standard mass action semantics where the propensity of a
    reaction is given by deterministic propensity =
    .. math::
    k \Prod_{inputs i} [S_i]^a_i
    stochastic propensity =
    .. math::
    k \Prod_{inputs i} (S_i)!/(S_i - a_i)!
    where a_i is the spectrometric coefficient of species i
    """

    def __init__(
        self,
        species: List[Species],
        reactions: List[Reaction],
        initial_concentration_dict: Dict[
            Species, Union[numbers.Real, Parameter]
        ] = None,
        show_warnings=False,
    ):
        self.species = species
        self.reactions = reactions
        self.initial_concentration_dict = None  # Create an unpopulated dictionary
        self.initial_concentration_dict = initial_concentration_dict  # update it

        ChemicalReactionNetwork.check_crn_validity(
            self._reactions, self._species, show_warnings=show_warnings
        )

    @property
    def species(self):
        return copy.deepcopy(self._species)

    @species.setter
    def species(self, species):
        """Sets the species of the CRN object. If the species is not set,
        it initializes an empty list and adds the species to it.
        If the species is already set, it raises an AttributeError.
        A _species_set is used to ensure that no duplicate species are added to the CRN.
        In earlier BioCRNPyler versions, a _species_dict was used.
        This dictionary had key as species, and value as "True", which is unintuitive.
        Since, sets are designed to keep unique elements, so we use a set to keep track of species.

        Args:
            species (_type_): _description_

        Raises:
            AttributeError: _description_
        """
        if not hasattr(self, "_species"):
            self._species = []
            self._species_set = set()
            self.add_species(species)
        else:
            raise AttributeError(
                "The species in a CRN cannot be removed or modified. New Species can be added with CRN.add_species(...)."
            )

    @property
    def reactions(self):
        return copy.deepcopy(self._reactions)

    @reactions.setter
    def reactions(self, reactions):
        if not hasattr(self, "_reactions"):
            self._reactions = []
            self.add_reactions(reactions)
        else:
            raise AttributeError(
                "The reactions in a CRN cannot be removed or modified. New reactions can be added with CRN.add_reactions(...)."
            )

    def add_species(self, species, copy_species=True, compartment=None):
        """Adds a Species or a list of Species to the CRN object

        :param species: Species instance or list of Species instances
        :param copy_species: whether to deep copy Species added to the CRN. Protects CRN validity at teh expense of speed.

        """
        if not isinstance(species, list):
            species = [species]

        species = Species.flatten_list(species)  # Flatten the list
        species = remove_bindloc(species)

        # Deepcopy the specied list here, to preserve its structure
        if copy_species:
            species = copy.deepcopy(species)

        for s in species:
            if not isinstance(s, Species):  # check species are Species
                raise ValueError("A non-species object was used as a species!")
            if s not in self._species_set:  # Do not add duplicate Species
                if compartment is not None and s.compartment.name == "default":
                    s.compartment = compartment
                self._species_set.add(s)
                # copy the species and add it to the CRN
                self._species.append(s)

    def add_reactions(
        self,
        reactions: Union[Reaction, List[Reaction]],
        copy_reactions=True,
        add_species=True,
        compartment=None,
    ) -> None:
        """Adds a reaction or a list of reactions to the CRN object

        :param reactions: Reaction instance or list of Reaction instances
        :param copy_reactions: whether to deep copy reactions before adding
                                them to the CRN. Protects CRN validity at the
                                expense of speed.
        :param add_species: whether to add species in reactions to the CRN. 
                            Prevents errors at the expense of speed.
        :return: None
        """
        if not isinstance(reactions, list):
            reactions = [reactions]

        # It is recommended to copy reactions before adding them to the CRN, so they are "protected"
        if copy_reactions:
            reactions = copy.deepcopy(reactions)  # deep copy all the reactions

        # Add the reactions to the CRN
        self._reactions += reactions

        # Add species from reactions into the CRN
        if add_species:
            for r in reactions:
                if not isinstance(r, Reaction):  # check reactions and Reactions
                    raise ValueError("A non-reaction object was used as a reaction!")

                # add all the Species in the reaction to the CRN
                reaction_species = list(set([w.species for w in r.inputs + r.outputs]))
                self.add_species(
                    reaction_species,
                    copy_species=copy_reactions,
                    compartment=compartment,
                )

    @property
    def initial_concentration_dict(self):
        return self._initial_concentration_dict

    @initial_concentration_dict.setter
    def initial_concentration_dict(self, initial_concentration_dict):
        if initial_concentration_dict is None:
            self._initial_concentration_dict = {}
        elif isinstance(initial_concentration_dict, dict):
            for s in initial_concentration_dict:
                if s not in self._species_set:
                    raise ValueError(
                        f"Trying to set the initial concentration of a Species {s} not in the CRN"
                    )
                elif parameter_to_value(initial_concentration_dict[s]) >= 0:
                    self.initial_concentration_dict[s] = initial_concentration_dict[s]
                else:
                    raise ValueError(
                        f"Trying to set a species {s} to a negative concentration {initial_concentration_dict[s]}"
                    )

    @staticmethod
    def check_crn_validity(
        reactions: List[Reaction], species: List[Species], show_warnings=True
    ) -> Tuple[List[Reaction], List[Species]]:
        """Checks that the given list of reactions and list of species can form a valid CRN.

        :param reactions: list of reaction
        :param species: list of species
        :param show_warnings: whether to show warning when duplicated reactions/species was found
        :return: tuple(reaction,species)
        """

        if not all(isinstance(r, Reaction) for r in reactions):
            raise ValueError("A non-reaction object was used as a reaction!")

        if not all(isinstance(s, Species) for s in species):
            raise ValueError("A non-species object was used as a species!")

        for r in reactions:
            if reactions.count(r) > 1 and show_warnings:
                warn(
                    f"Reaction {r} may be duplicated in CRN definitions. "
                    f"Duplicates have NOT been removed."
                )

        for s in species:
            if species.count(s) > 1 and show_warnings:
                warn(
                    f"Species {s} is duplicated in the CRN definition. "
                    f"Duplicates have NOT been removed."
                )

        # check that all species in the reactions are also in the species list and vice versa
        unique_species = set(species)
        all_species_in_reactions = set(
            Species.flatten_list([r.species for r in reactions])
        )
        if unique_species != all_species_in_reactions:
            species_without_reactions = unique_species - all_species_in_reactions
            if species_without_reactions and show_warnings:
                warn(
                    f"These Species {list(species_without_reactions)} are not part of any reactions in the CRN!"
                )
            unlisted_reactions = all_species_in_reactions - unique_species
            if unlisted_reactions and show_warnings:
                warn(
                    f"These Species {list(unlisted_reactions)} are not listed in the Species list, but part of the reactions!"
                )

        return reactions, species

    def __repr__(self):
        txt = "Species = "
        for s in self._species:
            txt += repr(s) + ", "
        txt = txt[:-2] + "\n"
        txt += "Reactions = [\n"

        for r in self._reactions:
            txt += "\t" + repr(r) + "\n"
        txt += "]"
        return txt

    def pretty_print(
        self,
        show_rates=True,
        show_material=True,
        show_attributes=True,
        show_initial_concentration=True,
        show_keys=True,
        show_compartment=False,
        **kwargs,
    ):
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string identifiers.
        `show_material` toggles whether species.material is printed.
        `show_attributes` toggles whether species.attributes is printed
        `show_rates` toggles whether reaction rate functions are printed
        `show_compartment` toggles whether species.compartment is printed
        """

        txt = "Species" + f"(N = {len(self._species)}) = " + "{\n"

        def ics(s):
            return (
                self.initial_concentration_dict[s]
                if s in self.initial_concentration_dict
                else 0
            )

        species_sort_list = [(parameter_to_value(ics(s)), s) for s in self._species]
        species_sort_list.sort()
        species_sort_list.reverse()
        for sind, (init_conc, s) in enumerate(species_sort_list):
            init_conc = ics(s)

            txt += s.pretty_print(
                show_material=show_material,
                show_compartment=show_compartment,
                show_attributes=show_attributes,
                **kwargs,
            )

            if show_initial_concentration:
                txt += f" (@ {parameter_to_value(init_conc)}),  "

                if show_keys:  # shows where the initial conditions came from
                    if isinstance(init_conc, ModelParameter):
                        txt += f"\n   found_key=(mech={init_conc.found_key.mechanism}, partid={init_conc.found_key.part_id}, name={init_conc.found_key.name}).\n   search_key=(mech={init_conc.search_key.mechanism}, partid={init_conc.search_key.part_id}, name={init_conc.search_key.name}).\n"

        txt += "\n}\n"
        txt += f"\nReactions ({len(self._reactions)}) = [\n"

        for rind in range(len(self._reactions)):
            r = self._reactions[rind]
            txt += (
                f"{rind}. "
                + r.pretty_print(
                    show_rates=show_rates,
                    show_material=show_material,
                    show_attributes=show_attributes,
                    show_keys=show_keys,
                    **kwargs,
                )
                + "\n"
            )
        txt += "]"
        return txt

    def initial_condition_vector(
        self, init_cond_dict: Union[Dict[str, float], Dict[Species, float]]
    ):
        x0 = [0.0] * len(self._species)
        for idx, s in enumerate(self._species):
            if s in init_cond_dict:
                x0[idx] = init_cond_dict[s]
        return x0

    def get_all_species_containing(self, species: Species, return_as_strings=False):
        """Returns all species (complexes and otherwise) containing a given species
        (or string).
        """
        return_list = []
        if not isinstance(species, Species):
            raise ValueError("species argument must be an instance of Species!")

        for s in self._species:
            if species in s.get_species(recursive=True):
                if return_as_strings:
                    return_list.append(repr(s))
                else:
                    return_list.append(s)
        return return_list

    def replace_species(self, species: Species, new_species: Species):
        """Replaces species with new_species in the entire CRN.

        Does not act in place: returns a new CRN.
        """

        if not isinstance(species, Species):
            raise ValueError("species argument must be an instance of Species!")

        if not isinstance(new_species, Species):
            raise ValueError("species argument must be an instance of Species!")

        new_species_list = []
        for s in self._species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_reaction_list = []
        for r in self._reactions:
            new_r = r.replace_species(species, new_species)
            new_reaction_list.append(new_r)

        return ChemicalReactionNetwork(new_species_list, new_reaction_list)

    def generate_sbml_model(
        self,
        stochastic_model=False,
        show_warnings=False,
        check_validity=True,
        **keywords,
    ):
        """Creates an new SBML model and populates with the species and
        reactions in the ChemicalReactionNetwork object

        :param stochastic_model: whether the model is stochastic
        :param show_warnings: of from check crn validity
        :param keywords: extra keywords pass onto create_sbml_model() and add_all_reactions()
        :return: tuple: (document,model) SBML objects
        """
        if check_validity:
            ChemicalReactionNetwork.check_crn_validity(
                self._reactions, self._species, show_warnings=show_warnings
            )

        document, model = create_sbml_model(**keywords)
        all_compartments = []
        for species in self._species:
            if species.compartment not in all_compartments:
                all_compartments.append(species.compartment)
        add_all_compartments(model=model, compartments=all_compartments, **keywords)

        add_all_species(
            model=model,
            species=self._species,
            initial_condition_dictionary=self.initial_concentration_dict,
        )
        add_all_reactions(
            model=model,
            reactions=self._reactions,
            stochastic_model=stochastic_model,
            **keywords,
        )

        if document.getNumErrors():
            warn(
                "SBML model generated has errors. Use document.getErrorLog() to print all errors."
            )
        return document, model

    def write_sbml_file(
        self, file_name=None, stochastic_model=False, check_validity=True, **keywords
    ) -> bool:
        """ "Writes CRN object to a SBML file

        :param file_name: name of the file where the SBML model gets written
        :param stochastic_model: export an SBML file which ready for stochastic simulations
        :param keywords: keywords that passed into generate_sbml_model()
        :return: bool, show whether the writing process was successful
        """
        document, _ = self.generate_sbml_model(
            stochastic_model=stochastic_model, check_validity=check_validity, **keywords
        )
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, "w") as f:
            f.write(sbml_string)
        return True

    def simulate_with_bioscrape(
        self,
        timepoints,
        initial_condition_dict=None,
        stochastic=False,
        return_dataframe=True,
        safe=False,
    ):
        """Simulate CRN model with bioscrape (https://github.com/biocircuits/bioscrape).
        Returns the data for all species as Pandas dataframe.
        """
        result = None
        warnings.warn(
            "simulate_with_bioscrape is depricated and will cease working in a future release. Instead, please use simulate_with_bioscrape_via_sbml."
        )

        result = self.simulate_with_bioscrape_via_sbml(
            timepoints,
            filename=None,
            initial_condition_dict=initial_condition_dict,
            return_dataframe=return_dataframe,
            safe=safe,
            stochastic=stochastic,
        )

        return result

    def simulate_with_bioscrape_via_sbml(
        self,
        timepoints,
        filename=None,
        initial_condition_dict=None,
        return_dataframe=True,
        stochastic=False,
        safe=False,
        return_model=False,
        check_validity=True,
        **kwargs,
    ):
        """Simulate CRN model with bioscrape via writing a SBML file temporarily.
        [Bioscrape on GitHub](https://github.com/biocircuits/bioscrape).

        Returns the data for all species as Pandas dataframe.
        """
        result = None
        m = None
        try:
            from bioscrape.simulator import py_simulate_model  # type: ignore
            from bioscrape.types import Model  # type: ignore

            if filename is None:
                self.write_sbml_file(
                    file_name="temp_sbml_file.xml",
                    stochastic_model=stochastic,
                    for_bioscrape=True,
                    check_validity=check_validity,
                )
                file_name = "temp_sbml_file.xml"
            elif isinstance(filename, str):
                file_name = filename
            else:
                raise ValueError(
                    f"filename must be None or a string. Recievied: {filename}"
                )

            if "sbml_warnings" in kwargs:
                sbml_warnings = kwargs.get("sbml_warnings")
            else:
                sbml_warnings = False
            m = Model(sbml_filename=file_name, sbml_warnings=sbml_warnings)
            # m.write_bioscrape_xml('temp_bs'+ file_name + '.xml') # Uncomment if you want a bioscrape XML written as well.
            if initial_condition_dict is not None:
                processed = process_initial_concentration_dict(initial_condition_dict)
                m.set_species(processed)
            result = py_simulate_model(
                timepoints,
                Model=m,
                stochastic=stochastic,
                safe=safe,
                return_dataframe=return_dataframe,
            )
        except ModuleNotFoundError:
            warnings.warn("bioscrape was not found, please install bioscrape")

        if return_model:
            return result, m
        else:
            return result

    def simulate_with_roadrunner(
        self,
        timepoints: List[float],
        initial_condition_dict: Dict[str, float] = None,
        return_roadrunner=False,
        check_validity=True,
    ):
        """To simulate using roadrunner.
        Arguments:
        timepoints: The array of time points to run the simulation for.
        initial_condition_dict:

        Returns the results array as returned by RoadRunner OR a Roadrunner model object.

        Refer to the libRoadRunner simulator library documentation
        for details on simulation results: (http://libroadrunner.org/)[http://libroadrunner.org/]
        NOTE : Needs roadrunner package installed to simulate.
        """
        res_ar = None
        try:
            import roadrunner  # type: ignore
            import io

            document, _ = self.generate_sbml_model(
                stochastic_model=False, check_validity=check_validity
            )
            sbml_string = libsbml.writeSBMLToString(document)
            # write the sbml_string into a temporary file in memory instead of a file
            string_out = io.StringIO()
            string_out.write(sbml_string)
            # use the temporary file in memory to load the model into libroadrunner
            rr = roadrunner.RoadRunner(string_out.getvalue())
            if initial_condition_dict:
                for species, value in initial_condition_dict.items():
                    rr.model[f"init([{species}])"] = value

            if return_roadrunner:
                return rr
            else:
                result = rr.simulate(timepoints[0], timepoints[-1], len(timepoints))
                res_ar = result
        except ModuleNotFoundError:
            warnings.warn("libroadrunner was not found, please install libroadrunner")
        return res_ar
