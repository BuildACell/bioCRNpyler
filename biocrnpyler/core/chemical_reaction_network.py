#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import numbers
import warnings
from typing import Dict, List, Tuple, Union
from warnings import warn

import libsbml  # type: ignore

from ..utils import (
    parameter_to_value,
    process_initial_concentration_dict,
    remove_bindloc,
)
from ..utils.sbmlutil import (
    add_all_compartments,
    add_all_reactions,
    add_all_species,
    create_sbml_model,
)
from .parameter import ModelParameter, Parameter
from .reaction import Reaction
from .species import Species


class ChemicalReactionNetwork(object):
    r"""Container for chemical species and their reactions.

    A ChemicalReactionNetwork (CRN) represents a biochemical system as a set
    of species and reactions between them. CRNs can be compiled to SBML format
    for simulation with various tools, or simulated directly with bioscrape or
    roadrunner.

    Parameters
    ----------
    species : list of Species
        List of chemical species in the network.
    reactions : list of Reaction
        List of reactions between species. Each reaction specifies inputs,
        outputs, and rate parameters.
    initial_concentration_dict : dict, optional
        Dictionary mapping Species to their initial concentrations. Values can
        be numbers or `Parameter` objects. If None, an empty dictionary is
        created.
    show_warnings : bool, default=False
        If True, shows warnings about duplicate species/reactions or
        inconsistencies during CRN validation.

    Attributes
    ----------
    species : list of Species
        Deep copy of the species list (use `add_species` to modify).
    reactions : list of Reaction
        Deep copy of the reactions list (use `add_reactions` to modify).
    initial_concentration_dict : dict
        Dictionary of initial concentrations for species in the CRN.

    See Also
    --------
    Species : Chemical species in a CRN.
    Reaction : Chemical reaction between species.
    Mixture : High-level interface for building CRNs from components.

    Notes
    -----
    Mass action reactions follow standard mass action kinetics:

    - Deterministic propensity: $k \prod_i [S$_i$]^{a_i}$
    - Stochastic propensity: $k \prod_i \frac{S_i!}{(S_i - a_i)!}$

    where $a_i$ is the stoichiometric coefficient of species $i$.

    A valid CRN requires:

    - All species in reactions must be in the species list
    - All species and reactions must be unique (duplicates trigger warnings)
    - Initial concentrations must be non-negative

    Once created, species and reactions cannot be removed, only added. This
    ensures CRN validity is maintained throughout its lifetime.

    Chemical reaction networks can be simulated by writing the output
    as SMBL using `write_sbml_file` and then loading into an external
    simulator, or by using the bioscrape package, which can be called
    directly using `simulate_with_bioscrape_via_sbml`.

    Examples
    --------
    Create a simple CRN manually:

    >>> # Define species
    >>> S = bcp.Species('S')
    >>> P = bcp.Species('P')
    >>> E = bcp.Species('protein_E')
    >>> C = bcp.Species('C')
    >>> # Define reactions
    >>> rxn1 = bcp.Reaction.from_massaction(
    ...     [S, E], [C], k_forward=0.1, k_reverse=1e-4)
    >>> rxn2 = bcp.Reaction.from_massaction([C], [E, P], k_forward=0.01)
    >>> # Create CRN
    >>> crn = bcp.ChemicalReactionNetwork(
    ...     species=[S, E, C, P],
    ...     reactions=[rxn1, rxn2]
    ... )

    Compile a CRN from a mixture:

    >>> enzyme = bcp.Enzyme('E', 'S', 'P')
    >>> mixture = bcp.Mixture(
    ...     components=[enzyme],
    ...     mechanisms={'catalysis': bcp.MichaelisMenten()},
    ...     parameters={'kb': 0.1, 'ku': 1e-4, 'kcat': 0.01})
    >>> crn = mixture.compile_crn()
    >>> print(
    ...     f"CRN has {len(crn.species)} species and "
    ...     f"{len(crn.reactions)} reactions")
    CRN has 4 species and 2 reactions

    Export to SBML and simulate:

    >>> crn.write_sbml_file('model.xml')
    >>> result = crn.simulate_with_bioscrape_via_sbml(
    ...     initial_condition_dict={S: 100, 'protein_E': 50, P: 0},
    ...     timepoints=np.linspace(0, 5))

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
        self.initial_concentration_dict = (
            None  # Create an unpopulated dictionary
        )
        self.initial_concentration_dict = (
            initial_concentration_dict  # update it
        )

        ChemicalReactionNetwork.check_crn_validity(
            self._reactions, self._species, show_warnings=show_warnings
        )

    @property
    def species(self):
        """list: List of Species : Deep copy of all species in the CRN."""
        return copy.deepcopy(self._species)

    @species.setter
    def species(self, species):
        """Set the initial species list for the CRN.

        Parameters
        ----------
        species : list of Species
            Initial species to add to the CRN. Additional species can be added
            later using `add_species`.

        Raises
        ------
        AttributeError
            If species list is already set. Species cannot be removed once
            added, only new species can be added with `add_species`.

        Notes
        -----
        This setter can only be called once during CRN initialization. A
        `_species_set` is maintained internally to ensure no duplicate species
        are added to the CRN.

        """
        if not hasattr(self, '_species'):
            self._species = []
            self._species_set = set()
            self.add_species(species)
        else:
            raise AttributeError(
                "The species in a CRN cannot be removed or modified. "
                "New Species can be added with CRN.add_species(...)."
            )

    @property
    def reactions(self):
        """List of Reaction: Deep copy of all reactions in the CRN."""
        return copy.deepcopy(self._reactions)

    @reactions.setter
    def reactions(self, reactions):
        """Set the initial reactions list for the CRN.

        Parameters
        ----------
        reactions : list of Reaction
            Initial reactions to add to the CRN. Additional reactions can be
            added later using `add_reactions`.

        Raises
        ------
        AttributeError
            If reactions list is already set. Reactions cannot be removed once
            added, only new reactions can be added with `add_reactions`.

        """
        if not hasattr(self, '_reactions'):
            self._reactions = []
            self.add_reactions(reactions)
        else:
            raise AttributeError(
                "The reactions in a CRN cannot be removed or modified. "
                "New reactions can be added with CRN.add_reactions(...)."
            )

    def add_species(self, species, copy_species=True, compartment=None):
        """Add species to the CRN.

        Parameters
        ----------
        species : Species or list of Species
            Species object(s) to add to the CRN. Lists are automatically
            flattened and binding locations are removed.
        copy_species : bool, default=True
            If True, deep-copies species before adding them to the CRN.
            Protects CRN validity at the expense of speed.
        compartment : Compartment, optional
            If provided, assigns this compartment to any species with default
            compartments.

        Raises
        ------
        ValueError
            If any element is not a Species object.

        Notes
        -----
        Duplicate species (based on equality) are automatically filtered out.
        Species are stored in both a list (`_species`) and a set
        (`_species_set`) for efficient duplicate checking.

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
                raise ValueError(
                    "A non-species object was used as a species!"
                )
            if s not in self._species_set:  # Do not add duplicate Species
                if (
                    compartment is not None
                    and s.compartment.name == 'default'
                ):
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
        """Add reactions to the CRN.

        Parameters
        ----------
        reactions : Reaction or list of Reaction
            Reaction object(s) to add to the CRN.
        copy_reactions : bool, default=True
            If True, deep-copies reactions before adding them to the CRN.
            Protects CRN validity at the expense of speed.
        add_species : bool, default=True
            If True, automatically adds any species appearing in the reactions
            to the CRN. Prevents missing species errors at the expense of
            speed.
        compartment : Compartment, optional
            If provided, assigns this compartment to any species with default
            compartments found in the reactions.

        Raises
        ------
        ValueError
            If any element is not a Reaction object.

        Notes
        -----
        Unlike species, reactions are not checked for duplicates when added.
        It is recommended to keep `copy_reactions=True` to protect the CRN
        from external modifications.

        """
        if not isinstance(reactions, list):
            reactions = [reactions]

        # It is recommended to copy reactions before adding them to
        # the CRN, so they are "protected"
        if copy_reactions:
            reactions = copy.deepcopy(
                reactions
            )  # deep copy all the reactions

        # Add the reactions to the CRN
        self._reactions += reactions

        # Add species from reactions into the CRN
        if add_species:
            for r in reactions:
                if not isinstance(
                    r, Reaction
                ):  # check reactions and Reactions
                    raise ValueError(
                        "A non-reaction object was used as a reaction!"
                    )

                # add all the Species in the reaction to the CRN
                reaction_species = list(
                    set([w.species for w in r.inputs + r.outputs])
                )
                self.add_species(
                    reaction_species,
                    copy_species=copy_reactions,
                    compartment=compartment,
                )

    @property
    def initial_concentration_dict(self):
        """dict: Dictionary mapping Species to initial concentrations."""
        return self._initial_concentration_dict

    @initial_concentration_dict.setter
    def initial_concentration_dict(self, initial_concentration_dict):
        """Set initial concentrations for species in the CRN.

        Parameters
        ----------
        initial_concentration_dict : dict or None
            Dictionary mapping Species objects to their initial
            concentrations.  Values can be numbers or `Parameter` objects. If
            None, an empty dictionary is created.

        Raises
        ------
        ValueError
            If a species in the dictionary is not in the CRN, or if any
            concentration is negative.

        """
        if initial_concentration_dict is None:
            self._initial_concentration_dict = {}
        elif isinstance(initial_concentration_dict, dict):
            for s in initial_concentration_dict:
                if s not in self._species_set:
                    raise ValueError(
                        "Trying to set the initial concentration of a "
                        f"Species {s} not in the CRN"
                    )
                elif parameter_to_value(initial_concentration_dict[s]) >= 0:
                    self.initial_concentration_dict[s] = (
                        initial_concentration_dict[s]
                    )
                else:
                    raise ValueError(
                        f"Trying to set a species {s} to a negative "
                        f"concentration {initial_concentration_dict[s]}"
                    )

    @staticmethod
    def check_crn_validity(
        reactions: List[Reaction], species: List[Species], show_warnings=True
    ) -> Tuple[List[Reaction], List[Species]]:
        """Validate that reactions and species can form a valid CRN.

        Checks for duplicate species/reactions and verifies that all species
        in reactions are present in the species list.

        Parameters
        ----------
        reactions : list of Reaction
            List of reactions to validate.
        species : list of Species
            List of species to validate.
        show_warnings : bool, default=True
            If True, issues warnings for duplicates or inconsistencies.

        Returns
        -------
        tuple of (list of Reaction, list of Species)
            The input reactions and species lists, unchanged.

        Raises
        ------
        ValueError
            If any reaction is not a Reaction object, or any species is not a
            Species object.

        Warns
        -----
        UserWarning
            - Duplicate reactions or species are found
            - Species exist without reactions
            - Reactions contain unlisted species

        """
        if not all(isinstance(r, Reaction) for r in reactions):
            raise ValueError("A non-reaction object was used as a reaction!")

        if not all(isinstance(s, Species) for s in species):
            raise ValueError("A non-species object was used as a species!")

        for r in reactions:
            if reactions.count(r) > 1 and show_warnings:
                warn(
                    f"Reaction {r} may be duplicated in CRN definitions. "
                    "Duplicates have NOT been removed."
                )

        for s in species:
            if species.count(s) > 1 and show_warnings:
                warn(
                    f"Species {s} is duplicated in the CRN definition. "
                    "Duplicates have NOT been removed."
                )

        # check that all species in the reactions are also in the
        # species list and vice versa
        unique_species = set(species)
        all_species_in_reactions = set(
            Species.flatten_list([r.species for r in reactions])
        )
        if unique_species != all_species_in_reactions:
            species_without_reactions = (
                unique_species - all_species_in_reactions
            )
            if species_without_reactions and show_warnings:
                warn(
                    f"These Species {list(species_without_reactions)} are "
                    "not part of any reactions in the CRN!"
                )
            unlisted_reactions = all_species_in_reactions - unique_species
            if unlisted_reactions and show_warnings:
                warn(
                    f"These Species {list(unlisted_reactions)} are not "
                    "listed in the Species list, but part of the reactions!"
                )

        return reactions, species

    def __repr__(self):
        txt = 'Species = '
        for s in self._species:
            txt += repr(s) + ', '
        txt = txt[:-2] + '\n'
        txt += 'Reactions = [\n'

        for r in self._reactions:
            txt += '\t' + repr(r) + '\n'
        txt += ']'
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
        """Generate detailed, human-readable string representation of the CRN.

        Parameters
        ----------
        show_rates : bool, default=True
            If True, displays reaction rate functions and parameters.
        show_material : bool, default=True
            If True, displays species material types (e.g., 'dna', 'protein').
        show_attributes : bool, default=True
            If True, displays species attributes.
        show_initial_concentration : bool, default=True
            If True, displays initial concentrations for each species.
        show_keys : bool, default=True
            If True, shows parameter database keys for initial concentrations
            (useful for debugging parameter lookup).
        show_compartment : bool, default=False
            If True, displays compartment information for each species.
        **kwargs
            Additional keyword arguments passed to species and reaction
            `pretty_print` methods.

        Returns
        -------
        str
            Formatted string with species (sorted by initial concentration)
            and reactions with detailed information.

        Notes
        -----
        This method provides much more detailed output than `__repr__`,
        making it useful for debugging and understanding CRN structure.
        Species are sorted by initial concentration (highest first) for easier
        analysis.

        """
        txt = 'Species' + f"(N = {len(self._species)}) = " + '{\n'

        def ics(s):
            return (
                self.initial_concentration_dict[s]
                if s in self.initial_concentration_dict
                else 0
            )

        species_sort_list = [
            (parameter_to_value(ics(s)), s) for s in self._species
        ]
        species_sort_list.sort()
        species_sort_list.reverse()
        for sind, (init_conc, s) in enumerate(species_sort_list):
            init_conc = ics(s)

            txt += '    ' + s.pretty_print(
                show_material=show_material,
                show_compartment=show_compartment,
                show_attributes=show_attributes,
                **kwargs,
            )

            if show_initial_concentration:
                txt += f" (@ {parameter_to_value(init_conc)}),  "

                if show_keys:  # shows where the initial conditions came from
                    if isinstance(init_conc, ModelParameter):
                        txt += (
                            "\n    found_key=("
                            f"mech={init_conc.found_key.mechanism}, "
                            f"partid={init_conc.found_key.part_id}, "
                            f"name={init_conc.found_key.name})."
                            "\n    search_key=("
                            f"mech={init_conc.search_key.mechanism}, "
                            f"partid={init_conc.search_key.part_id}, "
                            f"name={init_conc.search_key.name}).\n"
                        )
            txt += '\n'

        txt += '}\n'
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
                + '\n'
            )
        txt += ']'
        return txt

    def initial_condition_vector(
        self, init_cond_dict: Union[Dict[str, float], Dict[Species, float]]
    ):
        """Generate an initial condition vector for simulations.

        Parameters
        ----------
        init_cond_dict : dict
            Dictionary mapping species (or species names as strings) to their
            initial concentrations.

        Returns
        -------
        list of float
            Vector of initial concentrations matching the order of species in
            `self._species`. Species not in `init_cond_dict` are set to 0.0.

        """
        x0 = [0.0] * len(self._species)
        for idx, s in enumerate(self._species):
            if s in init_cond_dict:
                x0[idx] = init_cond_dict[s]
        return x0

    def get_all_species_containing(
        self, species: Species, return_as_strings=False
    ):
        """Find all species (complexes) that contain a given species.

        Searches recursively through all species in the CRN to find those that
        contain the target species as a component.

        Parameters
        ----------
        species : Species
            The species to search for within other species.
        return_as_strings : bool, default=False
            If True, returns species as string representations. If False,
            returns actual Species objects.

        Returns
        -------
        list
            List of Species objects (or strings if `return_as_strings=True`)
            that contain the target species.

        Raises
        ------
        ValueError
            If `species` is not a Species object.

        Examples
        --------
        >>> substrate = bcp.Species('S')
        >>> enzyme = bcp.Enzyme('E', substrate, 'P')
        >>> mixture = bcp.Mixture(
        ...     components=[enzyme],
        ...     mechanisms={'catalysis': bcp.MichaelisMenten()},
        ...     parameters={'kb': 0.1, 'ku': 1e-4, 'kcat': 0.01})
        >>> crn = mixture.compile_crn()
        >>> # Find all complexes containing S
        >>> crn.get_all_species_containing(substrate)
        [S, complex_S_protein_E_]

        """
        return_list = []
        if not isinstance(species, Species):
            raise ValueError(
                "species argument must be an instance of Species!"
            )

        for s in self._species:
            if species in s.get_species(recursive=True):
                if return_as_strings:
                    return_list.append(repr(s))
                else:
                    return_list.append(s)
        return return_list

    def replace_species(self, species: Species, new_species: Species):
        """Replace a species with another throughout the CRN.

        Creates a new CRN where all occurrences of a target species are
        replaced with a new species. The original CRN is not modified.

        Parameters
        ----------
        species : Species
            The species to be replaced.
        new_species : Species
            The species to replace with.

        Returns
        -------
        ChemicalReactionNetwork
            New CRN with the species replacement applied.

        Raises
        ------
        ValueError
            If either argument is not a Species object.

        Notes
        -----
        This method does not modify the original CRN. It creates and returns
        a new CRN with the replacement applied throughout all species and
        reactions.

        """
        if not isinstance(species, Species):
            raise ValueError(
                "species argument must be an instance of Species!"
            )

        if not isinstance(new_species, Species):
            raise ValueError(
                "species argument must be an instance of Species!"
            )

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
        **kwargs,
    ):
        """Generate an SBML model from the CRN.

        Creates SBML document and model objects containing all species,
        reactions, compartments, and parameters from the CRN.

        Parameters
        ----------
        stochastic_model : bool, default=False
            If True, generates an SBML model configured for stochastic
            simulation.
        show_warnings : bool, default=False
            If True, shows warnings from CRN validity checking.
        check_validity : bool, default=True
            If True, validates the CRN before generating SBML.
        **kwargs
            Additional keyword arguments passed to `create_sbml_model` and
            `add_all_reactions`.

        Returns
        -------
        tuple of (libsbml.SBMLDocument, libsbml.Model)
            The SBML document and model objects. The document can be written
            to a file or further manipulated.

        Warns
        -----
        UserWarning
            Issues a warning if the generated SBML model contains errors.

        See Also
        --------
        write_sbml_file : Write the SBML model directly to a file.

        """
        if check_validity:
            ChemicalReactionNetwork.check_crn_validity(
                self._reactions, self._species, show_warnings=show_warnings
            )

        document, model = create_sbml_model(**kwargs)
        all_compartments = []
        for species in self._species:
            if species.compartment not in all_compartments:
                all_compartments.append(species.compartment)
        add_all_compartments(
            model=model, compartments=all_compartments, **kwargs
        )

        add_all_species(
            model=model,
            species=self._species,
            initial_condition_dictionary=self.initial_concentration_dict,
        )
        add_all_reactions(
            model=model,
            reactions=self._reactions,
            stochastic_model=stochastic_model,
            **kwargs,
        )

        if document.getNumErrors():
            warn(
                "SBML model generated has errors. Use document.getErrorLog() "
                "to print all errors."
            )
        return document, model

    def write_sbml_file(
        self,
        file_name=None,
        stochastic_model=False,
        check_validity=True,
        **kwargs,
    ) -> bool:
        """Write the CRN to an SBML file.

        Generates an SBML model from the CRN and writes it to a file for use
        with simulators like COPASI, VCell, or bioscrape.

        Parameters
        ----------
        file_name : str
            Path where the SBML file will be written.
        stochastic_model : bool, default=False
            If True, exports an SBML file configured for stochastic
            simulations.
        check_validity : bool, default=True
            If True, validates the CRN before generating SBML.
        **kwargs
            Additional keyword arguments passed to `generate_sbml_model`.

        Returns
        -------
        bool
            True if the file was written successfully.

        See Also
        --------
        generate_sbml_model : Generate SBML objects without writing to file.

        Examples
        --------
        >>> crn.write_sbml_file('my_model.xml')
        >>> crn.write_sbml_file('stochastic_model.xml', stochastic_model=True)

        """
        document, _ = self.generate_sbml_model(
            stochastic_model=stochastic_model,
            check_validity=check_validity,
            **kwargs,
        )
        sbml_string = libsbml.writeSBMLToString(document)
        with open(file_name, 'w') as f:
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
        """Simulate CRN with bioscrape.

        .. deprecated:: 1.0.0
            This method is deprecated. Use
            `simulate_with_bioscrape_via_sbml` instead.

        Parameters
        ----------
        timepoints : array-like
            Time points for simulation.
        initial_condition_dict : dict, optional
            Dictionary of initial concentrations.
        stochastic : bool, default=False
            If True, runs stochastic simulation.
        return_dataframe : bool, default=True
            If True, returns results as pandas DataFrame.
        safe : bool, default=False
            Safe mode for bioscrape simulation.

        Returns
        -------
        DataFrame or array
            Simulation results.

        See Also
        --------
        simulate_with_bioscrape_via_sbml : Recommended simulation method.

        """
        result = None
        warnings.warn(
            "simulate_with_bioscrape is depricated and will cease working in "
            "a future release. Instead, please use "
            "simulate_with_bioscrape_via_sbml."
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
        """Simulate CRN with bioscrape via SBML export.

        Exports the CRN to an SBML file and simulates it using the bioscrape
        simulator. Bioscrape is a stochastic and deterministic simulator for
        biological circuits.

        Parameters
        ----------
        timepoints : array-like
            Array of time points at which to record simulation results.
        filename : str, optional
            Path to save the SBML file. If None, creates a temporary file
            'temp_sbml_file.xml'.
        initial_condition_dict : dict, optional
            Dictionary mapping species to initial concentrations. Overrides
            the CRN's `initial_concentration_dict`.
        return_dataframe : bool, default=True
            If True, returns results as a pandas DataFrame. If False, returns
            a numpy array.
        stochastic : bool, default=False
            If True, runs stochastic (Gillespie) simulation. If False, runs
            deterministic (ODE) simulation.
        safe : bool, default=False
            If True, uses bioscrape's safe mode which checks for errors.
        return_model : bool, default=False
            If True, returns a tuple of (results, bioscrape_model). If False,
            returns only results.
        check_validity : bool, default=True
            If True, validates the CRN before generating SBML.
        **kwargs
            Additional keyword arguments. 'sbml_warnings' can be set to True
            to show SBML parsing warnings.

        Returns
        -------
        DataFrame or array, or tuple
            Simulation results as DataFrame or array. If `return_model=True`,
            returns tuple of (results, bioscrape Model object).

        Warns
        -----
        UserWarning
            Issues a warning if bioscrape is not installed.

        Notes
        -----
        Requires bioscrape to be installed: `pip install bioscrape`

        Bioscrape GitHub: https://github.com/biocircuits/bioscrape

        Examples
        --------
        >>> result = crn.simulate_with_bioscrape_via_sbml(
        ...     timepoints=np.linspace(0, 10, 100)
        ... )
        >>> # Stochastic simulation
        >>> result = crn.simulate_with_bioscrape_via_sbml(
        ...     timepoints=np.linspace(0, 10, 100),
        ...     stochastic=True
        ... )

        """
        result = None
        m = None
        try:
            from bioscrape.simulator import py_simulate_model  # type: ignore
            from bioscrape.types import Model  # type: ignore

            if filename is None:
                self.write_sbml_file(
                    file_name='temp_sbml_file.xml',
                    stochastic_model=stochastic,
                    for_bioscrape=True,
                    check_validity=check_validity,
                )
                file_name = 'temp_sbml_file.xml'
            elif isinstance(filename, str):
                file_name = filename
            else:
                raise ValueError(
                    "filename must be None or a string. Recievied: "
                    f"{filename}"
                )

            if 'sbml_warnings' in kwargs:
                sbml_warnings = kwargs.get('sbml_warnings')
            else:
                sbml_warnings = False
            m = Model(sbml_filename=file_name, sbml_warnings=sbml_warnings)
            # Uncomment if you want a bioscrape XML written as well.
            # m.write_bioscrape_xml('temp_bs'+ file_name + '.xml')
            if initial_condition_dict is not None:
                processed = process_initial_concentration_dict(
                    initial_condition_dict
                )
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
        """Simulate CRN with libRoadRunner.

        Converts the CRN to SBML in memory and simulates it using the
        libRoadRunner deterministic simulator. libRoadRunner is a fast
        SBML simulator for deterministic (ODE) simulation.

        Parameters
        ----------
        timepoints : list of float
            Array of time points at which to record simulation results.
        initial_condition_dict : dict, optional
            Dictionary mapping species names (strings) to initial
            concentrations. Overrides the CRN's `initial_concentration_dict`.
        return_roadrunner : bool, default=False
            If True, returns the RoadRunner model object instead of simulation
            results. Useful for advanced control and analysis.
        check_validity : bool, default=True
            If True, validates the CRN before generating SBML.

        Returns
        -------
        array or RoadRunner
            If `return_roadrunner=False`, returns simulation results as a
            numpy array. If `return_roadrunner=True`, returns the RoadRunner
            model object.

        Warns
        -----
        UserWarning
            Issues a warning if libroadrunner is not installed.

        Notes
        -----
        Requires libroadrunner to be installed: `pip install libroadrunner`

        libRoadRunner documentation: https://libroadrunner.org/

        Examples
        --------
        >>> result = crn.simulate_with_roadrunner(
        ...     timepoints=np.linspace(0, 10, 100)
        ... )
        >>> # Get RoadRunner object for advanced control
        >>> rr = crn.simulate_with_roadrunner(
        ...     timepoints=np.linspace(0, 10, 100),
        ...     return_roadrunner=True
        ... )
        >>> # Run parameter scan with RoadRunner
        >>> result = rr.simulate(0, 10, 100)

        """
        res_ar = None
        try:
            import io

            import roadrunner  # type: ignore

            document, _ = self.generate_sbml_model(
                stochastic_model=False, check_validity=check_validity
            )
            sbml_string = libsbml.writeSBMLToString(document)
            # write the sbml_string into a temporary file in memory instead of
            # a file
            string_out = io.StringIO()
            string_out.write(sbml_string)
            # use the temporary file in memory to load the model into
            # libroadrunner
            rr = roadrunner.RoadRunner(string_out.getvalue())
            if initial_condition_dict:
                for species, value in initial_condition_dict.items():
                    rr.model[f"init([{species}])"] = value

            if return_roadrunner:
                return rr
            else:
                result = rr.simulate(
                    timepoints[0], timepoints[-1], len(timepoints)
                )
                res_ar = result
        except ModuleNotFoundError:
            warnings.warn(
                "libroadrunner was not found, please install libroadrunner"
            )
        return res_ar
