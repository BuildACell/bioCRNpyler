#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.


import copy
import itertools
from typing import List, Union
from warnings import warn

from ..utils.general import remove_bindloc
from .propensities import MassAction, Propensity
from .species import Species, WeightedSpecies


class Reaction(object):
    r"""Chemical reaction in a CRN with species and rate law.

    A `Reaction` represents a chemical transformation between species with
    an associated propensity (rate law). Reactions can be irreversible or
    reversible, and support various kinetic models through different
    propensity types.

    Parameters
    ----------
    inputs : list of Species or list of WeightedSpecies
        Reactant species. Can be Species objects (stoichiometry=1) or
        WeightedSpecies objects (with custom stoichiometry). Duplicates
        are automatically combined.
    outputs : list of Species or list of WeightedSpecies
        Product species. Can be Species objects (stoichiometry=1) or
        WeightedSpecies objects (with custom stoichiometry). Duplicates
        are automatically combined.
    propensity_type : Propensity
        Propensity object defining the rate law (e.g., MassAction, Hill).

    Attributes
    ----------
    inputs : list of WeightedSpecies
        Reactant species with stoichiometry.
    outputs : list of WeightedSpecies
        Product species with stoichiometry.
    propensity_type : Propensity
        The rate law for this reaction.
    is_reversible : bool
        True if the propensity supports reversible kinetics.
    species : list of Species
        All species involved in the reaction (inputs, outputs, and
        propensity species).

    See Also
    --------
    Species : Chemical species in a CRN.
    WeightedSpecies : Species with stoichiometric coefficient.
    Propensity : Base class for rate laws.
    MassAction : Mass action kinetics propensity.

    Notes
    -----
    A reaction has the form:

    .. math::
        \sum_i n_i I_i \rightarrow \sum_i m_i O_i

    where :math:`n_i` is the stoichiometry of reactant :math:`I_i` and
    :math:`m_i` is the stoichiometry of product :math:`O_i`.

    For reversible reactions:

    .. math::
        \sum_i n_i I_i \rightleftharpoons \sum_i m_i O_i

    Stoichiometry is handled as follows:

    - Species lists automatically combine duplicates
    - A + A --> B becomes 2A --> B
    - Stoichiometry affects rate calculations in mass action kinetics

    Different propensity types implement different rate laws:

    - MassAction: Standard mass action kinetics
    - Hill functions: Cooperative binding kinetics
    - GeneralPropensity: Custom formula

    Examples
    --------
    Create a simple irreversible reaction:

    >>> A = bcp.Species('A')
    >>> B = bcp.Species('B')
    >>> prop = bcp.MassAction(k_forward=0.1)
    >>> rxn = bcp.Reaction([A], [B], prop)

    Create a reversible reaction:

    >>> C = bcp.Species('C')
    >>> prop = bcp.MassAction(k_forward=100.0, k_reverse=0.01)
    >>> rxn = bcp.Reaction([A, B], [C], prop)
    >>> rxn.is_reversible
    True

    Create a reaction with stoichiometry:

    >>> rxn = bcp.Reaction([A, A], [B], prop)  # 2A <--> B

    Use the from_massaction class method:

    >>> rxn = bcp.Reaction.from_massaction([A, B], [C], k_forward=100.0)

    """

    def __init__(
        self,
        inputs: Union[List[Species], List[WeightedSpecies]],
        outputs: Union[List[Species], List[WeightedSpecies]],
        propensity_type: Propensity,
    ):
        if len(inputs) == 0 and len(outputs) == 0:
            warn("Reaction Inputs and Outputs both contain 0 Species")

        self.inputs = remove_bindloc(Species.flatten_list(inputs))
        self.outputs = remove_bindloc(Species.flatten_list(outputs))
        self.propensity_type = propensity_type

    @property
    def propensity_type(self) -> Propensity:
        """Propensity: The rate law for this reaction."""
        return self._propensity_type

    @propensity_type.setter
    def propensity_type(self, new_propensity_type: Propensity):
        """Set the propensity type for the reaction.

        Parameters
        ----------
        new_propensity_type : Propensity
            New propensity object. Must be a valid Propensity subclass
            instance.

        Raises
        ------
        ValueError
            If `new_propensity_type` is not a valid Propensity instance.

        """
        if not Propensity.is_valid_propensity(new_propensity_type):
            raise ValueError(
                f"unknown propensity type: {new_propensity_type} "
                f"({type(new_propensity_type)})!"
            )

        self._propensity_type = new_propensity_type

    @classmethod
    def from_massaction(
        cls,
        inputs: Union[List[Species], List[WeightedSpecies]],
        outputs: Union[List[Species], List[WeightedSpecies]],
        k_forward: float,
        k_reverse: float = None,
    ):
        """Create a Reaction with mass action kinetics.

        Convenience constructor for creating reactions with MassAction
        propensity.

        Parameters
        ----------
        inputs : list of Species or list of WeightedSpecies
            Reactant species.
        outputs : list of Species or list of WeightedSpecies
            Product species.
        k_forward : float
            Forward reaction rate constant. Must be positive.
        k_reverse : float, optional
            Reverse reaction rate constant. If None, reaction is
            irreversible. If provided, must be positive.

        Returns
        -------
        Reaction
            New Reaction object with MassAction propensity.

        Examples
        --------
        Create an irreversible reaction:

        >>> A = bcp.Species('A')
        >>> B = bcp.Species('B')
        >>> rxn = bcp.Reaction.from_massaction([A], [B], k_forward=0.1)

        Create a reversible reaction:

        >>> rxn = bcp.Reaction.from_massaction(
        ...     [A, B], [C], k_forward=100.0, k_reverse=0.01)

        """
        mak = MassAction(k_forward=k_forward, k_reverse=k_reverse)

        return cls(inputs=inputs, outputs=outputs, propensity_type=mak)

    @property
    def is_reversible(self) -> bool:
        """bool: True if the reaction has reversible kinetics.

        Determined by the propensity type's is_reversible property.

        """
        return self.propensity_type.is_reversible

    @property
    def inputs(self) -> List[WeightedSpecies]:
        """List of WeightedSpecies: Reactant species with stoichiometry."""
        return self._input_complexes

    @inputs.setter
    def inputs(self, new_input_complexes: List[WeightedSpecies]):
        """Set the reaction inputs.

        Parameters
        ----------
        new_input_complexes : list of Species or list of WeightedSpecies
            New reactant species. Species are automatically converted to
            WeightedSpecies. Duplicates are combined with adjusted
            stoichiometry.

        """
        self._input_complexes = Reaction._check_and_convert_complex_list(
            complexes=new_input_complexes
        )

    @property
    def outputs(self) -> List[WeightedSpecies]:
        """List of WeightedSpecies: Product species with stoichiometry."""
        return self._output_complexes

    @outputs.setter
    def outputs(self, new_output_complexes: List[WeightedSpecies]):
        """Set the reaction outputs.

        Parameters
        ----------
        new_output_complexes : list of Species or list of WeightedSpecies
            New product species. Species are automatically converted to
            WeightedSpecies. Duplicates are combined with adjusted
            stoichiometry.

        """
        self._output_complexes = Reaction._check_and_convert_complex_list(
            complexes=new_output_complexes
        )

    @staticmethod
    def _check_and_convert_complex_list(
        complexes: Union[List[Species], List[WeightedSpecies]],
    ) -> List[WeightedSpecies]:
        """Convert and validate species list to WeightedSpecies list.

        Converts Species to WeightedSpecies, validates all elements are
        proper types, and combines duplicate species by summing
        stoichiometry.

        Parameters
        ----------
        complexes : list of Species or list of WeightedSpecies
            Input species list to convert and validate.

        Returns
        -------
        list of WeightedSpecies
            Converted list with duplicates combined. Each unique species
            appears once with total stoichiometry.

        Raises
        ------
        TypeError
            If `complexes` contains elements that are neither Species nor
            WeightedSpecies.

        Notes
        -----
        Duplicate handling examples:

        - [A, A, B] --> [2A, B]
        - [WeightedSpecies(A, 2), A] --> [WeightedSpecies(A, 3)]

        """
        if all(
            [isinstance(one_complex, Species) for one_complex in complexes]
        ):
            # we wrap each Species object to WeightedSpecies
            complexes = [
                WeightedSpecies(species=species) for species in complexes
            ]
        else:
            if not all(
                [
                    isinstance(one_complex, WeightedSpecies)
                    for one_complex in complexes
                ]
            ):
                raise TypeError(
                    f"inputs must be list of Species or list of "
                    f"ChemicalComplexes! Recieved {complexes}"
                )

        # filter out duplicates and adjust stoichiometry
        out_list = []
        # Create a dictionary of unique species and their stoichiometry count
        stoichiometry_count = WeightedSpecies._count_weighted_species(
            complexes
        )
        for one_complex, stoichiometry in stoichiometry_count.items():
            new_complex = WeightedSpecies(
                species=one_complex.species, stoichiometry=stoichiometry
            )
            out_list.append(new_complex)

        return out_list

    # @property
    # def k_forward(self):
    #    return self.propensity_type.k_forward

    # @property
    # def k_reverse(self):
    #    return self.propensity_type.k_reverse

    def replace_species(self, species: Species, new_species: Species):
        """Create new reaction with a species replaced.

        Replaces all occurrences of a species throughout the reaction
        (inputs, outputs, and propensity species) with a new species.

        Parameters
        ----------
        species : Species
            Species to be replaced.
        new_species : Species
            Species to replace with.

        Returns
        -------
        Reaction
            New Reaction object with species replaced. The original
            reaction is not modified.

        Raises
        ------
        ValueError
            If either argument is not a Species object.

        Examples
        --------
        >>> A = bcp.Species('A')
        >>> B = bcp.Species('B')
        >>> C = bcp.Species('C')
        >>> rxn = bcp.Reaction.from_massaction([A, B], [C], k_forward=0.1)
        >>> D = bcp.Species('D')
        >>> rxn2 = rxn.replace_species(A, D)  # D + B --> C

        """
        if not isinstance(species, Species) or not isinstance(
            new_species, Species
        ):
            raise ValueError(
                "both species and new_species argument must be an "
                "instance of Species!"
            )

        new_inputs = []
        for s in self.inputs:
            new_s = s.replace_species(species, new_species)
            new_inputs.append(new_s)

        new_outputs = []
        for s in self.outputs:
            new_s = s.replace_species(species, new_species)
            new_outputs.append(new_s)

        # get a shallow copy of the parameters and species, so we can
        # replace some of them
        propensity_type_dict = copy.copy(self.propensity_type.propensity_dict)
        for key, prop_species in propensity_type_dict['species'].items():
            propensity_type_dict['species'][key] = (
                prop_species.replace_species(species, new_species)
            )

        new_propensity_type = self.propensity_type.from_dict(
            propensity_type_dict
        )

        new_r = Reaction(
            inputs=new_inputs,
            outputs=new_outputs,
            propensity_type=new_propensity_type,
        )
        return new_r

    def __repr__(self):
        """Return string representation of the reaction.

        Returns
        -------
        str
            Reaction equation showing inputs, outputs, material types, and
            attributes, but not rates or parameters.

        """
        # Helper function to print the text of a rate function.
        return self.pretty_print(
            show_rates=False,
            show_material=True,
            show_attributes=True,
            show_parameters=False,
        )

    def pretty_print(
        self,
        show_rates=True,
        show_material=True,
        show_attributes=True,
        show_parameters=True,
        **kwargs,
    ):
        """Generate detailed, formatted string representation of reaction.

        Parameters
        ----------
        show_rates : bool, default=True
            If True, includes rate law formula below reaction equation.
        show_material : bool, default=True
            If True, shows species material types (e.g., 'dna', 'protein').
        show_attributes : bool, default=True
            If True, shows species attributes.
        show_parameters : bool, default=True
            If True, shows parameter values in rate law.
        **kwargs
            Additional keyword arguments passed to species and propensity
            pretty_print methods. Can include 'stochastic' (bool) for
            stochastic vs deterministic rate display.

        Returns
        -------
        str
            Formatted reaction string. Format:
            'inputs --> outputs' or 'inputs <--> outputs' (reversible)
            Optionally followed by rate law and parameters.

        Examples
        --------
        >>> A = bcp.Species('A')
        >>> B = bcp.Species('B')
        >>> rxn = bcp.Reaction.from_massaction([A], [B], k_forward=0.1)
        >>> print(rxn.pretty_print())
        A --> B
         Kf=k_forward * A
          k_forward=0.1        Kf = 0.1 * A

        >>> print(rxn.pretty_print(show_rates=False))
        A --> B

        """
        kwargs['show_rates'] = show_rates
        kwargs['show_material'] = show_material
        kwargs['show_attributes'] = show_attributes

        txt = '+'.join([s.pretty_print(**kwargs) for s in self.inputs])

        if self.is_reversible:
            txt += ' <--> '
        else:
            txt += ' --> '

        txt += '+'.join([s.pretty_print(**kwargs) for s in self.outputs])
        if show_rates:
            # These kwargs are essential for massaction propensities
            kwargs['reaction'] = self
            if 'stochastic' not in kwargs:
                kwargs['stochastic'] = False

            txt += '\n' + f"{self.propensity_type.pretty_print(**kwargs)}"

        return txt

    def __eq__(self, other):
        """Test equality between reactions.

        Two reactions are equal if they have the same inputs, outputs,
        and propensity (in any order).

        Parameters
        ----------
        other : Reaction
            Another reaction to compare with.

        Returns
        -------
        bool
            True if reactions have identical inputs, outputs, and
            propensity.

        Raises
        ------
        TypeError
            If `other` is not a Reaction object.

        Notes
        -----
        Order of species in inputs/outputs doesn't matter:

        :math:`A + B \rightarrow C` equals :math:`B + A \rightarrow C`
        
        since species are compared using sets.

        """
        # Check if reactions are equivalent.
        #
        # Two reactions are equivalent if they have the same inputs,
        # outputs, and propensity.
        if not isinstance(other, Reaction):
            raise TypeError(
                f"Only reactions can be compared with reaction! "
                f"We got {type(other)}."
            )

        if len(self.inputs) != len(other.inputs) or len(self.outputs) != len(
            other.outputs
        ):
            return False

        return (
            set(self.inputs),
            set(self.outputs),
            self.propensity_type,
        ) == (set(other.inputs), set(other.outputs), other.propensity_type)

    def __contains__(self, item: Species):
        """Check if a species is involved in the reaction.

        Checks whether a species appears in the reaction's inputs,
        outputs, or propensity (e.g., Hill kinetics with regulatory
        species).

        Parameters
        ----------
        item : Species
            Species to check for.

        Returns
        -------
        bool
            True if species is in inputs, outputs, or propensity species.

        Raises
        ------
        NotImplementedError
            If `item` is not a Species object.

        Examples
        --------
        >>> A = bcp.Species('A')
        >>> B = bcp.Species('B')
        >>> C = bcp.Species('C')
        >>> rxn = bcp.Reaction.from_massaction([A], [B], k_forward=0.1)
        >>> A in rxn
        True
        >>> C in rxn
        False

        """
        if isinstance(item, Species):
            if (
                item in self.inputs
                or item in self.outputs
                or item in self.propensity_type.species
            ):
                return True
        else:
            raise NotImplementedError
        return False

    @property
    def species(self) -> List[Species]:
        """List of Species: All species involved in the reaction.

        Returns a flattened list of all species from inputs, outputs, and
        the propensity (e.g., Hill functions have regulatory species).

        Returns
        -------
        list of Species
            All species in the reaction. May contain duplicates if a
            species appears in multiple roles.

        Notes
        -----
        This property collects species from three sources:

        1. Input species (reactants)
        2. Output species (products)
        3. Propensity species (e.g., activators/repressors in Hill)

        Examples
        --------
        >>> A = bcp.Species('A')
        >>> B = bcp.Species('B')
        >>> rxn = bcp.Reaction.from_massaction([A], [B], k_forward=0.1)
        >>> rxn.species
        [A, B]

        """
        in_part = []
        for s in self.inputs:
            in_part.extend(Species.flatten_list(s.species))
        out_part = []
        for s in self.outputs:
            out_part.extend(Species.flatten_list(s.species))

        return list(
            itertools.chain(in_part, out_part, self.propensity_type.species)
        )
