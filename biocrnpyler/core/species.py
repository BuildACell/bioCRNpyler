#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
import warnings
from typing import List, Union

from .compartment import Compartment
from .polymer import MonomerCollection, OrderedMonomer, OrderedPolymer


class Species(OrderedMonomer):
    """A formal species object for a chemical reaction network (CRN).

    Represents a chemical species in a CRN with a name, material type,
    attributes, and compartment. Species inherits from `OrderedMonomer`,
    allowing it to be part of polymer structures while also functioning as
    an independent chemical entity in reactions.

    Parameters
    ----------
    name : str
        Name of the species. Must consist of letters, numbers, or
        underscores, cannot contain double underscores, and cannot
        begin/end with special characters.
    material_type : str, default=''
        Type of material (e.g., 'dna', 'rna', 'protein', 'complex').
        Required if name starts with a number. Must start with a letter.
    attributes : list of str or None, optional
        List of attribute tags for the species (e.g., 'degraded',
        'phosphorylated'). Each attribute must be alphanumeric.
    compartment : Compartment, str, or None, optional
        The compartment containing this species. If None, uses default
        compartment. If str, creates a new Compartment with that name.
    **kwargs
        Additional keyword arguments passed to `OrderedMonomer`.

    Attributes
    ----------
    name : str
        The name of the species.
    material_type : str
        The material type of the species.
    attributes : list of str
        List of attribute tags associated with the species.
    compartment : Compartment
        The compartment containing this species.
    direction : str, int, or None
        Directional orientation (inherited from `OrderedMonomer`). When
        set, the direction is also added as an attribute.

    See Also
    --------
    ComplexSpecies : Species formed from multiple bound species.
    OrderedPolymerSpecies : Polymer species for chemical reactions.
    WeightedSpecies : Species with stoichiometry coefficient.

    Notes
    -----
    Species names must:

    - Contain only letters, numbers, and underscores
    - Not contain double underscores ('__')
    - Not end with an underscore
    - Start with a letter or number (if starting with number, requires
      material_type)

    Species are represented as strings in the format:

    `material_type_name_attribute1_attribute2_compartment`

    Components are omitted if empty or default values.

    Two species are equal if they have the same name, material_type,
    attributes, compartment, parent, and position.

    Examples
    --------
    Create a simple species:

    >>> S = bcp.Species('S')
    >>> S.name
    'S'

    Create a protein with attributes:

    >>> GFP = bcp.Species(
    ...     name='GFP',
    ...     material_type='protein',
    ...     attributes=['fluorescent', 'degraded']
    ... )
    >>> repr(GFP)
    'protein_GFP_fluorescent_degraded'

    Create a species in a compartment:

    >>> cytoplasm = bcp.Compartment('cytoplasm')
    >>> enzyme = bcp.Species(
    ...     name='enzyme',
    ...     material_type='protein',
    ...     compartment=cytoplasm
    ... )

    """

    def __init__(
        self,
        name: str,
        material_type='',
        attributes: Union[List, None] = None,
        compartment=None,
        **kwargs,
    ):
        OrderedMonomer.__init__(self, **kwargs)

        self.name = name
        self.material_type = material_type
        self._attributes = []  # Set this to avoid errors
        self.attributes = attributes
        self.compartment = compartment

    @property
    def attributes(self):
        if not hasattr(self, '_attributes'):
            self._attributes = []
        return self._attributes

    @attributes.setter
    def attributes(self, attributes):
        if not hasattr(self, '_attributes'):
            self._attributes = []
        if attributes is not None:
            if not isinstance(attributes, list):
                attributes = list(attributes)
            for attribute in attributes:
                self.add_attribute(attribute)
        elif attributes is None:
            self._attributes = []

    def remove_attribute(self, attribute: str):
        """Remove an attribute from the species.

        Parameters
        ----------
        attribute : str
            The attribute to remove. Must be an alphanumeric string.

        Notes
        -----
        If the attribute is not present or is None, no action is taken.
        All occurrences of the attribute are removed from the attributes
        list.

        """
        if not hasattr(self, '_attributes') or attribute is None:
            return
        else:
            assert isinstance(attribute, str) and attribute.isalnum(), (
                f"Attribute: {attribute} must be an alpha-numeric string"
            )
            self._attributes = [a for a in self.attributes if a != attribute]

    def add_attribute(self, attribute: str):
        """Add an attribute to the species.

        Parameters
        ----------
        attribute : str
            The attribute to add. Must be an alphanumeric string and
            non-None.

        Raises
        ------
        AssertionError
            If attribute is not an alphanumeric string or is None.

        Notes
        -----
        Duplicate attributes are not added - each attribute appears only
        once in the attributes list.

        Examples
        --------
        >>> species = bcp.Species('MyProtein')
        >>> species.add_attribute('degraded')
        >>> species.attributes
        ['degraded']

        """
        if not hasattr(self, '_attributes'):
            self._attributes = []
        assert (
            isinstance(attribute, str)
            and attribute is not None
            and attribute.isalnum()
        ), f"Attribute: {attribute} must be an alpha-numeric string"
        if attribute not in self.attributes:
            self._attributes.append(attribute)

    @property
    def name(self):
        if self._name is None:
            return ''
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        if name is None:
            raise TypeError(f"Name must be a string. Recievied {name}")
        else:
            self._name = self._check_name(name)

    @property
    def compartment(self):
        return self._compartment

    @compartment.setter
    def compartment(self, compartment):
        if compartment is None:
            self._compartment = Compartment(name='default')
        else:
            if isinstance(compartment, str):
                self._compartment = Compartment(
                    name=self._check_name(compartment)
                )
            elif isinstance(compartment, Compartment):
                self._compartment = compartment

    @property
    def direction(self):
        if hasattr(self, '_direction'):
            return self._direction
        else:
            return None

    @direction.setter
    def direction(self, direction):
        """Set the directional orientation of the species.

        Overrides `OrderedMonomer.direction` to automatically add the
        direction as an attribute and remove the old direction attribute.

        Parameters
        ----------
        direction : str, int, or None
            The direction to assign. Common values include 'forward',
            'reverse', 0, 1, or None. When set, the direction is added as
            an attribute.

        Notes
        -----
        This is inherited from `OrderedMonomer`. A species with direction
        will use it as an attribute as well. This is overwritten to make
        direction an attribute.

        """
        # Remove old direction from attributes
        self.remove_attribute(self.direction)
        # set the new direction
        self._direction = direction
        # Add it to attributes
        if direction is not None:
            self.add_attribute(direction)

    def remove(self):
        """Remove the species from its parent polymer.

        Overrides `OrderedMonomer.remove` to also remove the direction
        attribute if present.

        Returns
        -------
        Species
            Returns self after removal for method chaining.

        """
        if self.direction is not None:
            self.remove_attribute(self.direction)
        return OrderedMonomer.remove(self)  # call the OrderedMonomer function

    # Note: this is used because properties can't be overwritten without
    # setters being overwritten in subclasses.
    def _check_name(self, name):
        """Validate that name follows proper formatting rules.

        Parameters
        ----------
        name : str or None
            The name to validate.

        Returns
        -------
        str or None
            The validated name, unchanged if valid.

        Raises
        ------
        ValueError
            If name violates formatting rules.
        TypeError
            If name is not a string or None.

        Notes
        -----
        Valid names must:

        - Contain only underscores and alphanumeric characters
        - Not contain double underscores ('__')
        - Not end with an underscore
        - Start with an alphanumeric character

        """
        if name is None:
            return name
        elif isinstance(name, str):
            no_underscore_string = name.replace('_', '')
            if (
                no_underscore_string.isalnum()
                and '__' not in name
                and name[len(name) - 1] != '_'
                and name[0].isalnum()
            ):
                return name
            else:
                raise ValueError(
                    f"name attribute {name} must consist of letters, "
                    "numbers, or underscores and cannot contain double "
                    "underscores or begin/end with a special character."
                )
        else:
            raise TypeError("Name must be a string.")

    @property
    def material_type(self):
        return self._material_type

    @material_type.setter
    def material_type(self, material_type: str):
        """Material type for the species.

        Check that the string contains is alpha-numeric characters or "_"
        and that the first character is a letter.  If the name is a starts
        with a number, there must be a material type.

        """
        if material_type in [None, ''] and self.name[0].isnumeric():
            raise ValueError(
                f"species name: {self.name} contains a number as the "
                f"first character and therefore requires a material_type."
            )
        elif material_type is None:
            self._material_type = None
        elif (
            material_type.replace('_', '').isalnum()
            and material_type.replace('_', '')[0].isalpha()
            and '__' not in material_type
            and material_type[len(material_type) - 1] != '_'
        ) or material_type == '':
            self._material_type = material_type
        else:
            raise ValueError(
                f"material_type {material_type} must be alpha-numeric "
                f"and start with a letter."
            )

    def __repr__(self):
        txt = ''
        if self.material_type not in ['', None]:
            txt = self.material_type + '_'

        txt += self.name

        if len(self.attributes) > 0 and self.attributes != []:
            for i in self.attributes:
                if i is not None:
                    txt += '_' + str(i)
        if (
            self.compartment is not None
            and self.compartment.name != 'default'
        ):
            # Only add a compartment name if it is not the default one.  if
            # compartment name is already there with an underscore remove it
            # from the string first to not repeat the compartment tag
            txt = txt.replace('_' + self.compartment.name, '')
            txt += '_' + self.compartment.name
        txt = txt.replace("'", '')
        return txt

    def replace_species(self, species, new_species):
        """Replace a species with another species.

        For a simple Species, returns `new_species` if this species equals
        `species`, otherwise returns self. For complex species, acts
        recursively.

        Parameters
        ----------
        species : Species
            The species to search for and replace.
        new_species : Species
            The species to replace with.

        Returns
        -------
        Species
            Either `new_species` (if self == species) or self.

        Raises
        ------
        ValueError
            If either argument is not a Species instance.

        """
        if not isinstance(species, Species):
            raise ValueError(
                "species argument must be an instance of Species!"
            )

        if not isinstance(new_species, Species):
            raise ValueError(
                "species argument must be an instance of Species!"
            )

        if self == species:
            return new_species
        else:
            return self

    def get_species(self, recursive=None):
        """Get a list containing this species.

        Returns
        -------
        list of Species
            A list containing only this species: [self].

        Notes
        -----
        This method is used in recursive calls where `ComplexSpecies`
        returns a list of constituent species while `Species` returns just
        itself in a list. The `recursive` parameter is accepted for
        compatibility but not used in the base Species class.

        """
        return [self]

    def pretty_print(
        self,
        show_material=True,
        show_compartment=False,
        show_attributes=True,
        show_initial_condition=False,
        **kwargs,  # TODO: allows spurious keywords; fix...
    ):
        """Generate a human-readable string representation of the species.

        Parameters
        ----------
        show_material : bool, default=True
            If True, includes material_type in brackets around the species.
        show_compartment : bool, default=False
            If True, shows the compartment name in the representation.
        show_attributes : bool, default=True
            If True, includes attributes in parentheses after the name.
        show_initial_condition : bool, default=False
            Placeholder for compatibility with CRN printing.
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        str
            Formatted string representation of the species.

        Notes
        -----
        This method provides more detailed output than `__repr__`,
        useful for understanding CRNs but does not return string
        identifiers compatible with parsers.

        Format: `material_type[name(attr1, attr2)-direction]`

        Examples
        --------
        >>> S = bcp.Species('S', material_type='protein',
        ...                 attributes=['active'])
        >>> S.pretty_print()
        'protein[S(active)]'

        """
        txt = ''
        if self.material_type not in ['', None] and show_material:
            txt = self.material_type + '['

        txt += self.name

        if self.compartment not in ['default', None] and show_compartment:
            txt += ' in ' + self.compartment.name + '.'

        if (
            len(self.attributes) > 0
            and self.attributes != []
            and show_attributes
        ):
            txt += '('
            for i in self.attributes:
                if i is not None:
                    txt += str(i) + ', '
            txt = txt[:-2] + ')'

        txt = txt.replace("'", '')
        if hasattr(self, 'direction') and self.direction is not None:
            txt += '-' + self.direction
        if self.material_type not in ['', None] and show_material:
            txt += ']'

        return txt

    def __eq__(self, other):
        """Check if two species are equivalent.

        Two species are equal if they have the same name, material_type,
        attributes (as sets), parent, compartment, and position.

        Parameters
        ----------
        other : Species
            The species to compare with.

        Returns
        -------
        bool
            True if species are equivalent, False otherwise.

        Notes
        -----
        Equality between parents and children can result in loops, so
        string equality is used for parent comparison.

        """
        if (
            isinstance(other, Species)
            and self.material_type == other.material_type
            and self.name == other.name
            and set(self.attributes) == set(other.attributes)
            and str(self.parent) == str(other.parent)
            and self.compartment == other.compartment
            and self.position == other.position
        ):
            return True
        else:
            return False

    def monomer_eq(self, other):
        """Check if two monomers are equal, ignoring parent and position.

        Parameters
        ----------
        other : Species
            The species to compare with.

        Returns
        -------
        bool
            True if species have the same name, material_type, attributes,
            and compartment, regardless of parent or position.

        Notes
        -----
        This is the same as normal equality but does not check for parents
        or positions. Useful for comparing species that may be in different
        polymer contexts.

        """
        if (
            isinstance(other, Species)
            and self.material_type == other.material_type
            and self.name == other.name
            and set(self.attributes) == set(other.attributes)
            and self.compartment == other.compartment
        ):
            return True
        else:
            return False

    def __gt__(self, Species2):
        return self.name > Species2.name

    def __lt__(self, Species2):
        return self.name < Species2.name

    def __hash__(self):
        return str.__hash__(repr(self))

    def __contains__(self, item):
        return item in self.get_species()

    def contains_species_monomer(self, s):
        """Check if the species contains a monomer, ignoring context.

        Parameters
        ----------
        s : Species
            The species monomer to search for.

        Returns
        -------
        bool
            True if the species contains a monomer equal to `s` (ignoring
            parent, position, and direction), False otherwise.

        Notes
        -----
        This is a less stringent version of `__contains__` that checks
        without considering Species.parent, Species.position, or direction.
        Useful for determining if a species is present regardless of its
        polymer context.

        """
        s_copy = copy.deepcopy(s)
        s_copy.remove()
        for ss in self.get_species(recursive=True):
            ss_copy = copy.copy(ss)
            ss_copy.remove()
            if ss_copy == s_copy:
                return True
        return False

    @staticmethod
    def flatten_list(in_list) -> List:
        """Recursively flatten a nested list of species.

        Parameters
        ----------
        in_list : list or Species
            A potentially nested list of species, or a single species.

        Returns
        -------
        list
            Flattened list containing all species. None elements are
            filtered out.

        Examples
        --------
        >>> S1 = bcp.Species('S1')
        >>> S2 = bcp.Species('S2')
        >>> nested = [S1, [S2, None]]
        >>> bcp.Species.flatten_list(nested)
        [S1, S2]

        """
        out_list = []
        if not isinstance(in_list, list):
            out_list.append(in_list)
        else:
            for element in in_list:
                if isinstance(element, list):
                    out_list += Species.flatten_list(element)
                elif element is None:
                    pass
                else:
                    out_list += [element]
        return out_list


class WeightedSpecies:
    """Container for a species with stoichiometric coefficient.

    Wraps a `Species` object together with its stoichiometry for use in
    reactions. This class is primarily used internally by the Reaction
    class to represent reactants and products with their coefficients.

    Parameters
    ----------
    species : Species
        The species object.
    stoichiometry : int, default=1
        The stoichiometric coefficient. Must be a positive integer.

    Attributes
    ----------
    species : Species
        The wrapped species object.
    stoichiometry : int
        The stoichiometric coefficient (positive integer).

    See Also
    --------
    Species : Base class for chemical species.
    Reaction : Chemical reaction containing weighted species.

    Examples
    --------
    Create a weighted species:

    >>> S = bcp.Species('S')
    >>> ws = bcp.WeightedSpecies(species=S, stoichiometry=2)
    >>> ws.stoichiometry
    2

    """

    def __init__(self, species: Species, stoichiometry: int = 1):
        self.species: Species = species
        self.stoichiometry: int = stoichiometry

    @property
    def stoichiometry(self):
        return self._stoichiometry

    @stoichiometry.setter
    def stoichiometry(self, new_stoichiometry):
        if new_stoichiometry <= 0:
            raise ValueError(
                f"Stoichiometry must be positive integer! "
                f"We got {new_stoichiometry}!"
            )
        self._stoichiometry = int(new_stoichiometry)

    def pretty_print(self, **kwargs):
        return (
            f"{self.stoichiometry if self.stoichiometry > 1 else ''}"
            + f"{self.species.pretty_print(**kwargs)}"
        )

    def replace_species(self, *args, **kwargs):
        return self.species.replace_species(*args, **kwargs)

    @staticmethod
    def _count_weighted_species(weighted_species: List):
        """Merge species in a list with different stoichiometry.

        Combines `WeightedSpecies` objects with the same species by summing
        their stoichiometric coefficients.

        Parameters
        ----------
        weighted_species : list of WeightedSpecies
            List of weighted species to merge.

        Returns
        -------
        dict
            Dictionary mapping unique `WeightedSpecies` to their total
            stoichiometry.

        Examples
        --------
        >>> s1 = bcp.Species(name='a')
        >>> ws1 = bcp.WeightedSpecies(species=s1, stoichiometry=2)
        >>> ws2 = bcp.WeightedSpecies(species=s1, stoichiometry=5)
        >>> ws_list = [ws1, ws2]
        >>> freq_dict = bcp.WeightedSpecies._count_weighted_species(ws_list)
        >>> len(freq_dict)
        1

        """
        # convert to set doesn't work because we need only species equality
        unique_species = []
        for w_species in weighted_species:
            if not any(
                w_species.species == u_s.species for u_s in unique_species
            ):
                unique_species.append(w_species)

        freq_dict = dict(zip(unique_species, [0] * len(unique_species)))
        for w_species in weighted_species:
            for key in freq_dict:
                if key.species == w_species.species:
                    freq_dict[key] += w_species.stoichiometry

        return freq_dict

    def __eq__(self, other):
        if other.__class__ is self.__class__:
            return (other.species, other.stoichiometry) == (
                self.species,
                self.stoichiometry,
            )
        return False

    def __hash__(self):
        return hash(self.species) + hash(self.stoichiometry)


class ComplexSpecies(Species):
    """Species formed from multiple bound species.

    A special kind of species representing a complex of two or more bound
    species. ComplexSpecies should always be created using the `Complex`
    function, not directly. Order of species in the list does not matter:
    ComplexSpecies([s1, s2]) == ComplexSpecies([s2, s1]).

    Parameters
    ----------
    species : list of Species or str
        List of species forming the complex. Must contain at least 2
        species.
    name : str or None, optional
        Custom name for the complex. If None, generates a name from
        constituent species.
    material_type : str, default='complex'
        Material type identifier for the complex.
    attributes : list of str or None, optional
        Attributes for the complex species.
    compartment : Compartment, str, or None, optional
        Compartment containing the complex.
    called_from_complex : bool, default=False
        Internal flag to enforce use of `Complex` function.

    Attributes
    ----------
    species : list of Species
        Sorted list of constituent species in the complex.
    species_set : list of Species
        Unique species in the complex, sorted by string representation.
    name : str
        Name of the complex (auto-generated if not provided).

    See Also
    --------
    Complex : Metaclass for creating ComplexSpecies.
    OrderedComplexSpecies : Complex where species order matters.
    Species : Base class for chemical species.

    Notes
    -----
    ComplexSpecies add an additional '_' at the end of their string
    representation to differentiate edge cases.

    Species order does not affect equality: ComplexSpecies([s1, s2])
    == ComplexSpecies([s2, s1])

    For ordered complexes, use `OrderedComplexSpecies`.

    If no name is provided, the complex is named by concatenating all
    constituent species names with counts for duplicates.

    Always use the `Complex` function to create `ComplexSpecies`:

    >>> # Correct
    >>> complex_species = bcp.Complex([S1, S2])

    >>> # Incorrect (will raise warning)
    >>> complex_species = bcp.ComplexSpecies([S1, S2])

    Examples
    --------
    Create a complex (using Complex function):

    >>> S1 = bcp.Species('S1')
    >>> S2 = bcp.Species('S2')
    >>> complex_species = bcp.Complex([S1, S2])

    Check if a species is in a complex:

    >>> S1 in complex_species
    True

    """

    def __init__(
        self,
        species: List[Union[Species, str]],
        name: Union[str, None] = None,
        material_type='complex',
        attributes=None,
        compartment=None,
        called_from_complex=False,
    ):
        # A little check to enforce use of Complex() to create ComplexSpecies
        if not called_from_complex:
            warnings.warn(
                "ComplexSpecies should be created using the "
                "Complex([List of Species]) function, not directly!"
            )

        # Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError(
                "chemical_reaction_network.complex requires 2 "
                "or more species in its constructor."
            )
        self.species = species

        # call super class
        Species.__init__(
            self,
            name=name,
            material_type=material_type,
            attributes=attributes,
            compartment=compartment,
        )

    def __repr__(self):
        """Generate string representation of ComplexSpecies.

        Returns
        -------
        str
            String representation with an additional '_' at the end to
            differentiate edge cases.

        Notes
        -----
        ComplexSpecies add an additional '_' onto the end of their string
        representation. This ensures that some edge cases are
        differentiated.

        """
        txt = Species.__repr__(self)
        txt += '_'
        return txt

    @property
    def name(self):
        if self._name is None:
            name = ''
            for s in self.species_set:
                count = self.species.count(s)
                name += str(s) + '_'
                if count > 1:
                    name += f"{count}x_"
            name = name[:-1]
            return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    def __contains__(self, item):
        """Check if a species is contained in the complex.

        Parameters
        ----------
        item : Species
            The species to search for.

        Returns
        -------
        bool
            True if the species is found in the complex or any nested
            complexes, False otherwise.

        Raises
        ------
        ValueError
            If `item` is not a Species instance.

        Notes
        -----
        This method searches recursively through nested ComplexSpecies to
        find the target species.

        """
        if not isinstance(item, Species):
            raise ValueError(
                "Operator 'in' requires chemical_reaction_network.Species "
                f"(or a subclass). Received: {str(item)}"
            )
        if item in self.species:
            # this is the base case
            return True
        else:
            # this is the recursive part. We want to check all
            # internal complexes for the thing we're looking for
            for content in self.species:
                if isinstance(content, ComplexSpecies):
                    if item in content:
                        return True
            # if we got here then we've failed to find it
            return False

    @property
    def species_set(self):
        species_set = list(set(self.species))
        list.sort(species_set, key=lambda s: repr(s))
        return species_set

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if not isinstance(species, list):
            raise TypeError(f"species must be a list: recieved {species}.")
        species = Species.flatten_list(species)
        if not all(isinstance(s, Species) for s in species):
            raise TypeError(
                f"recieved a non-species as a member of the list species: "
                f"{species}."
            )
        else:
            list.sort(species, key=lambda s: repr(s))
            self._species = species

    def replace_species(self, species: Species, new_species: Species):
        """Replace a species throughout the entire complex.

        Acts recursively on nested ComplexSpecies. Does not modify in
        place - returns a new ComplexSpecies.

        Parameters
        ----------
        species : Species
            The species to replace.
        new_species : Species
            The species to replace with.

        Returns
        -------
        ComplexSpecies
            A new ComplexSpecies with the replacement applied.

        Raises
        ------
        ValueError
            If either argument is not a Species instance.

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
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self._name is not None:
            new_name = self.name

        return Complex(
            species=new_species_list,
            name=new_name,
            material_type=self.material_type,
            attributes=self.attributes,
        )

    def get_species(self, recursive=False):
        """Get all species in the complex.

        Parameters
        ----------
        recursive : bool, default=False
            If True, returns species inside nested ComplexSpecies
            recursively. If False, returns only this ComplexSpecies.

        Returns
        -------
        list of Species
            List of species. If recursive=False, returns [self]. If
            recursive=True, returns [self] plus all constituent species.

        """
        if not recursive:
            species = [self]
        else:
            species = [self]
            for s in self.species:
                species += s.get_species(recursive=True)

        return species

    def pretty_print(
        self,
        show_material=True,
        show_compartment=False,
        show_attributes=True,
        show_initial_condition=False,
        **kwargs,  # TODO: allows spurious keywords; fix...
    ):
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string identifiers.
        show_material toggles whether species.material is printed.
        show_attributes toggles whether species.attributes is printed

        """
        txt = ''
        if self.material_type not in ['', None] and show_material:
            txt += self.material_type
        txt += '['
        for s in self.species_set:
            count = self.species.count(s)
            if count > 1:
                txt += f"{count}x_"
            txt += (
                s.pretty_print(
                    show_material=show_material, show_attributes=False
                )
                + ':'
            )
        txt = txt[:-1]

        if self.compartment not in ['default', None] and show_compartment:
            txt += ' in ' + self.compartment.name + '.'

        if (
            len(self.attributes) > 0
            and self.attributes != []
            and show_attributes
        ):
            txt += '('
            for i in self.attributes:
                if i is not None:
                    txt += str(i) + ', '
            txt = txt[:-2] + ')'

        txt = txt.replace("'", '')
        if hasattr(self, 'direction') and self.direction is not None:
            txt += '-' + self.direction
        txt += ']'

        return txt

    def monomer_count(self, m):
        """Count occurrences of a monomer in the complex.

        Parameters
        ----------
        m : Species
            The monomer to count.

        Returns
        -------
        int
            Number of times the monomer appears in the complex, using
            `monomer_eq` for equality comparison.

        Notes
        -----
        This effectively implements `self.species.count(m)` but uses
        `monomer_eq` for equality, which ignores parent and position.

        """
        return sum([s.monomer_eq(m) for s in self.species])


class OrderedComplexSpecies(ComplexSpecies):
    """Complex species where species order is significant.

    A special kind of species formed from a complex of two or more species
    where the order matters. OrderedComplexSpecies should always be created
    using the `Complex` function with `ordered=True`, not directly.
    Unlike ComplexSpecies, [s1, s2, s3] != [s1, s3, s2].

    Parameters
    ----------
    species : list of Species or str
        Ordered list of species forming the complex. Must contain at least
        2 species.
    name : str or None, optional
        Custom name for the complex. If None, generates a name from
        constituent species in order.
    material_type : str, default='ordered_complex'
        Material type identifier for the ordered complex.
    attributes : list of str or None, optional
        Attributes for the complex species.
    compartment : Compartment, str, or None, optional
        Compartment containing the complex.
    called_from_complex : bool, default=False
        Internal flag to enforce use of `Complex` function.

    Attributes
    ----------
    species : list of Species
        Ordered list of constituent species (NOT sorted).
    name : str
        Name of the complex (auto-generated if not provided).

    See Also
    --------
    Complex : Metaclass for creating ordered complexes.
    ComplexSpecies : Complex where species order doesn't matter.
    OrderedPolymerSpecies : Ordered polymer for chemical reactions.

    Notes
    -----
    Unlike `ComplexSpecies`, the order of species matters:
    OrderedComplexSpecies([s1, s2]) != OrderedComplexSpecies([s2, s1])

    Similar to ComplexSpecies, OrderedComplexSpecies add an additional
    '_' at the end.

    Always use `Complex` with `ordered=True`:

    >>> # Correct
    >>> complex_species = bcp.Complex([S1, S2], ordered=True)

    >>> # Incorrect (will raise warning)
    >>> complex_species = bcp.OrderedComplexSpecies([S1, S2])

    Examples
    --------
    Create an ordered complex:

    >>> S1 = bcp.Species('S1')
    >>> S2 = bcp.Species('S2')
    >>> ordered = bcp.Complex([S1, S2], ordered=True)
    >>> reversed = bcp.Complex([S2, S1], ordered=True)
    >>> ordered == reversed
    False

    """

    def __init__(
        self,
        species,
        name=None,
        material_type='ordered_complex',
        attributes=None,
        compartment=None,
        called_from_complex=False,
    ):
        # A little check to enforce use of Complex() to create ComplexSpecies
        if not called_from_complex:
            warnings.warn(
                "OrderedComplexSpecies should be created using the "
                "Complex([List of Species]) function, not directly!"
            )

        # Set species because it is used for default naming
        if len(Species.flatten_list(species)) <= 1:
            raise ValueError(
                "chemical_reaction_network.complex requires 2 "
                "or more species in its constructor."
            )
        self.species = species

        # Call the Species superclass constructor
        Species.__init__(
            self,
            name=name,
            material_type=material_type,
            attributes=attributes,
            compartment=compartment,
        )

    @property
    def name(self):
        if self._name is None:
            name = ''
            for s in self.species:
                if isinstance(s, str):
                    s = Species(name=s)
                if s.material_type not in ['']:
                    name += f"{s.material_type}_{s.name}_"
                else:
                    name += f"{s.name}_"
            name = name[:-1]
            return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    @property
    def species(self):
        return self._species

    @species.setter
    def species(self, species):
        if not isinstance(species, list):
            raise TypeError(f"species must be a list: recieved {species}.")
        species = Species.flatten_list(species)
        if not all(isinstance(s, Species) for s in species):
            raise TypeError(
                f"recieved a non-species as a member of the list species: "
                f"{species}."
            )
        else:
            self._species = species

    def replace_species(self, species: Species, new_species: Species):
        """Replaces species with new_species in the entire Complex Species.

        Acts recursively on nested ComplexSpecies.

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
        for s in self.species:
            new_s = s.replace_species(species, new_species)
            new_species_list.append(new_s)

        new_name = None
        if self._name is not None:
            new_name = self.name

        return Complex(
            species=new_species_list,
            name=new_name,
            material_type=self.material_type,
            attributes=self.attributes,
            ordered=True,
        )

    def pretty_print(
        self,
        show_material=True,
        show_compartment=False,
        show_attributes=True,
        show_initial_condition=False,
        **kwargs,  # TODO: allows spurious keywords; fix...
    ):
        """A more powerful printing function.

        Useful for understanding CRNs but does not return string
        identifiers.  show_material toggles whether species.material is
        printed.  show_attributes toggles whether species.attributes is
        printed.

        """
        txt = ''
        if self.material_type not in ['', None] and show_material:
            txt += self.material_type

        txt += '['

        for s in self.species:
            txt += (
                s.pretty_print(
                    show_material=show_material, show_attributes=False
                )
                + ':'
            )
        txt = txt[:-1]

        if self.compartment not in ['default', None] and show_compartment:
            txt += ' in ' + self.compartment.name + '.'

        if (
            len(self.attributes) > 0
            and self.attributes != []
            and show_attributes
        ):
            txt += '('
            for i in self.attributes:
                if i is not None:
                    txt += str(i) + ', '
            txt = txt[:-2] + ')'

        txt = txt.replace("'", '')
        txt += ']'

        return txt


class OrderedPolymerSpecies(OrderedComplexSpecies, OrderedPolymer):
    """Ordered polymer that can participate in chemical reactions.

    A polymer composed of Species (which are also OrderedMonomers) that can
    act as a reactant or product in chemical reactions. The internal
    species represent multiple binding sites and/or functional regions.

    Parameters
    ----------
    species : list of Species or list of [Species, direction]
        List of species monomers to form the polymer. Each element can be
        a Species or a [Species, direction] pair.
    name : str or None, optional
        Custom name for the polymer. If None, auto-generated from
        constituent species.
    material_type : str, default='ordered_polymer'
        Material type identifier for the polymer.
    compartment : Compartment, str, or None, optional
        Compartment containing the polymer.
    attributes : list of str or None, optional
        Attributes for the polymer species.
    circular : bool, default=False
        If True, the polymer has circular topology (e.g., plasmid).

    Attributes
    ----------
    polymer : tuple of Species
        Ordered tuple of species monomers in the polymer.
    species : tuple of Species
        Alias for `polymer` (inherited from OrderedPolymer).
    circular : bool
        Flag indicating circular topology.
    default_material : str
        Class attribute defining default material type.

    See Also
    --------
    OrderedPolymer : Base class for ordered polymers.
    OrderedComplexSpecies : Ordered complex base class.
    PolymerConformation : Set of polymers with connections.

    Notes
    -----
    When used as a reaction input, either the entire
    OrderedPolymerSpecies or one of its internal Species (with
    Species.parent = OrderedPolymerSpecies) can be passed to mechanisms.

    Species inside an `OrderedPolymerSpecies` model multiple binding
    sites and/or functional regions. `ComplexSpecies` can be formed at
    specific locations by passing the internal Species.

    The `circular` attribute indicates circular topology but does not
    automatically enforce circular constraints in operations.

    Examples
    --------
    Create a linear polymer:

    >>> S1 = bcp.Species('S1')
    >>> S2 = bcp.Species('S2')
    >>> polymer = bcp.OrderedPolymerSpecies(
    ...     species=[S1, S2],
    ...     name='my_polymer'
    ... )
    >>> len(polymer)
    2

    Create a circular polymer (plasmid):

    >>> plasmid = bcp.OrderedPolymerSpecies(
    ...     species=[S1, S2],
    ...     circular=True
    ... )
    >>> plasmid.circular
    True

    """

    default_material = 'ordered_polymer'

    def __init__(
        self,
        species,
        name=None,
        material_type=default_material,
        compartment=None,
        attributes=None,
        circular=False,
    ):
        self.material_type = material_type
        self.compartment = compartment
        self.parent = None
        self.position = None
        self.direction = None

        self.circular = circular

        if attributes is None:
            self.attributes = []
        else:
            self.attributes = attributes

        self._name = OrderedComplexSpecies._check_name(self, name)

        self.material_type = material_type
        # self.species = []
        monomers = []
        for specie in species:
            if isinstance(specie, OrderedPolymer) or isinstance(
                specie, PolymerConformation
            ):
                raise NotImplementedError(
                    "OrderedPolymer and PolymerConformation cannot be "
                    "used as a Monomer at this time."
                )
            elif isinstance(specie, Species) and isinstance(
                specie, OrderedMonomer
            ):
                monomers += [specie]
            elif (isinstance(specie, tuple) or isinstance(specie, list)) and (
                isinstance(specie[0], Species)
                and isinstance(specie[0], OrderedMonomer)
            ):
                monomers += [specie]
            else:
                raise ValueError(
                    f"{specie} should be a Species or list "
                    "[Species, 'direction']"
                )
                # only species are acceptable

        OrderedPolymer.__init__(self, monomers)
        self.material_type = material_type

    @classmethod
    def from_polymer_species(cls, ops, replace_dict, **kwargs):
        """Create OrderedPolymerSpecies with specific monomers replaced.

        Parameters
        ----------
        ops : OrderedPolymerSpecies
            The original polymer species to modify.
        replace_dict : dict
            Dictionary mapping monomer indices (int) to new Species to
            insert at those positions.
        **kwargs
            Additional keyword arguments for the new OrderedPolymerSpecies.
            Defaults are inherited from `ops` if not specified.

        Returns
        -------
        OrderedPolymerSpecies
            New polymer with specified monomers replaced.

        Notes
        -----
        If `replace_dict` is empty, returns a deep copy of `ops`.

        Examples
        --------
        Replace monomer at position 1:

        >>> S1 = bcp.Species('S1')
        >>> S2 = bcp.Species('S2')
        >>> S3 = bcp.Species('S3')
        >>> polymer = bcp.OrderedPolymerSpecies([S1, S2])
        >>> new_polymer = bcp.OrderedPolymerSpecies.from_polymer_species(
        ...     polymer, {1: S3}
        ... )

        """
        if replace_dict == {}:
            # nothing to replace!
            return copy.deepcopy(ops)
        monomers = []
        for i in range(len(ops.polymer)):
            if i in replace_dict:
                monomers.append(replace_dict[i])
            else:
                monomer = copy.copy(ops[i])
                direction = monomer.direction
                monomer.remove()
                monomers.append([monomer, direction])

        # Set keywords
        if 'circular' not in kwargs:
            kwargs['circular'] = ops.circular
        if 'material_type' not in kwargs:
            kwargs['material_type'] = ops.default_material
        if 'compartment' not in kwargs:
            kwargs['compartment'] = ops.compartment
        if 'attributes' not in kwargs:
            kwargs['attributes'] = ops.attributes
        # Produces a new OrderedPolymerSpecies
        return cls(monomers, **kwargs)

    @property
    def species_set(self):
        return set(self.polymer)

    @property
    def species(self):
        return self.polymer

    def get_species_list(self):
        return self.polymer

    @property
    def circular(self):
        if 'circular' in self.attributes:
            return True
        else:
            return False

    @circular.setter
    def circular(self, value: bool):
        if value:
            self.add_attribute('circular')
        else:
            self.remove_attribute('circular')

    def set_species_list(self, spec_tuple: tuple):
        OrderedPolymer.__init__(self, spec_tuple)

    @property
    def name(self):
        if self._name is None:
            name = ''
            for monomer in self.polymer:
                name += str(monomer) + '_'
            name = name[:-1]  # remove last underscore
        else:
            name = self._name

        return name

    def __hash__(self):
        ophash = OrderedPolymer.__hash__(self)
        ophash += hash(str(self))
        # hash(self.circular) + hash(self.name) + hash(self.material_type)
        # + hash(self.attributes)
        return ophash

    def replace(self, position, part, direction=None):
        # TODO only change the name if the part we are replacing is
        # actually different
        mydir = direction
        if (mydir is None) and (part.direction is not None):
            mydir = part.direction
        if (
            part == self.polymer[position]
            and self.polymer[position].direction == mydir
        ):
            # in this case we are replacing a part with the same thing, so
            # do nothing but it could be true that the reference changes?
            # That shouldnt be
            pass
        else:
            OrderedPolymer.replace(
                self, position=position, part=part, direction=mydir
            )

    def __contains__(self, item):
        for part in self.species:
            part = copy.copy(part)  # Only compare things which aren't None
            if item.parent is None:
                part.parent = None
            if item.direction is None:
                part.direciton = None
            if item.position is None:
                part.position = None
            if part == item:
                return True
        return False


class PolymerConformation(Species, MonomerCollection):
    """Set of polymers and their connections via ComplexSpecies.

    Represents a conformation of one or more PolymerSpecies connected by
    ComplexSpecies containing monomers from the polymers. This class
    provides unique naming for conformations and serves as a data structure
    for polymer hypergraphs.

    Parameters
    ----------
    complexes : list of ComplexSpecies, optional
        List of ComplexSpecies connecting monomers from
        OrderedPolymerSpecies. Must contain monomers from the polymers.
    polymer : OrderedPolymerSpecies or list of Species, optional
        Single polymer or list of species to form a polymer. Exactly one
        of `complexes` or `polymer` must be provided.
    material_type : str, default='conformation'
        Material type identifier.
    name : str or None, optional
        Custom name for the conformation. If None, auto-generated.
    **kwargs
        Additional keyword arguments passed to Species constructor.

    Attributes
    ----------
    polymers : list of OrderedPolymerSpecies
        List of polymers in this conformation.
    complexes : list of ComplexSpecies
        List of complexes connecting monomers in the polymers.
    name : str
        Auto-generated name encoding polymer and complex structure.

    See Also
    --------
    OrderedPolymerSpecies : Polymer species for chemical reactions.
    ComplexSpecies : Complex of multiple bound species.
    Complex : Metaclass for creating complexes.

    Notes
    -----
    Auto-generated names follow the format:
    `conformation__[Polymer1]_[Polymer2]_[indices]_[Complex1]_[Complex2]__`

    where indices encode which polymers each complex binds to and the list of
    `PolymerSpecies` and `ComplexSpecies` are in alphabetical order.

    A `PolymerConformation` represents a hypergraph where:

    - Monomers are vertices
    - `ComplexSpecies` are hyperedges connecting arbitrary numbers of
      vertices
    - Multiple edges between the same vertices are allowed

    Users typically do not create PolymerConformations directly. The
    `Complex` function automatically creates them when complexing
    monomers from OrderedPolymerSpecies.

    Examples
    --------
    Create from a single polymer:

    >>> S1 = bcp.Species('S1')
    >>> S2 = bcp.Species('S2')
    >>> polymer = bcp.OrderedPolymerSpecies([S1, S2])
    >>> conformation = bcp.PolymerConformation(polymer=polymer)

    """

    def __init__(
        self,
        complexes=None,
        polymer=None,
        material_type='conformation',
        name=None,
        **kwargs,
    ):
        Species.__init__(
            self, name=name, material_type=material_type, **kwargs
        )

        if isinstance(complexes, list) and len(complexes) == 0:
            complexes = None

        if (complexes is not None and polymer is not None) or (
            complexes is None and polymer is None
        ):
            raise ValueError(
                "PolymerConformation requires either: a list of "
                "ComplexSpecies which contain monomers from "
                "OrderedPolymerSpecies or a single OrderedPolymerSpecies "
                "in its constructor."
            )
        elif complexes is not None:
            self.complexes = complexes
        elif polymer is not None:
            # This order matters because complexes resets polymers
            self.complexes = []
            if isinstance(polymer, OrderedPolymerSpecies):
                self.polymers = [copy.deepcopy(polymer)]
            elif isinstance(polymer, list) and [
                isinstance(m, Species) for m in polymer
            ]:
                self.polymers = [OrderedPolymerSpecies(polymer)]
            else:
                raise ValueError(
                    f"polymer must be an OrderedPolymerSpecies or a list "
                    f"of Species. Received: {polymer}."
                )

    @classmethod
    def from_polymer_conformation(
        cls, pcs, complexes=None, complexes_to_remove=None, **kwargs
    ):
        """Create PolymerConformation from existing conformations.

        Produces a new PolymerConformation by merging complexes from
        previous PolymerConformations and adding new complexes.

        Parameters
        ----------
        pcs : list of PolymerConformation
            List of existing PolymerConformations to merge.
        complexes : list of ComplexSpecies, optional
            Additional complexes to add to the conformation. Default is an
            empty list.
        complexes_to_remove : list of ComplexSpecies, optional
            Complexes to exclude from the merged conformation. Default is
            an empty list.
        **kwargs
            Additional keyword arguments for the new PolymerConformation.

        Returns
        -------
        PolymerConformation
            New conformation merging all input conformations and complexes.

        Raises
        ------
        TypeError
            If `pcs` is not a list of PolymerConformations.

        """
        if not isinstance(pcs, list) or not any(
            [isinstance(pc, PolymerConformation) for pc in pcs]
        ):
            raise TypeError(
                f"pcs must be a list of PolymerConformations. Recieved {pcs}."
            )

        if complexes is None:
            complexes = []
        if complexes_to_remove is None:
            complexes_to_remove = []

        # generate a list of all complexes
        for pc in pcs:
            for c in pc.complexes:
                if not any([c == cc for cc in complexes_to_remove]):
                    complexes.append(c)

        return cls(complexes, **kwargs)

    @classmethod
    def from_polymer_replacement(
        cls, pc, old_polymers, new_polymers, **kwargs
    ):
        """Create PolymerConformation by replacing polymers.

        Produces a PolymerConformation from an existing one by replacing
        specified polymers with new ones, updating all complexes
        accordingly.

        Parameters
        ----------
        pc : PolymerConformation
            The conformation to modify.
        old_polymers : list of OrderedPolymerSpecies
            Polymers to replace. Must be the same instances (not just
            equal) as those in `pc.polymers`.
        new_polymers : list of OrderedPolymerSpecies
            New polymers to use as replacements. Must be the same length
            as `old_polymers`.
        **kwargs
            Additional keyword arguments for the new PolymerConformation.
            Defaults are inherited from `pc` if not specified.

        Returns
        -------
        PolymerConformation
            New conformation with polymers replaced.

        Raises
        ------
        TypeError
            If arguments are not the correct types.
        ValueError
            If `old_polymers` are not instances in `pc.polymers`, or if
            lists have different lengths.

        Notes
        -----
        This method updates all complexes to reference monomers from the
        new polymers at the same positions as in the old polymers.

        """
        if not isinstance(pc, PolymerConformation):
            raise TypeError(
                f"pc must be a PolymerConformation. Recieved {pc}."
            )

        if not isinstance(new_polymers, list) and all(
            [isinstance(p, OrderedPolymerSpecies) for p in new_polymers]
        ):
            raise TypeError(
                f"new_polymers must be a list of OrderedPolymerSpecies. "
                f"Recieved: {new_polymers}"
            )

        if not isinstance(old_polymers, list) and all(
            [isinstance(p, OrderedPolymerSpecies) for p in old_polymers]
        ):
            raise TypeError(
                f"new_polymers must be a list of OrderedPolymerSpecies. "
                f"Recieved: {old_polymers}"
            )

        if not all(
            [any([p is pp for pp in pc.polymers]) for p in old_polymers]
        ):
            raise ValueError(
                "All OrderedPolymerSpecies in old_polymers must be "
                "contained (as instances, not string equivalents) "
                "in pc.polymers."
            )

        if len(old_polymers) != len(new_polymers):
            raise ValueError(
                "old_polymers and new_polymers must be the same length."
            )

        complexes = pc.complexes
        new_complexes = []

        # Create a set of new ComplexeSpecies with the correct monomers inside
        for c in complexes:
            pol_inds = pc.get_polymer_indices(c)
            species = []
            for ci, s in enumerate(c.species):
                pi = pol_inds[ci]
                # If s is part of a polymer:
                if pi is not None:
                    # add the new_polymers[pi][s.position] to the new species
                    species.append(new_polymers[pi][s.position])
                # Otherwise, the species is not part of a polymer
                else:
                    species.append(s)

            # Create the new Complex
            if isinstance(
                c, OrderedComplexSpecies
            ):  # is the complex ordered?
                new_complexes += [
                    OrderedComplexSpecies(
                        species,
                        material_type=c.material_type,
                        attributes=c.attributes,
                        comparment=c.compartment,
                    )
                ]
            elif isinstance(c, ComplexSpecies):
                new_complexes += [
                    ComplexSpecies(
                        species,
                        material_type=c.material_type,
                        attributes=c.attributes,
                        compartment=c.compartment,
                    )
                ]
            else:
                raise TypeError(
                    f"Invalid object found in "
                    f"PolymerConformation.complexes: {c}."
                )

        # Set keywords
        if 'material_type' not in kwargs:
            kwargs['material_type'] = pc.material_type
        if 'compartment' not in kwargs:
            kwargs['compartment'] = pc.compartment
        if 'attributes' not in kwargs:
            kwargs['attributes'] = pc.attributes

        return cls(new_complexes, **kwargs)

    def copy_remove_complexes(self, complexes):
        """Returns a new PolymerConformation without these complexes."""
        if not isinstance(complexes, list):
            complexes = [complexes]

        # Check if the complexes are in the PolymerConformation
        for c in complexes:
            if c not in self.complexes:
                raise ValueError(
                    f"Complex {c} not in PolymerConformation {self}."
                )

        new_complexes = [c for c in self.complexes if c not in complexes]

        return PolymerConformation(complexes=new_complexes)

    # To be a valid MonomerCollection
    @property
    def monomers(self):
        return tuple(self._species)

    @property
    def polymers(self):
        return self._polymers

    @polymers.setter
    def polymers(self, polymers):
        assert all([isinstance(p, OrderedPolymerSpecies) for p in polymers])
        polymers.sort()
        for p in polymers:
            p.parent = self
        self._polymers = polymers

    @property
    def complexes(self):
        return self._complexes

    @complexes.setter
    def complexes(self, complexes):
        """Set complexes and connections.

        This setter copies the complexes and the polymers they connect into
        the PolymerConformation.  This is done in such a way to preserve
        references between parents and children without relying on hash
        functions (in case two polymers are identical).

        """
        self._complexes = []
        self._polymers = []
        if not isinstance(complexes, list) or not all(
            [isinstance(c, ComplexSpecies) for c in complexes]
        ):
            raise TypeError(
                f"complexes must be a list containing ComplexSpecies. "
                f"Recieved {complexes}."
            )

        complex_counts = [
            (
                sum(
                    [cc == c and cc.species == c.species for cc in complexes]
                ),
                c,
            )
            for c in complexes
        ]
        if any([i[0] > 1 for i in complex_counts]):
            raise ValueError(
                f"Complexes contains two or more identical complexes. "
                f"Recieved: {[c[1] for c in complex_counts if c[0] > 1]}."
            )

        # Sort the complexes by their name, and the ids of the polymers
        # inside of them to get a unique ordering for identically named
        # Complexes and Polymers this will produce a unique ordering of the
        # internal polymers as well
        def sort_func(c):
            return (c, tuple(id(s.parent) for s in c.species))

        complexes.sort(key=sort_func)

        polymers = []  # Original polymers stored here
        copied_polymers = []  # Polymers are copied here
        copied_complexes = []  # Complexes are copied here

        # Find all the polymers in the complexes. "is" is used instead of
        # equality to deal with the possibility of multiple instances
        # (copies) of the same polymer being bound together.
        for c in complexes:
            c_copy = copy.deepcopy(c)
            complex_contains_polymer = False
            for i, s in enumerate(c.species):
                parent = s.parent

                if parent is not None and isinstance(
                    parent, OrderedPolymerSpecies
                ):
                    complex_contains_polymer = True
                    # parent has not been copied
                    if not any([parent is p for p in polymers]):
                        polymers.append(parent)
                        parent_copy = copy.deepcopy(parent)
                        parent_copy.parent = self  # set the polymer's parent
                        copied_polymers.append(parent_copy)
                        # set s to be the new version from the deep copy
                        c_copy.species[i] = parent_copy[s.position]
                    # parent has already been copied
                    else:
                        index = [parent is p for p in polymers].index(
                            True
                        )  # get parent index
                        parent_copy = copied_polymers[index]
                        # set s to be the new version from the deep copy
                        c_copy.species[i] = parent_copy[s.position]

                # This edgecase should not occur.
                elif parent is not None and not isinstance(
                    parent, OrderedPolymerSpecies
                ):
                    raise ValueError(
                        f"Species {s} found inside complex {c} with a "
                        f"parent {parent} which is not an "
                        f"OrderedPolymerSpecies."
                    )

            if not complex_contains_polymer:
                raise ValueError(
                    f"Complex {c} does not contain any Species inside "
                    f"of Polymers."
                )

            # Set the parent of the ComplexSpecies to self
            c_copy.parent = self
            copied_complexes.append(c_copy)

        self.polymers = copied_polymers

        # Sort the complexes ordered by their polymer indices then name
        if len(copied_complexes) > 0:

            def complex_sort_func(c):
                return (
                    tuple(
                        [i if i else -1 for i in self.get_polymer_indices(c)]
                    ),
                    c,
                )

            copied_complexes.sort(key=complex_sort_func)
        self._complexes = copied_complexes

    def get_polymer_indices(self, c):
        # Takes a complex and returns a list of indices to the polymers
        # that complex contains monomers from this complex should be in
        # self.complexes.  If a Species in the Complex is not in a Polymer,
        # None is added to the list.
        indices = []
        for j, s in enumerate(c.species):
            parent = s.parent
            for i, p in enumerate(self.polymers):
                if parent is p:
                    indices.append(i)
            if len(indices) <= j:
                indices.append(None)

        return indices

    def get_polymer_positions(self, c, polymer_ind):
        # Takes a complex and the index of a polymer in the conformation
        # and returns a list of positions that ComplexSpecies is bound at
        p = self.polymers[polymer_ind]
        positions = []
        for s in c.species:
            if s.parent is p:
                positions.append(s.position)
            else:
                positions.append(None)
        return tuple(positions)

    def get_polymer(self, p):
        # Takes a polymer and returns a matching instance inside
        # self.polymers (or None)
        p = copy.copy(p)
        p.parent = self
        if p in self.polymers:
            i = self.polymers.index(p)
            return self.polymers[i]
        else:
            return None

    def count_polymer(self, p):
        p = copy.copy(p)
        p.parent = self
        count = 0
        for pp in self.polymers:
            if p == pp:
                count += 1
        return count

    def get_complex(self, c):
        # takes a ComplexSpecies and returns the instance of this Species
        # in the PolymerConformation (or None)
        c = copy.copy(c)
        c.parent = self
        # if c in self.complexes:
        #    i = self.complexes.index(c)
        #    return self.complexes[i]
        for cc in self.complexes:
            if cc == c and [
                (s.position, str(s.parent)) for s in c.species
            ] == [(s.position, str(s.parent)) for s in cc.species]:
                return cc
        return None

    def remove_complex(self, c):
        # Removes a copy of this complex from PolymerConformation if possible
        if c in self.complexes:
            self.complexes.remove(c)
        else:
            raise ValueError(
                f"Complex {c} not in PolymerConformation {self}."
            )

    def get_complexes_at(self, p_ind, p_pos):
        # returns all complexes in polymer p_ind bound to location p_pos
        found_complexes = []
        for c in self.complexes:
            for s in c.species:
                if s.parent == self.polymers[p_ind] and s.position == p_pos:
                    found_complexes.append(c)
        return found_complexes

    @property
    def name(self):
        if self._name is None:
            # If there is nothing in the PolymerConformation, use the name
            # of the internal OrderedPolymerSpecies
            if len(self.complexes) == 0 and len(self.polymers) == 1:
                return str(self.polymers[0])
            else:
                name = ''
                for p in self.polymers:
                    name += '_' + str(p)
                for c in self.complexes:
                    parent_inds = self.get_polymer_indices(c)
                    name += '_'
                    for i, ind in enumerate(parent_inds):
                        if ind is None:
                            name += 'n'
                        else:
                            locs = list(self.get_polymer_positions(c, ind))
                            loc = locs[i]

                            name += f"p{ind}l{loc}"
                    name += '_' + str(c)
                return name
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        self._name = self._check_name(name)

    def __repr__(self):
        """String representation of a polymer conformation.

        PolymerConformations add an additional "_" onto the end of their
        string representation This ensures that some edge cases are
        differentiated.

        A PolymerConformation with no attributes or complexes consisting of
        just one OrderedPolymerSpecies uses the same representation as the
        OrderedPolymerSpecies.

        """
        if len(self.complexes) == 0 and len(self.polymers) == 1:
            txt = Species.__repr__(self)
            txt = txt.replace('conformation_', '')
        else:
            txt = Species.__repr__(self)
            txt += '_'
        return txt


class Complex:
    """Metaclass for creating chemical complexes.

    `Complex` is not a class that gets instantiated directly - it creates
    instances of `ComplexSpecies`, `OrderedComplexSpecies`,
    `OrderedPolymerSpecies`, or `PolymerConformation` based on the input
    species and their parent relationships.

    Parameters
    ----------
    species : list of Species
        List of species to combine into a complex. Can include standalone
        Species, Species with parents (monomers in polymers), or entire
        OrderedPolymerSpecies.
    ordered : bool, default=False
        If True, creates OrderedComplexSpecies where species order
        matters. If False, creates ComplexSpecies where order is
        irrelevant.
    **kwargs
        Additional keyword arguments passed to the created species class.

    Returns
    -------
    ComplexSpecies, OrderedComplexSpecies, OrderedPolymerSpecies, or
    PolymerConformation
        The type of species returned depends on the input structure:

        - Simple species list -> ComplexSpecies or OrderedComplexSpecies
        - Monomers from one polymer -> OrderedPolymerSpecies
        - Monomers from multiple polymers/conformations ->
          PolymerConformation

    See Also
    --------
    ComplexSpecies : Unordered complex of multiple species.
    OrderedComplexSpecies : Ordered complex of multiple species.
    OrderedPolymerSpecies : Polymer species for reactions.
    PolymerConformation : Multiple polymers with connections.

    Notes
    -----
    The `__new__` method implements logic for different scenarios:

    1. No parents: Creates ComplexSpecies or OrderedComplexSpecies
    2. Single polymer parent: Creates OrderedPolymerSpecies with
       complex at binding site
    3. Multiple polymer parents or conformations: Creates
       PolymerConformation merging all complexes
    4. Error cases: Raises exceptions for invalid combinations

    The correct species type is automatically determined from the input,
    allowing flexible complex formation without explicit type selection.

    Examples
    --------
    Create a simple complex:

    >>> S1 = bcp.Species('S1')
    >>> S2 = bcp.Species('S2')
    >>> complex = bcp.Complex([S1, S2])
    >>> type(complex)
    biocrnpyler.core.species.ComplexSpecies

    Create an ordered complex:

    >>> ordered = bcp.Complex([S1, S2], ordered=True)
    >>> type(ordered)
    biocrnpyler.core.species.OrderedComplexSpecies

    Create a complex at a polymer binding site:

    >>> S3 = bcp.Species('S3')
    >>> polymer = bcp.OrderedPolymerSpecies([S1, S2])
    >>> # S1 is now inside the polymer at position 0
    >>> complex = bcp.Complex([polymer[0], S3])
    >>> type(complex.parent)
    biocrnpyler.core.species.OrderedPolymerSpecies

    """

    def __new__(cls, *args, **kwargs):
        """Create an instance of the appropriate species type.

        This method analyzes the input species and their parent
        relationships to determine which type of complex to create.

        Parameters
        ----------
        *args
            Positional arguments, first should be the species list.
        **kwargs
            Keyword arguments including 'species' and 'ordered'.

        Returns
        -------
        ComplexSpecies, OrderedComplexSpecies, OrderedPolymerSpecies, or
        PolymerConformation
            The appropriate species type based on input structure.

        Raises
        ------
        TypeError
            If species argument is not a list, or if trying to complex
            entire OrderedPolymerSpecies that are already in
            PolymerConformations, or if invalid parent types are found.
        ValueError
            If trying to form complexes between monomers from multiple
            OrderedPolymerSpecies without PolymerConformations.

        Notes
        -----
        Cases handled:

        1. No Species have parents -> `ComplexSpecies` or
           1OrderedComplexSpecies`
        2. Single Species has parent `OrderedPolymerSpecies` (no parent) ->
           `OrderedPolymerSpecies` with complex at binding site
        3. Multiple Species with OrderedPolymerSpecies1` parents (no
           parents) -> Error (must use PolymerConformations)
        4. Entire OrderedPolymerSpecies in PolymerConformations -> Error
        5. One or more `Species` from polymer Conformations ->
           `PolymerConformation` merging all complexes

        """
        species = []
        # below is extracting the "species" keywords from the args
        if 'species' in kwargs:
            species = kwargs.pop('species')
        elif len(args) >= 1:
            species = args[0]
            args = args[1:]

        if not isinstance(species, list):
            raise TypeError(
                f"First argument to Complex (or species keyword), "
                f"must be a list of Species; recieved {species}."
            )

        # Check whether ot make a ComplexSpecies or OrderedComplexSpecies
        if kwargs.pop('ordered', False):
            ComplexClass = OrderedComplexSpecies
        else:
            ComplexClass = ComplexSpecies

        # Use to supress errors in ComplexSpecies and OrderedComplexSpecies
        kwargs['called_from_complex'] = True

        # parent_species is a list of OrderedPolymerSpecies and/or
        # PolymerConformations
        parent_species = []
        # bindloc is location of a Species is bound to parent_species
        bindlocs = []

        # insertloc is the location of the species inside the
        # Complex.species list (important to maintain order in
        # OrderedComplexSpecies)
        insertlocs = []
        child_species = []  # the species with parents
        other_species = []  # Other species in the Complex

        # Below cycle through species and see if one has a parent. If it
        # does, that means the species is in an OrderedPolymerSpecies and
        # the Complex should be formed around it.
        for i, specie in enumerate(species):
            if hasattr(specie, 'parent') and (specie.parent is not None):
                # This adds to a growing list of parents, which will be
                # placed in an PolymerConformation. It is very important to
                # deepcopy here because the underlying
                # OrderedPolymerSpecies or PolymerConformation will be
                # modified.  insertlocs stores the order of the species for
                # creating OrderedComplexSpecies.
                parent_species.append(specie.parent)
                bindlocs.append(specie.position)
                insertlocs.append(i)
                child_species.append(specie)
            else:
                other_species += [specie]

        # Case 1: If no OrderedPolymerSpecies is found, just call the
        # regular constructor.
        if len(parent_species) == 0:
            return ComplexClass(species, *args, **kwargs)

        # Case 2: the Complex is being formed inside an
        # OrderedPolymerSpecies (which is not in a PolymerConformation).
        elif (
            len(parent_species) == 1
            and isinstance(parent_species[0], OrderedPolymerSpecies)
            and (parent_species[0].parent is None)
        ):
            parent_species = parent_species[0]
            bindloc = bindlocs[0]

            # Create an OrderedcomplexSepcies or ComplexSpecies
            child = copy.copy(child_species[0])
            child_direction = child_species[0].direction
            # Remove the child species, the new_complex will replace it
            child.remove()
            species.remove(child_species[0])
            # place the new child in the list in the proper location to
            # preserve order
            species.insert(insertlocs[0], child)
            new_complex = ComplexClass(
                species, *args, **kwargs
            )  # create the Complex

            # Create a new OrderedPolymerSpecied which is copied from the
            # parent with the new complex replacing bindloc (inheriting the
            # same direction).
            new_polymer_species = OrderedPolymerSpecies.from_polymer_species(
                parent_species, {bindloc: [new_complex, child_direction]}
            )

            # Case 2: OrderedPolymerSpecies has no parent
            if parent_species.parent is None:
                return new_polymer_species[bindloc]

            # # (Old Case 3 now case 4) OrderedPolymerSpecies is inside a
            # # PolymerConformation has been moved to case
            # elif isinstance(parent_species.parent, PolymerConformation):
            #    raise RuntimeError(
            #        "This case should no longer occur as defined")
            # # create a new PolymerConformation that replaces the
            # # appropriate monomer in the polymer
            #
            # new_pc = PolymerConformation.from_polymer_replacement(
            #    parent_species.parent, [parent_species],
            #    [new_polymer_species])
            #
            # # return the newly created Monomer attached to the new
            # # OrderedPolymerSpecies inside the new PolymerConformation
            # return new_pc.get_polymer(
            #     new_polymer_species)[bindloc]
            # else:
            #    #This error should never occur
            #    raise TypeError(
            #        f"Unknown parent type {type(parent_species)} recieved "
            #        f"for {parent_species} .parent {parent_species.parent}.")

        # Case 3-4: In the following cases, multiple species may have parents
        elif all(
            [isinstance(p, OrderedPolymerSpecies) for p in parent_species]
        ) and any([p.parent is None for p in parent_species]):
            raise TypeError(
                "In order to form Complexes between Monomers inside "
                "OrderedPolymerSpecies, each OrderedPolymerSpecies must be "
                "placed in one or more PolymerConformations."
            )

        elif any(
            [
                isinstance(s, OrderedPolymerSpecies) and s.parent is not None
                for s in species
            ]
        ):
            raise ValueError(
                "OrderedPolymerSpecies cannot be combined into a Complex "
                "if they are already part of an PolymerConformation. "
                "Maybe you meant to Complex specific monomers?"
            )

        # Case 5: Multiple species in one more more PolymerConformations
        # are being Complexed Together
        else:
            pcs = []
            # these Species will go inside the ComplexSpecies later
            merged_species = other_species

            # Cycle through Species with parents
            insert_loc_offset = 0
            for i, p in enumerate(parent_species):
                s = child_species[i]

                # if the parent is a OrderedPolymerSpecies
                if isinstance(p, OrderedPolymerSpecies):
                    merged_species.insert(
                        insertlocs[i] + insert_loc_offset, s
                    )

                    # if the Polymer is already in a Conformation...
                    if p.parent is not None and not any(
                        [p.parent is pp for pp in pcs]
                    ):
                        pcs.append(p.parent)

                # if the parent is a PolymerConformation and child is a
                # ComplexSpecies
                elif (
                    isinstance(p, PolymerConformation)
                    and isinstance(s, ComplexSpecies)
                    and not isinstance(s, OrderedPolymerSpecies)
                ):
                    # Store all the unique PolymerConformations
                    if not any([p is pp for pp in pcs]):
                        p = copy.copy(p)
                        p.remove_complex(s)
                        pcs.append(p)

                    # Merge the species lists
                    merged_species = (
                        merged_species[: insertlocs[i] + insert_loc_offset]
                        + s.species
                        + merged_species[insertlocs[i] + insert_loc_offset :]
                    )
                    # this takes care of ordering offsets during the merge
                    insert_loc_offset += len(s.species) - 1
                else:
                    raise TypeError(
                        f"Cannot form a complex from {species}. "
                        f"Invalid Parent Species {p} for child {s}."
                    )

            # Create a Complex and merged PolymerConformation
            new_complex = ComplexClass(merged_species, *args, **kwargs)
            new_pc = PolymerConformation.from_polymer_conformation(
                pcs, [new_complex]
            )
            return new_pc.get_complex(new_complex)
