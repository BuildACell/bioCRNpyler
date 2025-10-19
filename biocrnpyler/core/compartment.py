# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


class Compartment:
    """Spatial compartment for organizing species in a CRN model.

    Compartments represent physically distinct regions where chemical species
    can exist, such as the cytoplasm, nucleus, extracellular space, or
    organelles. Each compartment has a name, size, and spatial dimensionality.
    Species in different compartments are treated as distinct, even if they
    have the same molecular identity.

    Parameters
    ----------
    name : str
        Name of the compartment. Must consist of letters, numbers, or
        underscores. Cannot contain double underscores, and cannot begin or
        end with special characters. Must start with a letter. The name
        'default' is reserved by BioCRNpyler.
    size : float or int, default=1e-6
        Size of the compartment in the units specified by `unit`. Default is
        1 microliter (1e-6 liters).
    spatial_dimensions : int, default=3
        Number of spatial dimensions (0 for point, 1 for line, 2 for surface,
        3 for volume). Must be non-negative.
    unit : str, optional
        Unit identifier for the compartment size (e.g., 'L', 'mL', 'µL').
        Must be a supported unit in BioCRNpyler. See documentation for
        supported units or add custom units in 'core/units.py'.

    Attributes
    ----------
    name : str
        Name of the compartment.
    size : float
        Size of the compartment.
    spatial_dimensions : int
        Number of spatial dimensions.
    unit : str or None
        Unit identifier for the compartment size.

    Raises
    ------
    TypeError
        If `name` is None.
    ValueError
        If `name` is not a string, contains invalid characters, or if `size`
        or `spatial_dimensions` are invalid.

    See Also
    --------
    Species : Chemical species that can be assigned to compartments.
    Mixture : Container that can have a default compartment.

    Notes
    -----
    The reserved name 'default' is used internally by BioCRNpyler for species
    that have not been explicitly assigned to a compartment. User-defined
    compartments should use other names.

    Two compartments are considered equal if they have the same name. If two
    compartments have the same name but different sizes or spatial dimensions,
    a ValueError is raised to prevent inconsistencies.

    Examples
    --------
    Create a cytoplasm compartment:

    >>> cytoplasm = bcp.Compartment(
    ...     name="cytoplasm",
    ...     size=1e-15,  # 1 femtoliter (bacterial cell volume)
    ...     spatial_dimensions=3,
    ...     unit="L"
    ... )

    Create a membrane compartment (2D):

    >>> membrane = bcp.Compartment(
    ...     name="membrane",
    ...     size=1e-12,  # 1 square micrometer
    ...     spatial_dimensions=2,
    ...     unit="m^2"
    ... )

    Use compartments with species:

    >>> species_cyto = bcp.Species("Protein_X", compartment=cytoplasm)
    >>> species_mem = bcp.Species("Protein_X", compartment=membrane)
    >>> species_cyto == species_mem  # False - different compartments

    """

    def __init__(self, name: str, size=1e-6, spatial_dimensions=3, unit=None):
        self.name = name
        self.spatial_dimensions = spatial_dimensions
        self.size = size
        self.unit = unit

    @property
    def name(self):
        """str: Name of the compartment."""
        return self._name

    @name.setter
    def name(self, name: str):
        """Set the compartment name with validation.

        Parameters
        ----------
        name : str
            Name for the compartment. Must consist of letters, numbers, or
            underscores. Cannot contain double underscores ('__'), and cannot
            begin or end with underscores. Must start with a letter.

        Raises
        ------
        TypeError
            If `name` is None.
        ValueError
            If `name` is not a string or contains invalid characters.

        """
        if name is None:
            raise TypeError("Compartment name must be a string.")
        elif isinstance(name, str):
            no_underscore_string = name.replace('_', '')
            if (
                no_underscore_string.isalnum()
                and '__' not in name
                and name[len(name) - 1] != '_'
                and name[0].isalpha()
            ):
                self._name = name
            else:
                raise ValueError(
                    f"name attribute {name} must consist of letters, "
                    f"numbers, or underscores and cannot contained double "
                    "underscores or begin/end with a special character."
                )
        else:
            raise ValueError("Compartment name must be a string.")

    @property
    def spatial_dimensions(self):
        """int: Number of spatial dimensions.

        0 for point, 1 for line, 2 for surface, 3 for volume.

        """
        return self._spatial_dimensions

    @spatial_dimensions.setter
    def spatial_dimensions(self, spatial_dimensions: int):
        """Set the spatial dimensions with validation.

        Parameters
        ----------
        spatial_dimensions : int
            Number of spatial dimensions. Must be a non-negative integer.
            Common values: 0 (point), 1 (line), 2 (surface), 3 (volume).

        Raises
        ------
        ValueError
            If `spatial_dimensions` is not an integer or is negative.

        """
        if not isinstance(spatial_dimensions, int):
            raise ValueError(
                "Compartment spatial dimension must be an integer."
            )
        elif spatial_dimensions < 0:
            raise ValueError(
                "Compartment spatial dimension must be non-negative."
            )
        else:
            self._spatial_dimensions = spatial_dimensions

    @property
    def size(self):
        """float: Size of compartment in units specified by unit attribute."""
        return self._size

    @size.setter
    def size(self, size: float):
        """Set the compartment size with validation.

        Parameters
        ----------
        size : float or int
            Size of the compartment. Must be non-negative. Units are specified
            by the `unit` attribute.

        Raises
        ------
        ValueError
            If `size` is not a float or int, or is negative.

        """
        if not isinstance(size, (float, int)):
            raise ValueError("Compartment size must be a float or int.")
        elif size < 0:
            raise ValueError("Compartment size must be non-negative.")
        else:
            self._size = size

    @property
    def unit(self):
        """str: Unit identifier for compartment size (e.g., 'mL', 'uL')."""
        return self._unit

    @unit.setter
    def unit(self, unit: str):
        """Set the unit identifier with validation.

        Parameters
        ----------
        unit : str or None
            Unit identifier for the compartment size. Must be a supported unit
            in BioCRNpyler (e.g., 'L', 'mL', 'µL' for volumes). See
            documentation for supported units. Can be None for unitless sizes.

        Raises
        ------
        ValueError
            If `unit` is not a string (when not None).

        """
        if unit is not None:
            if not isinstance(unit, str):
                raise ValueError(
                    "Unit of compartment must be a string representing "
                    "compartment size."
                )
            self._unit = unit
        else:
            self._unit = None

    def __eq__(self, other):
        """Check equality of compartments by name.

        Two compartments are considered equal if they have the same name.
        If two compartments have the same name but different sizes or spatial
        dimensions, a ValueError is raised to prevent inconsistencies.

        Parameters
        ----------
        other : Compartment
            Another compartment to compare with.

        Returns
        -------
        bool
            True if compartments have the same name (and consistent
            attributes), False otherwise.

        Raises
        ------
        ValueError
            If compartments have the same name but different sizes or spatial
            dimensions.

        Notes
        -----
        This comparison is based solely on the compartment name. If two
        compartments share a name, they must also share the same physical
        properties (size and spatial dimensions) to maintain model
        consistency.

        """
        if isinstance(other, Compartment) and self.name == other.name:
            # Now if two compartments have same name but other attributes
            # are different, throw an error:
            if (
                self.size != other.size
                or self.spatial_dimensions != other.spatial_dimensions
            ):
                raise ValueError(
                    "Compartments with same names must have the same size "
                    "and spatial dimensions."
                )
            return True
        else:
            return False
