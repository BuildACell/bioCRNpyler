# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import warnings


class Compartment:

    """ A formal Compartment object for a Species in a CRN.
     A Compartment must have a name. They may also have a spatial dimension (such as 2 for two-dimensional,
     or 3 for three-dimensional) and the size in litres. 
     Note: The "default" keyword is reserved for BioCRNpyler allotting a default compartment. 
     Users must choose a different string. 
     The unit attribute for Compartment can be used to set unit for the Compartment size. Make sure
     that the string identifier used for the unit is a supported unit in BioCRNpyler. Check the 
     documentation to find a list of supported units. Add your own custom units in units.py if needed.
    """

    def __init__(self, name: str, size=1e-6, spatial_dimensions=3, unit=None, **keywords):
        self.name = name
        self.spatial_dimensions = spatial_dimensions
        self.size = size
        self.unit = unit

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name: str):
        if name is None:
            raise TypeError("Compartment name must be a string.")
        elif type(name) is str:
            no_underscore_string = name.replace("_", "")
            if no_underscore_string.isalnum() and "__" not in name and name[len(name)-1] != "_" and name[0].isalpha():
                self._name = name
            else:
                raise ValueError(
                    f"name attribute {name} must consist of letters, numbers, or underscores and cannot contained double underscores or begin/end with a special character.")
        else:
            raise ValueError('Compartment name must be a string.')

    @property
    def spatial_dimensions(self):
        return self._spatial_dimensions

    @spatial_dimensions.setter
    def spatial_dimensions(self, spatial_dimensions: int):
        if type(spatial_dimensions) is not int:
            raise ValueError(
                'Compartment spatial dimension must be an integer.')
        elif spatial_dimensions < 0:
            raise ValueError(
                'Compartment spatial dimension must be non-negative.')
        else:
            self._spatial_dimensions = spatial_dimensions

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, size: float):
        if type(size) not in [float, int]:
            raise ValueError('Compartment size must be a float or int.')
        elif size < 0:
            raise ValueError('Compartment size must be non-negative.')
        else:
            self._size = size

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self, unit: str):
        if unit is not None:
            if not isinstance(unit, str):
                raise ValueError(
                    "Unit of compartment must be a string representing compartment size.")
            self._unit = unit
        else:
            self._unit = None

    def __eq__(self, other):
        """
        Overrides the default implementation
        Two compartments are equivalent if they have the same name, spatial dimension, and size
        :param other: Compartment instance
        :return: boolean
        """
        if isinstance(other, Compartment) and self.name == other.name:
            # Now if two compartments have same name but other attributes are different, throw an error:
            if self.size != other.size or self.spatial_dimensions != other.spatial_dimensions:
                raise ValueError(
                    'Compartments with same names must have the same size and spatial dimensions.')
            return True
        else:
            return False
