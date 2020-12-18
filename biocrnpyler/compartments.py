
import warnings

class Compartment():

    """ A formal Compartment object for a Species in a CRN.
     A Compartment must have a name. They may also have a spatial dimension (such as 2 for two-dimensional,
     or 3 for three-dimensional) and the volume in litres. 
    """

    def __init__(self, name: str, spatial_dimensions = 3, volume=1e-6, **keywords):
        self.name = name
        self.spatial_dimensions = spatial_dimensions 
        self.volume = volume 

    @property
    def name(self):
        if self._name is None:
            return ""
        else:
            return self._name

    @name.setter
    def name(self, name: str):
        if name is None:
            raise TypeError("Name must be a string.")
        else:
            no_underscore_string = name.replace("_", "")
            if no_underscore_string.isalnum() and "__" not in name and name[len(name)-1] != "_" and not name[0].isnumeric():
                self._name = name
            else:
                raise ValueError(f"name attribute {name} must consist of letters, numbers, or underscores and cannot contained double underscores or end in an underscore.")

    @property
    def spatial_dimensions(self):
        return self._spatial_dimensions

    @spatial_dimensions.setter
    def spatial_dimensions(self, spatial_dimensions: int):
        self._spatial_dimensions = spatial_dimensions

    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, volume: float):
        self._volume = volume
        