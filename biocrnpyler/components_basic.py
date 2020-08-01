from .component import *
from .species import Species, ComplexSpecies, OrderedComplexSpecies
from .mechanism import *
from .mechanisms_binding import *


class DNA(Component):
    """The DNA class is used to represent a DNA sequence that has a given length.

    Produces no reactions.
    """

    def __init__(self, name: str, length=0, attributes=None, **keywords):
        """Initialize a DNA object to store DNA related information.

        :param name: Name of the sequence (str)
        :param length: length of the basepairs (int)
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        self.species = self.set_species(name, material_type="dna", attributes=attributes)
        self.length = length
        Component.__init__(self=self, name=name, **keywords)

    def get_species(self) -> Species:
        return self.species

    def update_species(self) -> List[Species]:
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        return []

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, dna_length):
        if dna_length is None:
            self._length = dna_length
        elif dna_length >= 0 and isinstance(dna_length, int):
            self._length = dna_length
        else:
            raise ValueError("Length cannot be negative!")


class RNA(Component):
    """A class to represent Components made of RNA. Produces no reactions."""

    def __init__(self, name: str, length=0, attributes=None, **keywords):
        """Initialize a RNA object to store RNA related information

        :param name: name of the rna
        :param length: number of basepairs (int)
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """

        self.length = length
        self.species = self.set_species(name, material_type="rna",
                                        attributes=attributes)
        Component.__init__(self=self, name=name, **keywords)

    def get_species(self) -> Species:
        return self.species

    def update_species(self) -> List[Species]:
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        return []


class Protein(Component):
    """A class to represent Components made of Protein. Produces no reactions."""

    def __init__(self, name: str, length=0, attributes=None, **keywords):
        """Initialize a Protein object to store Protein related information.

        :param name: name of the protein
        :param length: length of the protein in number of amino acids
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """

        self.length = length
        self.species = self.set_species(name, material_type="protein",
                                        attributes=attributes)

        Component.__init__(self=self, name=name, **keywords)

    def get_species(self) -> Species:
        return self.species

    def update_species(self) -> List[Species]:
        species = [self.get_species()]
        return species

    def update_reactions(self) -> List:
        return []


class ChemicalComplex(Component):
    """A complex forms when two or more species bind together Complexes inherit the attributes of their species."""

    def __init__(self, species: List[Species], name: str=None, material_type="complex",
                 attributes=None, **keywords):
        """Initialize a ChemicalComplex object to store ChemicalComplex related information.

        :param species: list of species inside a complex
        :param name: name of the complex
        :param material_type: option to rename the material_type, default: complex
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        if not isinstance(species, list) or len(species) < 2:
            raise ValueError(f"Invalid Species {species}. Species must be a list of Species, strings, Component objects.")

        self.internal_species = []  # a list of species inside the complex

        for s in species:
            self.internal_species.append(self.set_species(s))
        if attributes is None:
            attributes = []
        self.species = Complex(species=self.internal_species, name=name, material_type=material_type, attributes=attributes)
        
        if name is None:
            name = self.species.name

        Component.__init__(self=self, name=name, **keywords)

    def get_species(self) -> List[Species]:
        return self.species

    def update_species(self) -> List[Species]:

        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        species = mech_b.update_species(binder, bindee, complex_species=self.get_species(), component=self, part_id=self.name)

        return species

    def update_reactions(self) -> List[Reaction]:

        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        reactions = mech_b.update_reactions(binder, bindee, complex_species=self.get_species(), component=self, part_id=self.name)
        
        return reactions


class Enzyme(Component):
    """A class to represent Enzymes.

    Assumes the enzyme converts a single substrate to a single product.
    Uses a mechanism called "catalysis"
    """
    def __init__(self, enzyme: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 product: Union[Species,str, Component], attributes=None, **keywords):
        """Initialize an Enzyme object to store Enzyme related information.

        :param enzyme: name of the enzyme, reference to an Species or Component
        :param substrate: name of the enzyme, reference to an Species or Component
        :param product: name of the product, reference to an Species or Component
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
      
        # ENZYME NAME
        self.enzyme = self.set_species(enzyme, material_type='protein', attributes=attributes)
    
        # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            self.substrate = self.set_species(substrate)

        # PRODUCT
        if product is None:
            self.product = None
        else:
            self.product = self.set_species(product)

        Component.__init__(self=self, name=self.enzyme.name, **keywords)
    
    def get_species(self):
        return self.enzyme

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(self.enzyme, self.substrate, self.product) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_reactions(self.enzyme, self.substrate, self.product, component=self,  part_id=self.name)


class MultiEnzyme(Component):
    """A class to represent Enzymes with multiple substrates and products.

    Assumes the enzyme converts all substrates to a all products at once.
    For example: S1 + S2 + E --> P1 + P2 + E.
    For enzymes with multiple enzymatic reactions, create multiple Enzyme Components with the same internal species.
    Uses a mechanism called "catalysis"
    """
    def __init__(self, enzyme: Union[Species, str, Component],
                 substrates: List[Union[Species, str, Component]],
                 products: List[Union[Species, str, Component]], attributes=None, **keywords):
        """Initialize an MultiEnzyme object to store MultiEnzyme related information.

        :param enzyme: name of the enzyme, reference to an Species or Component
        :param substrate: list of (name of the enzyme, reference to an Species or Component)
        :param product: list of (name of the product, reference to an Species or Component)
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
      
        # ENZYME NAME
        self.enzyme = self.set_species(enzyme, material_type='protein', attributes=attributes)

        # SUBSTRATE(s)
        self.substrates = []
        for substrate in substrates:
            self.substrates.append(self.set_species(substrate))

        # PRODUCT(s)
        self.products = []
        for product in products:
            self.products.append(self.set_species(product))
      
        Component.__init__(self=self, name=self.enzyme.name, **keywords)
    
    def get_species(self):
        return self.enzyme

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(self.enzyme, self.substrates, self.products)

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_reactions(self.enzyme, self.substrates, self.products, component=self,  part_id=self.name)
