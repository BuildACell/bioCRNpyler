from .component import *
from .species import Species, ComplexSpecies, OrderedComplexSpecies
from .mechanism import *
from .mechanisms_binding import *

# These subclasses of Component represent different kinds of biomolecules.
class DNA(Component):
    """DNA class

    The DNA class is used to represent a DNA sequence that has a given
    length.  Its main purpose is as the parent object for DNA
    fragments and DNA assemblies.

    Note: for initialization of members of this class, the arguments
    should be as follows:

      DNA(name, length, [mechanisms], [config_file], [prefix])

        DNAtype(name, required_arguments, [length], [mechanisms],
                [config_file], [prefix], [optional_arguments])

          DNAelement(name, required_arguments, [length], [mechanisms],
                     [config_file], [optional_arguments])


    Data attributes
    ---------------
    name        Name of the sequence (str)
    length      Length of the sequence (int)
    mechanisms  Local mechanisms for this component (overrides defaults)
    parameters  Parameter dictionary for the DNA element

    """

    def __init__(
            self, name: str, length=0, attributes = None,  # positional arguments
            **keywords
    ):
        self.species = self.set_species(name, material_type="dna",
                               attributes=attributes)
        self._length = length
        Component.__init__(self=self, name=name, **keywords)

    def get_species(self):
        return self.species

    def update_species(self):
        species = [self.get_species()]
        return species

    def update_reactions(self):
        return []

    @property
    def length(self):
        return  self._length

    @length.setter
    def length(self, dna_length):
        if dna_length >= 0:
            self._length = dna_length
        else:
            raise ValueError("Length cannot be negative!")


class RNA(Component):
    """
    A class to represent Components made of RNA. Produces no reactions.
    """
    def __init__(
            self, name: str, length=0, attributes = None,  # positional arguments,
            **keywords
    ):
        self.length = length
        self.species = self.set_species(name, material_type="rna",
                               attributes=attributes)
        Component.__init__(self=self, name=name, **keywords)

    def get_species(self):
        return self.species

    def update_species(self):
        species = [self.get_species()]
        return species

    def update_reactions(self):
        return []


class Protein(Component):
    """
    A class to represent Components made of Protein. Produces no reactions.
    """
    def __init__(
            self, name: str, length=0, attributes = None,  # positional arguments
            **keywords
    ):
        self.length = length
        self.species = self.set_species(name, material_type="protein",
                               attributes=attributes)

        Component.__init__(self=self, name=name, **keywords)

    def get_species(self):
        return self.species

    def update_species(self):
        species = [self.get_species()]
        return species

    def update_reactions(self):
        return []


class ChemicalComplex(Component):
    """
    A complex forms when two or more species bind together
    Complexes inherit the attributes of their species
    """
    def __init__(
            self, species: List[Species],  # positional arguments
            name=None,  # Override the default naming convention for a complex
            material_type="complex",
            attributes = None,
            **keywords
    ):

        if not isinstance(species, list) or len(species) < 2:
            raise ValueError(f"Invalid Species {species}. Species must be a list of Species, strings, Component objects.")

        self.internal_species = [] #a list of species inside the complex

        for s in species:
            self.internal_species.append(self.set_species(s))
        if(attributes is None):
            attributes = []
        self.species = Complex(species = self.internal_species, name = name, material_type=material_type, attributes=attributes)
        print(self.species)
        
        if name is None:
            name = self.species.name

        Component.__init__(self=self, name=name, **keywords)

    def get_species(self):
        return self.species

    def update_species(self) -> List[Species]:

        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        species = mech_b.update_species(binder,bindee, complex_species = self.get_species(), component = self, part_id = self.name)

        return species

    def update_reactions(self) -> List[Reaction]:

        mech_b = self.get_mechanism('binding')
        bindee = self.internal_species[0]
        binder = self.internal_species[1:]
        reactions = mech_b.update_reactions(binder,bindee, complex_species = self.get_species(), component = self, part_id = self.name)
        
        return reactions

class Enzyme(Component):
    """
    A class to represent Enzymes.
    Assumes the enzyme converts a single substrate to a single product.
    Uses a mechanism called "catalysis"
    """
    def __init__(self, enzyme, substrate, product, attributes = None, **keywords):
      
        # ENZYME NAME
        self.enzyme = self.set_species(enzyme, material_type = 'protein', attributes = attributes)
    
        # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            self.substrate = self.set_species(substrate)
        if product is None:
            self.product = None
        else:
            self.product = self.set_species(product)

      
        Component.__init__(self = self, name = self.enzyme.name, **keywords)
    
    def get_species(self):
        return self.enzyme

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(self.enzyme, self.substrate, self.product) 
                                                                                           
    
    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_reactions(self.enzyme, self.substrate, self.product, component = self,  part_id = self.name)


class MultiEnzyme(Component):
    """
    A class to represent Enzymes with multiple substrates and products.
    Assumes the enzyme converts all substrates to a all products at once.
    For example: S1 + S2 + E --> P1 + P2 + E.
    For enzymes with multiple enzymatic reactions, create multiple Enzyme Components with the same internal species.
    Uses a mechanism called "catalysis"
    """
    def __init__(self, enzyme, substrates, products, attributes = None, **keywords):
      
        # ENZYME NAME
        self.enzyme = self.set_species(enzyme, material_type = 'protein', attributes = attributes)
    
        # SUBSTRATE
        self.substrates = []
        for substrate in substrates:
            self.substrates.append(self.set_species(substrate))

        self.products = []
        for product in products:
            self.products.append(self.set_species(product))
      
        Component.__init__(self = self, name = self.enzyme.name, **keywords)
    
    def get_species(self):
        return self.enzyme

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(self.enzyme, self.substrates, self.products) 
                                                                                           
    
    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_reactions(self.enzyme, self.substrates, self.products, component = self,  part_id = self.name)
