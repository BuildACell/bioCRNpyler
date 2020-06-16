from .component import *
from .chemical_reaction_network import Species, ComplexSpecies, OrderedComplexSpecies
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
    assy        DNA assembly that we are part of
    mechanisms  Local mechanisms for this component (overrides defaults)
    parameters  Parameter dictionary for the DNA element

    """

    def __init__(
            self, name: str, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            initial_conc=None,
            parameter_warnings = True,
            **keywords
    ):
        self.species = Species(name, material_type="dna",
                               attributes=list(attributes))
        self._length = length
        Component.__init__(self=self, name=name, mechanisms=mechanisms,
                           parameters=parameters, attributes=attributes,
                           initial_conc=initial_conc,
                           parameter_warnings = parameter_warnings, **keywords)

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
    def __init__(
            self, name: str, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            initial_conc=None,
            **keywords
    ):
        self.length = length
        self.species = Species(name, material_type="rna",
                               attributes=list(attributes))
        Component.__init__(self=self, name=name, mechanisms=mechanisms,
                           parameters=parameters, attributes=attributes,
                           initial_conc=initial_conc, **keywords)

    def get_species(self):
        return self.species

    def update_species(self):
        species = [self.get_species()]
        return species

    def update_reactions(self):
        return []


class Protein(Component):
    def __init__(
            self, name: str, length=0,  # positional arguments
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters
            attributes=[],
            initial_conc=None,
            **keywords
    ):
        self.length = length
        self.species = Species(name, material_type="protein",
                               attributes=attributes)

        Component.__init__(self=self, name=name, mechanisms=mechanisms,
                           parameters=parameters, attributes=attributes,
                           initial_conc=initial_conc, **keywords)

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
            self, species,  # positional arguments
            name = None, #Override the default naming convention for a complex
            mechanisms={},  # custom mechanisms
            parameters={},  # customized parameters,
            attributes=[],
            initial_conc=None,
            material_type = "complex",
            **keywords
    ):

        if len(species) < 2 or not isinstance(species, list):
            raise ValueError("Species must be a list of Species, strings, Component objects.")

        self.internal_species = [] #a list of species inside the complex

        for s in species:
            self.internal_species.append(self.set_species(s))

        self.species = ComplexSpecies(species = self.internal_species, name = name, material_type=material_type, attributes=list(attributes))

        if name is None:
            name = self.species.name

        self.default_mechanisms = {"binding": One_Step_Binding()}

        Component.__init__(self=self, name=name, mechanisms=mechanisms,
                           parameters=parameters, attributes=attributes,
                           initial_conc=initial_conc, **keywords)

    def get_species(self):
        return self.species

    def update_species(self) -> List[Species]:

        mech_b = self.mechanisms['binding']

        species = mech_b.update_species(self.internal_species, complex_species = self.get_species(), component = self, part_id = self.name)

        return species

    def update_reactions(self) -> List[Reaction]:

        mech_b = self.mechanisms['binding']

        reactions = mech_b.update_reactions(self.internal_species, complex_species = self.get_species(), component = self, part_id = self.name)
        
        return reactions


class Enzyme(Component):
    def __init__(self, enzyme, substrate, product, **keywords):
      
        # ENZYME NAME
        self.enzyme = self.set_species(enzyme, material_type = 'protein')
    
        # SUBSTRATE
        self.substrate = self.set_species(substrate)
        self.product = self.set_species(product)
        self.substrate_list = []
      
        Component.__init__(self = self, name = self.enzyme.name, **keywords)
        
    def update_species(self):
        mech_cat = self.mechanisms['catalysis']

        return mech_cat.update_species(self.enzyme, self.fuel_list, self.substrate_list, self.product_list, self.waste_list) 
                                                                                           
    
    def update_reactions(self):
        mech_cat = self.mechanisms['catalysis']

        return mech_cat.update_reactions(self.enzyme, self.fuel_list, self.substrate_list, self.product_list, self.waste_list, component = None,  part_id = None)
    