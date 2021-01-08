
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from .component import Component
from .reaction import Reaction
from .species import Complex, Species


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

class Metabolite(Component):
    def __init__(self, name: str, attributes=None, precursors = None, products = None, **keywords):
        """Initialize a Metabolite object to store Metabolite related information.
        Metabolites look for "production" and "degredation" mechanisms, but will not throw
        an error if none are found.
            
        :param name: name of the protein
        :param attributes: Species attribute
        :param precursors: list of chemical species which are directly transformed into this metabolite via the production mechanism
        :param products: list of chemical species directly produced from this metabolite via the degredation mechanism
        :param keywords: pass into the parent's (Component) initializer
        """

        self.species = self.set_species(name, material_type="metabolite",
                                        attributes=attributes)


        #Set percursor species list
        self.precursors = []
        if precursors is not None:
            for p in precursors:
                if p is None:
                    self.precursors.append(None) #None is a valid precursor representing constuitive production
                else:
                    self.precursors.append(self.set_species(p))

        #Set product species list
        self.products = []
        if products is not None:
            for p in products:
                if p is None:
                    self.products.append(None) #None is a valid product representing total degredation
                else:
                    self.products.append(self.set_species(p))


        Component.__init__(self=self, name=name, **keywords)

    def get_species(self) -> Species:
        return self.species

    def update_species(self) -> List[Species]:
        species = [self.get_species()]
        mech_pathway = self.get_mechanism('metabolic_pathway', optional_mechanism = True)

        if mech_pathway is not None:
            if len(self.precursors) > 0:
                species += mech_pathway.update_species(precursor = self.precursors, product = [self.get_species()], component = self, part_id = self.name+"_production")
            if len(self.products) > 0:
                species += mech_pathway.update_species(precursor = [self.get_species()], product = self.products, component = self, part_id = self.name+"_degredation")

        return species

    def update_reactions(self) -> List:
        reactions = []
        mech_pathway = self.get_mechanism('metabolic_pathway', optional_mechanism = True)

        if mech_pathway is not None:
            if len(self.precursors) > 0:
                reactions += mech_pathway.update_reactions(precursor = self.precursors, product = [self.get_species()], component = self, part_id = self.name+"_production")
            if len(self.products) > 0:
                reactions += mech_pathway.update_reactions(precursor = [self.get_species()], product = self.products, component = self, part_id = self.name+"_degredation")
        return reactions


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
    """A class to represent Enzymes with multiple substrates and products.

    Assumes the enzyme converts all substrates to a all products at once.
    For example: S1 + S2 + ... + S_N + E --> P1 + P2 + ... + P_M + E.
    For enzymes with multiple enzymatic reactions, create multiple Enzyme Components with the same internal species.
    Uses a mechanism called "catalysis"
    """
    def __init__(self, enzyme: Union[Species, str, Component],
                 substrates: List[Union[Species, str, Component]],
                 products: List[Union[Species, str, Component]], attributes=None, **keywords):
        """Initialize an MultiEnzyme object to store MultiEnzyme related information.

        :param enzyme: name of the enzyme or reference to an Species or Component
        :param substrates: list of (name of the substrate or reference to an Species or Component)
        :param products: list of (name of the product or reference to an Species or Component)
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        self.enzyme = self.set_species(enzyme, material_type='protein', attributes=attributes)
        self.substrates = substrates
        self.products = products

        Component.__init__(self=self, name=self.enzyme.name, **keywords)

    @property
    def substrates(self) -> List:
        return self._substrates

    @substrates.setter
    def substrates(self, new_substrates: List[Union[Species, str, Component]]):
        if not isinstance(new_substrates, list):
            new_substrates = [new_substrates]
        # convert the new substrates to Species
        self._substrates = [self.set_species(s) for s in new_substrates]

    @property
    def products(self) -> List:
        return self._products

    @products.setter
    def products(self, new_products: List[Union[Species, str, Component]]):
        if not isinstance(new_products, list):
            new_products = [new_products]
        # convert the new products to Products
        self._products = [self.set_species(p) for p in new_products]
    
    def get_species(self) -> Species:
        return self.enzyme

    def update_species(self) -> List[Species]:
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(Enzyme=self.enzyme, Sub=self.substrates, Prod=self.products)

    def update_reactions(self) -> List[Reaction]:
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_reactions(Enzyme=self.enzyme, Sub=self.substrates, Prod=self.products,
                                         component=self,  part_id=self.name)
