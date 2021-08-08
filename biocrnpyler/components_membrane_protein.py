
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from .component import Component
from .reaction import Reaction
from .species import Complex, Species


class MembraneChannel(Component):
    """A class to represent Membrane Channels.
    This membrane class is to classify a membrane channel that will intergrate into the membrane.
    Uses a mechanism called "catalysis".
    Size is used to indicate number of repeating components to create oligomer. Dimer =2, Trimers = 3, etc."
    """
    def __init__(self, membrane_protein: Union[Species, str, Component],
                 product: Union[Species,str, Component],
                 size=None, attributes=None, **keywords):   
        """Initialize a MembraneChannel object to store membrane channel related information.
        :param product: name of the membrane channel, reference to an Species or Component
        :param size: number of monomers needed for channel
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """    
        
        # PROTEIN
        self.membrane_protein = self.set_species(membrane_protein, material_type='protein', attributes=attributes)
    
        # PRODUCT
        if product is None:
            self.product = self.set_species('membrane_channel')
        else:
            self.product = self.set_species(product)
        
        if size is None:
            self.membrane_protein.size = 1
        else:
            self.membrane_protein.size = size
        
        Component.__init__(self=self, name=self.membrane_protein.name, **keywords)
        
    def get_species(self):
        return self.membrane_protein
    
    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')
        
        return mech_cat.update_species(self.membrane_protein, self.product)

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')
        print(self)
        
        return mech_cat.update_reactions(self.membrane_protein, self.product, component=self,  part_id=self.name)

class Transporter(Component):
    """A class to represent Membrane Channels.
    Assumes the membrane channel transport substrates in a specific direction across the membrane
    Uses a mechanism called "catalysis"
    """
    def __init__(self, membrane_channel: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 attributes=None, direction= None, **keywords):
        """Initialize a Transporter object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """

   # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            product=substrate
            self.substrate = self.set_species(substrate, compartment='Internal',attributes=attributes)
            self.product= self.set_species(product, compartment='External',attributes=attributes)
      
        if direction is None:
            # Protein NAME
            self.membrane_channel = self.set_species(membrane_channel, material_type='Passive', attributes=attributes)

        elif direction== 'Exporter':
            # Protein NAME
            self.membrane_channel = self.set_species(membrane_channel, material_type='Exporter', attributes=attributes)

        elif direction== 'Importer':
            # Protein NAME
            self.membrane_channel = self.set_species(membrane_channel, material_type='Importer', attributes=attributes)
            
            if substrate is None:
                self.substrate = None
            else:
                product=substrate
                self.substrate = self.set_species(substrate, compartment='External',attributes=attributes)
                self.product= self.set_species(product, compartment='Internal',attributes=attributes)
        
        else:
            print('help')


        Component.__init__(self=self, name=self.membrane_channel.name, **keywords)
        
        #######################
        print(self.attributes)
        
    def get_species(self):
        return self.membrane_channel

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_species(self.membrane_channel, self.substrate, self.product) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')

        return mech_cat.update_reactions(self.membrane_channel, self.substrate, self.product, component=self,  part_id=self.name)      