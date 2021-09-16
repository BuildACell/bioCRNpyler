
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from .component import Component
from .reaction import Reaction
from .species import Complex, Species

class IntegralMembraneProtein(Component):
    """A class to represent transmembrane proteins or integral membrane proteins.
    This membrane class is to classify a membrane channel that will intergrate into the membrane.
    Uses a mechanism called "catalysis".
    Size is used to indicate number of repeating components to create oligomer. Dimer =2, Trimers = 3, etc."
    """
    def __init__(self, membrane_protein: Union[Species, str, Component],
                 product: Union[Species,str, Component],
                 direction= None, size=None, attributes=None, **keywords):       
        """Initialize a MembraneChannel object to store membrane channel related information.
        :param product: name of the membrane channel, reference to an Species or Component
        :param direction: transport direction, taken as "Passive"-undirectional unless specified 
        :param size: number of monomers needed for channel
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        
        # PROTEIN
        self.membrane_protein = self.set_species(membrane_protein, material_type='protein', attributes=attributes)
    
        # PRODUCT is an integrated membrane protein (transmembrane_protein)
        if product is None:
            if direction is None:
                self.product = self.set_species(product, material_type= 'Passive')
            else:
                self.product = self.set_species(product, material_type= direction)

        else:
            if direction is None:
                self.product = self.set_species(product, material_type= 'Passive')
            else:
                self.product = self.set_species(product, material_type= direction)
        
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

  class MembraneChannel(Component):
    """A class to represent Membrane Channels.
    Assumes the membrane channel transport substrates in a specific direction across the membrane
    Uses a mechanism called "catalysis"
    """
    def __init__(self, integral_membrane_protein: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 direction=None, attributes=None, **keywords):
        """Initialize a Transporter object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param direction: give direction of transport ref to vesicle
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
        
    # PROTEIN
        if type(integral_membrane_protein) == str:
            if direction == None:
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='Passive', attributes=attributes)
            else:
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type=direction, attributes=attributes)
                if direction == 'Importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product=substrate
                        self.substrate = self.set_species(substrate, compartment='External',attributes=attributes)
                        self.product= self.set_species(product, compartment='Internal',attributes=attributes)
        else:
            if integral_membrane_protein.material_type == 'Passive':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='Passive', attributes=attributes)
            elif integral_membrane_protein.material_type == 'Exporter':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='Exporter', attributes=attributes)
            elif integral_membrane_protein.material_type == 'Importer':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='Importer', attributes=attributes)

                if substrate is None:
                    self.substrate = None
                else:
                    product=substrate
                    self.substrate = self.set_species(substrate, compartment='External',attributes=attributes)
                    self.product= self.set_species(product, compartment='Internal',attributes=attributes)

            else:
                print('Membrane channel direction not found.')


        Component.__init__(self=self, name=self.integral_membrane_protein.name, **keywords)
        
        #######################
        print(self.attributes)
        
    def get_species(self):
        return self.integral_membrane_protein

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_species(self.integral_membrane_protein, self.substrate, self.product) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_reactions(self.integral_membrane_protein, self.substrate, self.product, component=self,  part_id=self.name)    


class MembranePump(Component):
    """A class to represent Membrane Channels.
    Assumes the membrane channel transport substrates in a specific direction across the membrane
    Uses a mechanism called "catalysis"
    """
    def __init__(self, membrane_pump: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 direction=None, ATP=None,
                 attributes=None, **keywords):
        """Initialize a Transporter object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param direction: give direction of transport ref to vesicle
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

# PROTEIN
        if type(membrane_pump) == str:
            if ATP is None:
                ATP= 1
            else: 
                ATP = ATP
                
            if direction is  None:
                print('None')
                self.membrane_pump = self.set_species(membrane_pump, material_type='Passive', attributes=attributes)
                self.membrane_pump.ATP= ATP
            else:
                self.membrane_pump = self.set_species(membrane_pump, material_type=direction, attributes=attributes)
                self.membrane_pump.ATP= ATP
                if direction == 'Importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product=substrate
                        self.substrate = self.set_species(substrate, compartment='External',attributes=attributes)
                        self.product= self.set_species(product, compartment='Internal',attributes=attributes)

                
        else:
            if integral_membrane_protein.material_type == 'Passive':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='Passive', attributes=attributes)
            elif membrane_pump.material_type == 'Exporter':
                self.membrane_pump = self.set_species(membrane_pump, material_type='Exporter', attributes=attributes)
            elif membrane_pump.material_type == 'Importer':
                self.membrane_pump = self.set_species(membrane_pump, material_type='Importer', attributes=attributes)

                if substrate is None:
                    self.substrate = None
                else:
                    product=substrate
                    self.substrate = self.set_species(substrate, compartment='External',attributes=attributes)
                    self.product= self.set_species(product, compartment='Internal',attributes=attributes)

            else:
                print('Membrane channel direction not found.')
                
            if ATP is None:
                self.membrane_pump.ATP= 1
            else: 
                self.membrane_pump.ATP = ATP

            
        self.energy = self.set_species('ATP',  material_type='small_molecule', compartment='Internal',attributes=attributes)
        self.waste = self.set_species('ADP',  material_type='small_molecule', compartment='Internal',attributes=attributes)

        Component.__init__(self=self, name=self.membrane_pump.name, **keywords)
        
        #######################
        print(self.attributes)
        
    def get_species(self):
        return self.membrane_pump

    def update_species(self):
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_species(self.membrane_pump, self.substrate, self.product, self.energy, self.waste) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('catalysis')
        return mech_cat.update_reactions(self.membrane_pump, self.substrate, self.product, self.energy, self.waste, component=self,  part_id=self.name)    