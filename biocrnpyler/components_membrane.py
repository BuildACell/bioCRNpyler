
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union

from .component import Component
from .reaction import Reaction
from .species import Complex, Species

class IntegralMembraneProtein(Component):
    """A class to represent transmembrane proteins or integral membrane proteins.
    This membrane class is to classify a membrane channel that will intergrate into the membrane.
    Uses a mechanism called "membrane_insertion".
    Size is used to indicate number of repeating components to create oligomer. Dimer =2, Trimers = 3, etc."
    """
    def __init__(self, membrane_protein: Union[Species, str, Component],
                 product: Union[Species,str, Component],
                 direction= None, size=None, attributes=None, **keywords):       
        """Initialize a MembraneChannel object to store membrane channel related information.
        :param product: name of the membrane channel, reference to an Species or Component
        :param direction: transport direction, taken as "Passive"-undirectional unless specified 
        :param size: number of monomers needed for channel used in Membrane_Protein_Integration(Mechanism)
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        
    # PROTEIN
        self.membrane_protein = self.set_species(membrane_protein, material_type='protein', attributes=attributes)
    
    # PRODUCT is an integrated membrane protein (transmembrane_protein)
        if product is None:
            if direction is None:
                self.product = self.set_species(product, material_type= 'protein', attributes=['Passive'])
            else:
                self.product = self.set_species(product, material_type= 'protein', attributes=[direction])

        else:
            if direction is None:
                self.product = self.set_species(product, material_type='protein', attributes=['Passive'])
            else:
                self.product = self.set_species(product, material_type= 'protein', attributes=[direction])

    #Indicates the number of monomers that compose the channel, will be used in Membrane_Protein_Integration(Mechanism)
        if size is None:
            self.membrane_protein.size = 1
        else:
            self.membrane_protein.size = size
        
        Component.__init__(self=self, name=self.membrane_protein.name, **keywords)
        
    def get_species(self):
        return self.membrane_protein
    
    def update_species(self):
        mech_cat = self.get_mechanism('membrane_insertion')
        
        return mech_cat.update_species(self.membrane_protein, self.product)

    def update_reactions(self):
        mech_cat = self.get_mechanism('membrane_insertion')
        print(self)
        
        return mech_cat.update_reactions(self.membrane_protein, self.product, component=self,  part_id=self.name)

class MembraneChannel(Component):
    """A class to represent Membrane Channels.
    Assumes the membrane channel transport substrates in a specific direction across the membrane
    Uses a mechanism called "transport"
    """
    def __init__(self, integral_membrane_protein: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 direction=None, internal_compartment='Internal', external_compartment='External',
                 cell=None, attributes=None, **keywords):
        """Initialize a Transporter object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param direction: direction of transport based on Transporter action
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :parm cell: indicated the cell, identifier can be name or number
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
    #Additional information on the identity of the specific cell/vesicle (if needed).
        if cell is not None:
            if type(cell) is str:
                internal_compartment=internal_compartment+'_'+cell
                external_compartment=external_compartment+'_'+cell
            else:
                internal_compartment=internal_compartment+'_'+str(cell)
                external_compartment=external_compartment+'_'+str(cell) 

    #Substrate and product assignments.
        """In the case of membrane components, the substrate is the substance on which the transporter/channel acts without distinction of compartment. 
        The substrate and product are the same substance and the substance does not change as a result expect for the compartment.
        Therefore, the product here is based on the action of the transporter."""
        if substrate is None:
            self.substrate = None
        else:
            product=substrate
            self.substrate = self.set_species(substrate, compartment=internal_compartment)
            self.product= self.set_species(product, compartment=external_compartment)
    
    # PROTEIN
        if type(integral_membrane_protein) == str:
            if direction == None:
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='protein', attributes='Passive')
            else:
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='protein', attributes=direction)
                if direction == 'Importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product=substrate
                        self.substrate = self.set_species(substrate, compartment=external_compartment,attributes=attributes)
                        self.product= self.set_species(product, compartment=internal_compartment,attributes=attributes)
        else:
            if integral_membrane_protein.attributes[0] == 'Passive':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='protein', attributes='Passive')
            elif integral_membrane_protein.attributes[0] == 'Exporter':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='protein', attributes='Exporter')
            elif integral_membrane_protein.attributes[0] == 'Importer':
                self.integral_membrane_protein = self.set_species(integral_membrane_protein, material_type='protein', attributes='Importer')

                if substrate is None:
                    self.substrate = None
                else:
                    product=substrate
                    self.substrate = self.set_species(substrate, compartment=external_compartment,attributes=attributes)
                    self.product= self.set_species(product, compartment=internal_compartment,attributes=attributes)

            else:
                print('Membrane channel direction not found.')

        Component.__init__(self=self, name=self.integral_membrane_protein.name, **keywords)
        
    def get_species(self):
        return self.integral_membrane_protein

    def update_species(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_species(self.integral_membrane_protein, self.substrate, self.product) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_reactions(self.integral_membrane_protein, self.substrate, self.product, component=self,  part_id=self.name)    

class MembranePump(Component):
    """A class to represent Membrane Channels.
    Assumes the membrane channel transport substrates in a specific direction across the membrane
    Uses a mechanism called "transport"
    """
    def __init__(self, membrane_pump: Union[Species, str, Component],
                 substrate: Union[Species, str, Component],
                 direction=None, internal_compartment='Internal', external_compartment='External',
                 ATP=None, cell=None, attributes=None, **keywords):
        """Initialize a Transporter object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param direction: give direction of transport ref to vesicle
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :parm cell: indicated the cell, identifier can be name or number
        :param ATP: indicates the number of ATP required for transport
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
    #Additional information on the identity of the specific cell/vesicle (if needed).
        if cell is not None:
            if type(cell) is str:
                internal_compartment=internal_compartment+'_'+cell
                external_compartment=external_compartment+'_'+cell
            else:
                internal_compartment=internal_compartment+'_'+str(cell)
                external_compartment=external_compartment+'_'+str(cell) 
   
    # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            product=substrate
            self.substrate = self.set_species(substrate, compartment=internal_compartment,attributes=attributes)
            self.product= self.set_species(product, compartment=external_compartment,attributes=attributes)

    # PROTEIN
        if type(membrane_pump) == str:
            if ATP is None:
                ATP= 1
            else: 
                ATP = ATP
                
            if direction is  None:
                print('None')
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', attributes='Passive')
                self.membrane_pump.ATP= ATP
            else:
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', attributes=direction)
                self.membrane_pump.ATP= ATP
                if direction == 'Importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product=substrate
                        self.substrate = self.set_species(substrate, compartment=external_compartment,attributes=attributes)
                        self.product= self.set_species(product, compartment=internal_compartment,attributes=attributes)

        else:
            if membrane_pump.attributes[0] == 'Passive':
                self.integral_membrane_protein = self.set_species(membrane_pump, material_type='protein', attributes='Passive')
            elif membrane_pump.attributes[0] == 'Exporter':
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', attributes='Exporter')
            elif membrane_pump.attributes[0] == 'Importer':
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', attributes='Importer')

                if substrate is None:
                    self.substrate = None
                else:
                    product=substrate
                    self.substrate = self.set_species(substrate, compartment=external_compartment,attributes=attributes)
                    self.product= self.set_species(product, compartment=internal_compartment,attributes=attributes)

            else:
                print('Membrane channel direction not found.')
                
            if ATP is None:
                self.membrane_pump.ATP= 1
            else: 
                self.membrane_pump.ATP = ATP
            
        self.energy = self.set_species('ATP',  material_type='small_molecule', compartment=internal_compartment,attributes=attributes)
        self.waste = self.set_species('ADP',  material_type='small_molecule', compartment=internal_compartment,attributes=attributes)

        Component.__init__(self=self, name=self.membrane_pump.name, **keywords)
        
    def get_species(self):
        return self.membrane_pump

    def update_species(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_species(self.membrane_pump, self.substrate, self.product, self.energy, self.waste) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_reactions(self.membrane_pump, self.substrate, self.product, self.energy, self.waste, component=self,  part_id=self.name)    