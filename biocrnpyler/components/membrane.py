
#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List, Union
import warnings

from ..core.component import Component
from ..core.reaction import Reaction
from ..core.species import Complex, Species
from ..core.compartment import Compartment

class DiffusibleMolecule(Component):
    """A class to represent transmembrane proteins or integral membrane proteins.
    This membrane class is to classify a membrane channel that will intergrate into the membrane.
    By default, a DiffusibleMolecule uses a mechanism called "diffusion".
    """
    def __init__(self, substrate:Union[Species, str, Component], 
                 internal_compartment:Union[str, Compartment]='Internal', 
                 external_compartment:Union[str, Compartment]='External',
                 attributes=None, **keywords):       
        """Initialize a DiffusibleMolecule object to store molecule/substance related information.
        :param substrate: name of the diffusible substrate, reference to an Species or Component
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :param attributes: Species attribute, passed to Component
        :param keywords: pass into the parent's (Component) initializer
        """
        # Creates compartment object if compartment is a str
        if type(internal_compartment) is str:
            internal_compartment = Compartment(name=internal_compartment)
        if type(external_compartment) is str:
            external_compartment = Compartment(name=external_compartment)
        
        # Substrate
        self.substrate = self.set_species(substrate, compartment=internal_compartment)
        self.product = self.set_species(substrate, compartment=external_compartment)
        
        # Name the component
        name = self.substrate.name+'_'+self.substrate.compartment.name 

        Component.__init__(self=self, name=name, attributes=attributes, **keywords)
        
    def get_species(self):
        return self.substrate
    
    def update_species(self):
        mech_diff = self.get_mechanism('diffusion')
        return mech_diff.update_species(self.substrate, self.product)

    def update_reactions(self):
        mech_diff = self.get_mechanism('diffusion')
        return mech_diff.update_reactions(self.substrate, self.product, component=self,  part_id=self.name)

class IntegralMembraneProtein(Component):
    """A class to represent transmembrane proteins or integral membrane proteins.
    This membrane class is to classify a membrane channel that will intergrate into the membrane.
    Uses a mechanism called "membrane_insertion".
    Size is used to indicate number of repeating components to create oligomer. Dimer = 2, Trimers = 3, etc.
    """
    def __init__(self, membrane_protein:Union[Species, str, Component],
                 product:Union[Species,str, Component], direction:str=None, size:int=None, 
                 compartment:Union[str, Compartment]= 'Internal', 
                 membrane_compartment:Union[str, Compartment]='Membrane', 
                 attributes=None, **keywords):
        """Initialize a IntegralMembraneProtein object to store membrane channel related information.
        :param product: name of the membrane channel, reference to an Species or Component
        :param direction: transport direction (str), set to "Passive" by default, undirectional unless specified 
        :param size: number of monomers needed for channel used in Membrane_Protein_Integration(Mechanism)
        :param internal_compartment: name of internal compartment 
        :param membrane_compartment: name of membrane compartment 
        :param attributes: Species attribute.
        :param keywords: pass into the parent's (Component) initializer
        """      
        # Creates compartment object if compartment is a str
        if type(compartment) is str:
            compartment = Compartment(name=compartment)
        if type(membrane_compartment) is str:
            membrane_compartment = Compartment(name=membrane_compartment)
        
        # Logic for prioritizing compartments
        self.membrane_protein = self.set_species(membrane_protein)
        if self.membrane_protein.compartment.name == 'default':
          self.membrane_protein.compartment = compartment
        elif self.membrane_protein.compartment.name != compartment.name and compartment.name == 'Internal':
          warnings.warn("Inconsistent compartments, prioritizing membrane protein compartment.", UserWarning)
          compartment = self.membrane_protein.compartment
        else:
          warnings.warn("Inconsistent compartments, prioritizing integral membrane protein compartment.", UserWarning)
          self.membrane_protein.compartment = compartment

        # PROTEIN
        self.membrane_protein = self.set_species(membrane_protein, material_type='protein', 
                                                 compartment=compartment,
                                                 attributes=attributes)

        # PRODUCT is an integrated membrane protein (transmembrane_protein)
        if product is None:
            if direction is None:
                self.product = self.set_species(product, material_type= 'protein', 
                                                compartment=membrane_compartment, 
                                                attributes=['Passive'])
            else:
                self.product = self.set_species(product, material_type= 'protein', 
                                                compartment=membrane_compartment, 
                                                attributes=[direction])
        else:
            if direction is None:
                self.product = self.set_species(product, material_type='protein', 
                                                compartment=membrane_compartment, 
                                                attributes=['Passive'])
            else:
                self.product = self.set_species(product, material_type= 'protein', 
                                                compartment=membrane_compartment, 
                                                attributes=[direction])

        # Indicates the number of monomers that compose the channel, will be used in Membrane_Protein_Integration(Mechanism)
        if size is None:
            self.membrane_protein.size = 1
        else:
            self.membrane_protein.size = size
        
        # Name the component
        name = self.membrane_protein.name+'_'+self.membrane_protein.compartment.name

        Component.__init__(self=self, name=name, **keywords)
        
    def get_species(self):
        return self.membrane_protein
    
    def update_species(self):
        mech_ins = self.get_mechanism('membrane_insertion')
        return mech_ins.update_species(self.membrane_protein, self.product)

    def update_reactions(self):
        mech_ins = self.get_mechanism('membrane_insertion')
        return mech_ins.update_reactions(self.membrane_protein, self.product, component=self,  part_id=self.name)

class MembraneChannel(Component):
    """A class to represent membrane channels.
    
    The membrane channel transports substrates across the membrane
    following the concentration gradient.  Direction and mechanism will be
    based on the specific transporter.
    
    Uses a mechanism called "transport".

    """
    def __init__(
            self, integral_membrane_protein: Union[Species, str, Component],
            substrate: Union[Species, str, Component],direction:str=None, 
            internal_compartment:Union[str, Compartment]='Internal', 
            external_compartment:Union[str, Compartment]='External',
            attributes=None, **keywords):
        """Initialize a MembraneChannel object to store transport membrane
        related information.
        
        :param substrate: substrate to be transported (str, Species, Component)
        :param direction: direction of transport based on transporter action
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer

        """
        # Creates compartment object if compartment is a str
        if type(internal_compartment) is str:
            internal_compartment = Compartment(name=internal_compartment)
        if type(external_compartment) is str:
            external_compartment = Compartment(name=external_compartment)

        # Set up the integral membrane protein
        if isinstance(integral_membrane_protein, str):
            integral_membrane_protein = self.set_species(
                integral_membrane_protein, material_type='protein', 
                attributes='Passive' if direction is None else direction)
        self.integral_membrane_protein = integral_membrane_protein

        # Get the direction from the integral_membrane_protein, if not given
        # TODO: need more complete check for conflicting information
        if direction is None:
            if 'Importer' in integral_membrane_protein.attributes:
                direction = 'Importer'
            elif 'Exporter' in integral_membrane_protein.attributes:
                direction = 'Exporter'

        # Substrate and product assignments.
        #
        # In the case of membrane components, the substrate is the
        # substance on which the transporter/channel acts without
        # distinction of compartment.  The substrate and product are the
        # same substance and the substance does not change as a result
        # except for the compartment.  Therefore, the product here is based
        # on the action of the transporter.
        #
        # TODO: if the substrate is not passed as a string, we need more
        # sophisticated logic to set up the product, since the `set_species`
        # method will just return the existing species, without changing
        # the compartment.

        if substrate is None:
            self.substrate = None

        else:
            if direction == 'Importer':
                self.substrate = self.set_species(
                    substrate, compartment=external_compartment,
                    attributes=attributes)
                self.product= self.set_species(
                    self.substrate.name, compartment=internal_compartment,
                    attributes=attributes)
            else:
                self.substrate = self.set_species(
                    substrate, compartment=internal_compartment,
                    attributes=attributes)
                self.product= self.set_species(
                    self.substrate.name, compartment=external_compartment,
                    attributes=attributes)
        
        # Name the component
        name = self.integral_membrane_protein.name + '_' + \
            self.integral_membrane_protein.compartment.name
        
        Component.__init__(self=self, name=name, **keywords)
        
    def get_species(self):
        return self.integral_membrane_protein

    def update_species(self):
        mech_tra = self.get_mechanism('transport')
        return mech_tra.update_species(self.integral_membrane_protein, self.substrate, self.product) 

    def update_reactions(self):
        mech_tra = self.get_mechanism('transport')
        return mech_tra.update_reactions(self.integral_membrane_protein, self.substrate, self.product, 
                                         component=self,  part_id=self.name)    

class MembranePump(Component):
    """A class to represent membrane pumps or transporters that require ATP.
    The membrane pump transports substrates unidirectionally across the membrane, independent of the concentration gradient.
    Uses a mechanism called "transport".
    """
    def __init__(self, membrane_pump: Union[Species, str, Component],
                 substrate: Union[Species, str, Component], direction:str=None, 
                 internal_compartment:Union[str, Compartment]='Internal', 
                 external_compartment:Union[str, Compartment]='External',
                 ATP:int=None, attributes=None, **keywords):
        """Initialize a MembranePump object to store Transport membrane related information.
        :param substrate: name of the substrate, reference to an Species or Component
        :param direction: give direction of transport ref to vesicle
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :param ATP: indicates the number of ATP required for transport
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        # Creates compartment object if compartment is a str
        if type(internal_compartment) is str:
            internal_compartment = Compartment(name=internal_compartment)
        if type(external_compartment) is str:
            external_compartment = Compartment(name=external_compartment)

        # SUBSTRATE
        if substrate is None:
            self.substrate = None
        else:
            product=substrate
            self.substrate = self.set_species(substrate, 
                                              compartment=internal_compartment,
                                              attributes=attributes)
            self.product= self.set_species(product,
                                           compartment=external_compartment,
                                           attributes=attributes)

        # ENERGY and WASTE
        self.energy = self.set_species('ATP',  material_type='small_molecule', 
                                       compartment=internal_compartment,
                                       attributes=attributes)
        self.waste = self.set_species('ADP',  material_type='small_molecule', 
                                       compartment=internal_compartment,
                                       attributes=attributes)

        # PROTEIN
        if type(membrane_pump) == str:
            if ATP is None:
                ATP= 1
            else: 
                ATP = ATP
                
            if direction is  None:
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', 
                                                      attributes='Passive')
                self.membrane_pump.ATP= ATP
            else:
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein', 
                                                      attributes=direction)
                self.membrane_pump.ATP= ATP
                if direction == 'Importer':
                    if substrate is None:
                        self.substrate = None
                    else:
                        product=substrate
                        self.substrate = self.set_species(substrate, 
                                                          compartment=external_compartment,
                                                          attributes=attributes)
                        self.product= self.set_species(product,
                                                       compartment=internal_compartment,
                                                       attributes=attributes)
        else:
            if membrane_pump.attributes[0] == 'Passive':
                self.integral_membrane_protein = self.set_species(membrane_pump, material_type='protein', 
                                                                  attributes='Passive')
            elif membrane_pump.attributes[0] == 'Exporter':
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein',
                                                      attributes='Exporter')
            elif membrane_pump.attributes[0] == 'Importer':
                self.membrane_pump = self.set_species(membrane_pump, material_type='protein',
                                                      attributes='Importer')
                self.energy = self.set_species('ATP',  material_type='small_molecule', 
                                                compartment=external_compartment,
                                                attributes=attributes)
                self.waste = self.set_species('ADP',  material_type='small_molecule', 
                                                compartment=external_compartment,
                                                attributes=attributes)
                if substrate is None:
                    self.substrate = None
                else:
                    product=substrate
                    self.substrate = self.set_species(substrate, 
                                                      compartment=external_compartment,
                                                      attributes=attributes)
                    self.product= self.set_species(product, 
                                                   compartment=internal_compartment,
                                                   attributes=attributes)
            else:
                print('Membrane channel direction not found.')
                
            if ATP is None:
                self.membrane_pump.ATP= 1
            else: 
                self.membrane_pump.ATP = ATP
        
        # Name the component
        name = self.membrane_pump.name+'_'+self.membrane_pump.compartment.name

        Component.__init__(self=self, name=name, **keywords)
        
    def get_species(self):
        return self.membrane_pump

    def update_species(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_species(self.membrane_pump, self.substrate, self.product, self.energy, self.waste) 

    def update_reactions(self):
        mech_cat = self.get_mechanism('transport')
        return mech_cat.update_reactions(self.membrane_pump, self.substrate, self.product, self.energy, self.waste, 
                                         component=self,  part_id=self.name)    

class MembraneSensor(Component):
    """A class to represent a two-component system (TCS) membrane sensor.
    The membrane sensor protein senses the signal substrate and added the assigned substrate to the response protein.
    Uses a mechanism called "membrane_sensor".
    """
    def __init__(self, membrane_sensor_protein: Union[Species, str, Component],
                 response_protein: Union[Species, str, Component],
                 assigned_substrate: Union[Species, str, Component],
                 signal_substrate: Union[Species, str, Component],
                 product: Union[Species, str, Component]=None,
                 internal_compartment:Union[str, Compartment]='Internal', 
                 external_compartment:Union[str, Compartment]='External',
                 ATP:int=2, attributes=None, **keywords):
        """Initialize a MembraneSensor object to store Transport membrane related information.
        :param membrane_sensor_protein: name of the membrane protein in the TCS, reference to an Species or Component
        :param response_protein: name of the response protein in the TCS, reference to an Species or Component
        :param assigned_substrate: name of the assigned substrate in the TCS, reference to an Species or Component
        :param signal_substrate: name of the signal substrate in the TCS, reference to an Species or Component
        :param product: name of the product in the TCS, reference to an Species or Component
        :param internal_compartment: name of internal compartment 
        :param external_compartment: name of external compartment 
        :param ATP: indicates the number of ATP required for transport
        :param attributes: Species attribute
        :param keywords: pass into the parent's (Component) initializer
        """
        # Creates compartment object if compartment is a str
        if type(internal_compartment) is str:
            internal_compartment = Compartment(name=internal_compartment)
        if type(external_compartment) is str:
            external_compartment = Compartment(name=external_compartment)
                    
        # RESPONSE PROTEIN
        if response_protein is None:
            self.response_protein = None
        else:
            self.response_protein = self.set_species(response_protein, compartment=internal_compartment,
                                                     attributes=attributes)

        # PRODUCT PROTEIN
        if product is None:
            self.product = self.set_species(str(response_protein)+'active', compartment=internal_compartment,
                                            attributes=attributes)
        else:
            self.product = self.set_species(product, compartment=internal_compartment,
                                            attributes=attributes)

        # ASSIGNED SUBSTRATE
        if assigned_substrate is None:
            self.assigned_substrate = None
        else:
            self.assigned_substrate = self.set_species(assigned_substrate, compartment=internal_compartment,
                                                       attributes=attributes)
        # SIGNAL SUBSTRATE
        if signal_substrate is None:
            self.signal_substrate = None
        else:
            self.signal_substrate = self.set_species(signal_substrate, compartment=internal_compartment,
                                                     attributes=attributes)
    
        # PROTEIN
        if membrane_sensor_protein is None:
            self.membrane_sensor_protein = None
        else:
            self.membrane_sensor_protein = self.set_species(membrane_sensor_protein, material_type='protein', 
                                                            attributes=attributes)
        # ENERGY: ATP
        if ATP is None:
            self.membrane_sensor_protein.ATP= 1
        else: 
            self.membrane_sensor_protein.ATP = ATP

        self.energy = self.set_species('ATP',  material_type='small_molecule', compartment=internal_compartment,
                                       attributes=attributes)
        self.waste = self.set_species('ADP',  material_type='small_molecule', compartment=internal_compartment,
                                      attributes=attributes)
        # Name the component
        name = self.membrane_sensor_protein.name+'_'+self.membrane_sensor_protein.compartment.name

        Component.__init__(self=self, name=name, **keywords)

    def get_species(self):
        return self.membrane_sensor_protein

    def update_species(self):
        mech_sen = self.get_mechanism('membrane_sensor')
        return mech_sen.update_species(self.membrane_sensor_protein, self.response_protein, 
                                       self.assigned_substrate, self.signal_substrate, self.product,
                                       self.energy, self.waste) 

    def update_reactions(self):
        mech_sen = self.get_mechanism('membrane_sensor')
        return mech_sen.update_reactions(self.membrane_sensor_protein, self.response_protein, 
                                         self.assigned_substrate, self.signal_substrate, self.product,
                                         self.energy, self.waste, component=self,  part_id=self.name)   
