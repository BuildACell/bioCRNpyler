from .txtl import BasicExtract
from .mixture import Mixture
from .dna_assembly import DNAassembly
import csv
import libsbml

class CRNLab(object):
    '''
    Implements high level modeling akin to experimental TX-TL experiments
    '''
    def __init__(self, name = ''):
        self.Mixture = Mixture
        self.name = name
        self.volume = 0
        self.crn = None
    
    def mixture(self, name, **kwargs):        
        '''
        Create a Mixture of a given name into the CRNLab object
        Specify the extract and buffer optionally 
        Specify extra parameters to be loaded as dictionaries optionally
        '''
        extract = kwargs.get('extract') 
        buffer = kwargs.get('buffer')
        extract_parameters = kwargs.get('extract_parameters')
        extract_volume = kwargs.get('extract_volume')
        buffer_parameters = kwargs.get('buffer_parameters')
        buffer_volume = kwargs.get('buffer_volume')
        if kwargs.get('mixture_volume'):
            self.volume += kwargs.get('mixture_volume')
        if kwargs.get('final_volume'):
            self.volume = kwargs.get('final_volume')

        if kwargs.get('final_volume') and kwargs.get('mix_volume'):
            raise ValueError('Either set initial volume or the final volume')
        

        filename = name + str('.csv')
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            # Read off the parameters
            params = {}
            for row in reader:
                temp = row[1].replace('.','',1).replace('e','',1).replace('-','',1)
                if temp.isdigit():
                    params[str(row[0])] = float(row[1])
 
        if extract:
            self.extract(extract, parameters = extract_parameters, volume = extract_volume)
        if buffer:
            self.buffer(extract, parameters = buffer_parameters, volume = buffer_volume)
        return self.Mixture

    def extract(self, name, parameters = {}, volume = 0):
        '''
        Create BasicExtract with the given name 
        (Searches for the name.csv in the current folder to load parameters)
        Optionally load other parameters as dictionary.
        '''
        if volume:
            self.volume += volume
        # Look for extract config file of the given name
        filename = name + str('.csv')
        import csv
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            # Read off the parameters
            params = {}
            for row in reader:
                temp = row[1].replace('.','',1).replace('e','',1).replace('-','',1)
                if temp.isdigit():
                    params[str(row[0])] = float(row[1])
        if parameters:
            # If manually given
            params = parameters
            extract_mix = BasicExtract(self.name, parameters = params)
        else:
            extract_mix = BasicExtract(self.name, parameters = params)
        self.Mixture = extract_mix
        return self.Mixture

    def buffer(self, name = "", components = [], parameters = {}, volume = 0):
        '''
        TODO : To be implemented using energy models
        '''
        self.name = name
        if volume:
            self.volume += volume
        return 
   
    def add_dna(self, dna = None, name = "", promoter = "", rbs = "", protein = "", initial_conc ="", final_conc = "", volume = 0):
        if volume:
            self.volume += volume
        if dna:
            self.Mixture.add_components(dna)
        elif name and promoter and rbs and protein and initial_conc:
            if final_conc:
                dna = DNAassembly(name, promoter = promoter, rbs = rbs, protein = protein, initial_conc = final_conc)
            elif initial_conc:
                dna = DNAassembly(name, promoter = promoter, rbs = rbs, protein = protein, initial_conc = initial_conc*volume)
            self.add_dna(dna)
        else:
            raise ValueError('add_dna either needs a DNAassembly object OR a name, promoter, rbs, and protein etc.')
        return self.Mixture
    
    def add_component(self, component):
        self.Mixture.add_components(component)
        return self.Mixture
    
    def set_volumes(self):
        final_volume = self.volume
        # Not implemented yet
        if final_volume and False:
            print(self.Mixture.crn_species)
            for species in self.Mixture.crn_species:
                # Set the final concentration for all species
                # TODO
                species.initial_concentration = species.initial_concentration/final_volume
            self.crn = self.Mixture.compile_crn()
            return self.volume    
        else:
            return


    def combine_tubes(self):
        self.crn = self.Mixture.compile_crn()
        return self.crn

    def write_sbml_file(self, filename, **kwargs):
        self.set_volumes()
        if self.volume:
            kwargs['volume'] = self.volume
            document, _ = self.crn.generate_sbml_model(**kwargs)
        else:
            document, _ = self.crn.generate_sbml_model(**kwargs)
        sbml_string = libsbml.writeSBMLToString(document)
        f = open(filename, 'w')
        f.write(sbml_string)
        f.close()
        return f


