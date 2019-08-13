from .extracts import BasicExtract
from .mixture import Mixture
from .dna_assembly import DNAassembly
import csv
import libsbml
import warnings

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
        extract_parameters = kwargs.get('extract_parameters')
        extract_volume = kwargs.get('extract_volume')
        if kwargs.get('mixture_volume'):
            self.volume += kwargs.get('mixture_volume')
        if kwargs.get('final_volume'):
            self.volume = kwargs.get('final_volume')

        if kwargs.get('final_volume') and kwargs.get('mix_volume'):
            raise ValueError('Either set initial volume or the final volume')

        try:
            filename = name + str('.csv')
            with open(filename, 'r') as f:
                reader = csv.reader(f)
                # Read off the parameters
                params = {}
                for row in reader:
                    temp = row[1].replace('.','',1).replace('e','',1).replace('-','',1)
                    if temp.isdigit():
                        params[str(row[0])] = float(row[1])
        except:
            if kwargs.get('parameter_warnings'):
                warnings.warn('{0}.csv does not exist. Mixture will use default parameter files available.'.format(name))
            
        # if extract:
        self.extract(extract, parameters = extract_parameters, volume = extract_volume, **kwargs)
        # if buffer:
        #     self.buffer(extract, parameters = buffer_parameters, volume = buffer_volume, **kwargs)
        return self.Mixture

    def extract(self, name, volume = 0, **kwargs):
        '''
        Create BasicExtract with the given name 
        (Searches for the name.csv in the current folder to load parameters)
        Optionally load other parameters as dictionary.
        '''
        initial_concentration_dict = kwargs.get('initial_concentration_dict')
        if not initial_concentration_dict:
            # Set default values here
            initial_concentration_dict = {"protein_Ribo":10, "protein_RNAP":5, "protein_RNAase":2.5}
        if volume:
            self.volume += volume
        # Look for extract config file of the given name
        if kwargs.get('parameters'):
            parameters = kwargs.get('parameters')
            # If manually given
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
            params = parameters
            extract_mix = BasicExtract(self.name, parameters = params, 
                                        init = initial_concentration_dict, **kwargs)
        else:
            extract_mix = BasicExtract(self.name, init = initial_concentration_dict, **kwargs)
        self.Mixture = extract_mix
        return self.Mixture

   
    def add_dna(self, dna = None, name = "", promoter = "", rbs = "", protein = "", initial_conc ="", final_conc = "", volume = 0):
        if volume:
            self.volume += volume
        if dna:
            self.Mixture.add_components(dna)
        elif name and promoter and rbs and protein and initial_conc:
            if final_conc:
                dna = DNAassembly(name, promoter = promoter, rbs = rbs,
                                  protein = protein, initial_conc = final_conc)
            elif initial_conc:
                dna = DNAassembly(name, promoter = promoter, rbs = rbs,
                                  protein = protein,
                                  initial_conc = initial_conc*volume)
            self.add_dna(dna)
        else:
            raise ValueError('add_dna either needs a DNAassembly object '
                             'OR a name, promoter, rbs, and protein etc.')
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
                species.initial_concentration = \
                                    species.initial_concentration/final_volume
            self.crn = self.Mixture.compile_crn()
            return self.volume
        else:
            return


    def get_model(self):
        self.crn = self.Mixture.compile_crn()
        return self.crn

    def write_sbml_file(self, filename, **kwargs):
        self.set_volumes()
        if self.volume:
            kwargs['volume'] = self.volume
            document, _ = self.crn.generate_sbml_model(**kwargs)
        else:
            document, _ = self.crn.generate_sbml_model(**kwargs)
        if document.getNumErrors():
            warnings.warn('SBML document has errors. It is recommended that you fix them before generating this model.')
        status = libsbml.writeSBML(document, filename)
        if status == libsbml.LIBSBML_OPERATION_SUCCESS:
            print('SBML file written successfully to {0}'.format(filename))
        return document 

    def validate_sbml_generated(self, **kwargs):
        document, _ = self.crn.generate_sbml_model(**kwargs)
        if document.getNumErrors():
            warnings.warn('SBML document has errors. It is recommended that you fix them before generating this model.')
            print('Here are the errors')
            print(document.getErrorLog())
        else:
            print('The SBML model generated is a valid document according to the SBML Level 3 Version 1 specifications. Use sbml.org/validator for further troubleshooting.')


# TODO : 
# Need to test extensively with Parameter class and if its working.
# Need to create more CRNLab examples to create models for other stuff. 
# Try to copy stuff over from MATLAB txtl to implement examples there and recreate the results.  