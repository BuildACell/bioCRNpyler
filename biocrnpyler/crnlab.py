
#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import inspect
import sys
import warnings

import libsbml

from .dna_assembly import DNAassembly
from .mixtures_extract import Mixture


class CRNLab(object):
    '''
    Implements high level modeling akin to experimental TX-TL experiments
    '''
    def __init__(self, name = '', **kwargs):
        self.Mixture = Mixture
        self.name = name
        self.volume = 0
        self.crn = None
        self.warning_print = True

        warnings.warn("CRNLab is deprecated and will cease to function with future releases.")


    def mixture(self, name, **kwargs):
        """
        Create a Mixture of a given name into the CRNLab object
        Specify the extract and buffer optionally
        Specify extra parameters to be loaded as dictionaries optionally
        """

        warnings.warn("CRNLab is deprecated and will cease to function with future releases: to produce a Mixture please use that Mixture's class constructor.")

        extract = kwargs.get('extract') 
        if 'warning_print' in kwargs:
            self.warning_print = kwargs['warning_print']
        if not extract:
            if self.warning_print or kwargs['parameter_warnings']:
                warnings.warn('The extract argument not given, using the BasicExtract by default.')
            extract = 'BasicExtract'
        mixture_parameters = kwargs.get('mixture_parameters')
        if not mixture_parameters:
            if self.warning_print or kwargs['parameter_warnings']:
                warnings.warn('Using default parameters for the mixture {0}'.format(self.name))
        elif 'mixture_parameters' in kwargs:
            kwargs['parameter_file'] = mixture_parameters
        if kwargs.get('mixture_volume'):
            self.volume += kwargs.get('mixture_volume')
        if kwargs.get('final_volume'):
            self.volume = kwargs.get('final_volume')

        if kwargs.get('final_volume') and kwargs.get('mixture_volume'):
            raise ValueError('Either set initial volume or the final volume')
        elif not kwargs.get('final_volume') and not kwargs.get('mixture_volume'):
            if self.warning_print or kwargs['parameter_warnings']:
                warnings.warn('Default volume of 1 uL will be set for {0}'.format(self.name))
            self.volume += 1e-6
        
        extracts_classes = inspect.getmembers(sys.modules[__name__], inspect.isclass)
        for class_t in extracts_classes:
            if class_t[0] == extract:
                # Call the appropriate class with the arguments
                # and get the Mixture object back to update `self.Mixture`
                extract_mix = class_t[1](self.name, **kwargs)
        self.Mixture = extract_mix
        return self.Mixture

   
    def add_dna(self, dna = None, name = "", promoter = "", 
                rbs = "", protein = "", initial_conc ="", 
                final_conc = "", volume = 0, **kwargs):

        warnings.warn("CRNLab is deprecated and will cease to function with future releases: to add DNA to a mixture, please use Mixture.add_component.")

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
        warnings.warn("CRNLab is deprecated and will cease to function with future releases: to add a Component to a mixture, please use Mixture.add_component.")

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
        warnings.warn("CRNLab is deprecated and will cease to function with future releases: to generate a CRN model, please use Mixture.combile_crn()")
        self.crn = self.Mixture.compile_crn()
        return self.crn

    def write_sbml_file(self, filename, **kwargs):
        warnings.warn("CRNLab is deprecated and will cease to function with future releases: to save a CRN to SBML, please use ChemicalReactionNetwork.write_sbml_file()")
        self.set_volumes()
        if self.volume:
            kwargs['volume'] = self.volume
            document, _ = self.crn.generate_sbml_model(**kwargs)
        else:
            document, _ = self.crn.generate_sbml_model(**kwargs)
        if document.getNumErrors() and self.warning_print:
            warnings.warn('SBML document has errors. It is recommended that you fix them before generating this model.')
        status = libsbml.writeSBML(document, filename)
        if status == libsbml.LIBSBML_OPERATION_SUCCESS:
            print('SBML file written successfully to {0}'.format(filename))
        return document 

    def validate_sbml_generated(self, **kwargs):
        document, _ = self.crn.generate_sbml_model(**kwargs)
        if document.getNumErrors() and self.warning_print:
            warnings.warn('SBML document has errors. It is recommended that you fix them before generating this model.')
            print('Here are the errors')
            print(document.getErrorLog())
        else:
            print('The SBML model generated is a valid document according to the SBML Level 3 Version 1 specifications. Use sbml.org/validator for further troubleshooting.')

# TODO : 
# Need to create more CRNLab examples to create models for other stuff. 
# Try to copy stuff over from MATLAB txtl to implement examples there and recreate the results.  
# Add unit tests for CRNLab
