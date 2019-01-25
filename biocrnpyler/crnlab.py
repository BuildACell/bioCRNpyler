from .txtl import BasicExtract
from .mixture import Mixture
from .dna_assembly import DNAassembly

class CRNLab(Mixture):
    '''
    Implements high level modeling akin to experimental TX-TL experiments
    '''
    def __init__(self, name):
        self.Mixture = Mixture
        self.name = name

    def extract(self, name, parameters = {}):
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

    def buffer(self, name = "", components = [], parameters = {}, **kwargs):
        '''
        TODO : To be implemented using energy models
        '''
        self.name = name
        return 
   
    def add_dna(self, dna = None, name = "", promoter = "", rbs = "", protein = "", initial_conc =""):
        if dna:
            self.Mixture.add_components(dna)
        elif name and promoter and rbs and protein and initial_conc:
            dna = DNAassembly(name, promoter = promoter, rbs = rbs, protein = protein, initial_conc = initial_conc)
            self.add_dna(dna)
        else:
            raise ValueError('add_dna either needs a DNAassembly object OR a name, promoter, rbs, and protein etc.')
        return self.Mixture

    def combine_tubes(self):
        crn = self.Mixture.compile_crn()
        return crn


