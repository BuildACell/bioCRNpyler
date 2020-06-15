from .mechanism import *
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer
from .mechanisms_enzyme import *


class Transcription_MM(MichalisMentenCopyRXN):
    """Michalis Menten Transcription
        G + RNAP <--> G:RNAP --> G+RNAP+mRNA
    """

    def __init__(self, name="transcription_mm", rnap="RNAP", **keywords):
        if isinstance(rnap, Species):
            self.rnap = rnap
        elif isinstance(rnap, str):
            self.rnap = Species(name=rnap, material_type="protein")
        elif isinstance(rnap, Component) and rnap.get_species() != None:
            self.rnap = rnap.get_species()
        else:
            raise ValueError(
                "'rnap' parameter must be a string or a Component with defined "
                "get_species(), or a chemical_reaction_network.Species object")

        MichalisMentenCopyRXN.__init__(self=self, name=name, enzyme=self.rnap,
                                       mechanism_type="transcription")

    def update_species(self, dna, transcript=None, return_rnap=True,
                       **keywords):
        species = [dna]
        if return_rnap:
            species += [self.rnap]

        species += MichalisMentenCopyRXN.update_species(self, dna)
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")

        species += [transcript]

        return species

    def update_reactions(self, dna, component, part_id = None, complex=None, transcript=None,
                         **keywords):

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        rxns = []
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")
        rxns += MichalisMentenCopyRXN.update_reactions(self, dna, transcript,
                                                       complex=complex, kb=kb,
                                                       ku=ku, kcat=ktx)

        return rxns


class Translation_MM(MichalisMentenCopyRXN):
    """ Michalis Menten Translation
        mRNA + Rib <--> mRNA:Rib --> mRNA + Rib + Protein
    """

    def __init__(self, name="translation_mm", ribosome="Ribo", **keywords):
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        elif isinstance(ribosome, str):
            self.ribosome = Species(name=ribosome, material_type="ribosome")
        elif isinstance(ribosome, Component) and ribosome.get_species() != None:
            self.ribosome = ribosome.get_species()
        else:
            raise ValueError(
                "'ribosome' parameter must be a string, a Component with defined "
                "get_species, or a chemical_reaction_network.species")
        MichalisMentenCopyRXN.__init__(self=self, name=name,
                                       enzyme=self.ribosome,
                                       mechanism_type="translation")

    def update_species(self, transcript, protein=None,
                       return_ribosome=True, **keywords):
        species = [transcript]
        if return_ribosome:
            species += [self.ribosome]
        if protein is None:
            protein = Species(transcript.name, material_type="protein")
        species += [protein]

        species += MichalisMentenCopyRXN.update_species(self, transcript)

        return species

    def update_reactions(self, transcript, component, part_id = None, complex=None,
                         protein=None, **keywords):
        rxns = []

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktl = component.get_parameter("ktl", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        if protein is None:
            protein = Species(transcript.name, material_type="protein")
        rxns += MichalisMentenCopyRXN.update_reactions(self, transcript,
                                                       protein, complex=complex,
                                                       kb=kb, ku=ku,
                                                       kcat=ktl)
        return rxns


class Degredation_mRNA_MM(MichalisMentenRXN):
    """Michalis Menten mRNA Degredation by Endonucleases
       mRNA + Endo <--> mRNA:Endo --> Endo
    """
    def __init__(self, name="rna_degredation_mm", nuclease="RNAase",
                 **keywords):
        if isinstance(nuclease, Species):
            self.nuclease = nuclease
        elif isinstance(nuclease, str):
            self.nuclease = Species(name=nuclease, material_type="protein")
        else:
            raise ValueError("'nuclease' parameter requires a "
                             "chemical_reaction_network.species or a string")
        MichalisMentenRXN.__init__(self=self, name=name, enzyme=self.nuclease,
                                   mechanism_type="rna_degredation")

    def update_species(self, rna, return_nuclease=True, **keywords):
        species = [rna]
        if return_nuclease:
            species += [self.nuclease]
        species += MichalisMentenRXN.update_species(self, rna)
        return species

    def update_reactions(self, rna, component, part_id = None, complex=None, **keywords):

        #Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        kdeg = component.get_parameter("kdeg", part_id = part_id, mechanism = self)
        kb = component.get_parameter("kb", part_id = part_id, mechanism = self)
        ku = component.get_parameter("ku", part_id = part_id, mechanism = self)

        rxns = []
        rxns += MichalisMentenRXN.update_reactions(self, rna, Prod=None, complex=complex, kb=kb, ku=ku, kcat=kdeg)
        return rxns

class SimpleTranscription(Mechanism):
    def __init__(self, name = "simple_transcription", mechanism_type = "transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, transcript = None, **keywords):
        if transcript is None:
            transcript = Species(dna.name, material_type="rna")

        return [dna, transcript]

    def update_reactions(self, dna, component = None, ktx = None, part_id = None, transcript = None, **keywords):

        if ktx == None and Component != None:
            ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        elif component == None and ktx == None:
            raise ValueError("Must pass in component or a value for ktx")

        if transcript is None:
            transcript = Species(dna.name, material_type="rna")

        rxns = [Reaction(inputs = [dna], outputs = [dna, transcript], k = ktx)]
        return rxns

class SimpleTranslation(Mechanism):
    def __init__(self, name = "simple_translation", mechanism_type = "translation"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, transcript, protein = None,  **keywords):
        if protein is None:
            protein = Species(transcript.name, material_type="protein")

        return [transcript, protein]

    def update_reactions(self, transcript, component = None, ktl = None, part_id = None, protein = None, **keywords):

        if ktl == None and Component != None:
            ktl = component.get_parameter("ktl", part_id = part_id, mechanism = self)
        elif component == None and ktl == None:
            raise ValueError("Must pass in component or a value for ktl")

        if protein is None:
            protein = Species(transcript.name, material_type="protein")

        rxns = [Reaction(inputs = [transcript], outputs = [transcript, protein], k = ktl)]
        return rxns

class OneStepGeneExpression(Mechanism):
    def __init__(self, name="gene_expression",
                 mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, protein=None, transcript=None, **keywords):
        species = [dna]
        if protein == None:
            protein = Species(dna.name, material_type="protein")

        species += [protein]
        return species

    def update_reactions(self, dna, component = None, kexpress = None,
                         protein=None, transcript = None, part_id = None, **keywords):

        if kexpress == None and Component != None:
            kexpress = component.get_parameter("kexpress", part_id = part_id, mechanism = self)
        elif component == None and kexpress == None:
            raise ValueError("Must pass in component or a value for kexpress")

        if protein is None:
            protein = Species(dna.name, material_type="protein")
        rxns = [Reaction(inputs=[dna], outputs=[dna, protein], k = kexpress)]
        return rxns


class PositiveHillTranscription(Mechanism):
    #Set the name and mechanism_type
    def __init__(self, name="positivehill_transcription", mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    #Overwrite update_species
    def update_species(self, dna, regulator, transcript = None, **keywords):

        if transcript is None: #Species names can be automatically created
            transcript = Species(dna.name, material_type = "rna")

        return [dna, transcript, regulator] #it is best to return all species that will be involved in the reactions


    #Overwrite update_reactions
    #This always requires the inputs component and part_id to find the relevant parameters
    def update_reactions(self, dna, regulator, component, part_id, transcript = None, **keywords):

        if transcript is None: #Species names should be automatically created the same here as above
            transcript = Species(dna.name, material_type = "rna")

        ktx = component.get_parameter("k", part_id = part_id, mechanism = self)
        n = component.get_parameter("n", part_id = part_id, mechanism = self)
        K = component.get_parameter("K", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)

        params = {"k":ktx, "n":n, "K":K, "s1":regulator, "d":dna}

        reaction = Reaction(inputs = [dna], outputs = [dna, transcript],
                            propensity_type = "proportionalhillpositive", propensity_params = params)

        reaction_leak = Reaction(inputs = [dna], outputs = [dna, transcript], k = kleak)

        #In this case, we just return one reaction
        return [reaction, reaction_leak]

class NegativeHillTranscription(Mechanism):
    #Set the name and mechanism_type
    def __init__(self, name="negativehill_transcription", mechanism_type="transcription"):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    #Overwrite update_species
    def update_species(self, dna, regulator, transcript = None, **keywords):

        if transcript is None: #Species names can be automatically created
            transcript = Species(dna.name, material_type = "rna")

        return [dna, transcript, regulator] #it is best to return all species that will be involved in the reactions

    #Overwrite update_reactions
    #This always requires the inputs component and part_id to find the relevant parameters
    def update_reactions(self, dna, regulator, component, part_id, transcript = None, **keywords):

        if transcript is None: #Species names should be automatically created the same here as above
            transcript = Species(dna.name, material_type = "rna")

        ktx = component.get_parameter("k", part_id = part_id, mechanism = self)
        n = component.get_parameter("n", part_id = part_id, mechanism = self)
        K = component.get_parameter("K", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)

        params = {"k":ktx, "n":n, "K":K, "s1":regulator, "d":dna}

        reaction = Reaction(inputs = [dna], outputs = [dna, transcript],
                            propensity_type = "proportionalhillnegative", propensity_params = params)

        reaction_leak = Reaction(inputs = [dna], outputs = [dna, transcript], k = kleak)

        #In this case, we just return one reaction
        return [reaction, reaction_leak]

class multi_tx(Mechanism):
    '''
    Multi-RNAp Transcription w/ Isomerization:
    Detailed transcription mechanism accounting for each individual
    RNAp occupancy states of gene.

    n ={0, max_occ}
    DNA:RNAp_n + RNAp <--> DNA:RNAp_n_c --> DNA:RNAp_n+1
    DNA:RNAp_n --> DNA:RNAp_0 + n RNAp + n mRNA
    DNA:RNAp_n_c --> DNA:RNAp_0_c + n RNAp + n mRNA

    n --> number of open configuration RNAp on DNA
    max_occ --> Physical maximum number of RNAp on DNA (based on RNAp and DNA dimensions)
    DNA:RNAp_n --> DNA with n open configuration RNAp on it
    DNA:RNAp_n_c --> DNA with n open configuration RNAp and 1 closed configuration RNAp on it

    For more details, see examples/MultiTX_Demo.ipynb
    '''

    # initialize mechanism subclass
    def __init__(self, pol, name='multi_tx', mechanism_type='transcription', **keywords):

        if isinstance(pol,str):
            self.pol = Species(name=pol, material_type='protein')

        elif isinstance(pol,Species):
            self.pol = pol

        else:
            raise ValueError("'pol' must be a string or Species")


        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type, **keywords)

    # species update
    def update_species(self, dna, transcript, component, part_id, **keywords):
        max_occ = int(component.get_parameter("max_occ", part_id = part_id, mechanism = self))
        cp_open = []
        cp_closed = []
        for n in range(1,max_occ + 1):
            name_open = self.pol.name + 'x' + dna.name + '_' + str(n)
            cp_open.append(ComplexSpecies([dna]+[self.pol for i in range(n)],name=name_open))
            if n > 1:
                name_closed = self.pol.name + 'x' + dna.name + '_closed' + '_' + str(n-1)
                cp_closed.append(ComplexSpecies([dna]+[self.pol for i in range(n-1)],name=name_closed))
            else:
                name_closed = self.pol.name + 'x' + dna.name + '_closed' + '_' + str(0)
                cp_closed.append(ComplexSpecies([dna]+[self.pol for i in range(1)],name=name_closed))

        cp_misc = [self.pol,dna,transcript]


        return cp_open + cp_closed + cp_misc

    def update_reactions(self, dna, transcript, component, part_id, **keywords):

        '''
    DNA:RNAp_n + RNAp <--> DNA:RNAp_n_c --> DNA:RNAp_n+1
    kf1 = k1, kr1 = k2, kf2 = k_iso
    DNA:RNAp_n --> DNA:RNAp_0 + n RNAp + n mRNA
    kf = ktx_solo
    DNA:RNAp_n_c --> DNA:RNAp_0_c + n RNAp + n mRNA
    kf = ktx_solo

    max_occ =  maximum occupancy of gene (physical limit)
        '''

        # parameter loading
        k1 = component.get_parameter("k1", part_id = part_id, mechanism = self)
        k2 = component.get_parameter("k2", part_id = part_id, mechanism = self)
        k_iso = component.get_parameter("k_iso", part_id = part_id, mechanism = self)
        ktx_solo = component.get_parameter("ktx_solo", part_id = part_id, mechanism = self)
        max_occ = int(component.get_parameter("max_occ", part_id = part_id, mechanism = self))

        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(1,max_occ + 1):
            name_open = self.pol.name + 'x' + dna.name + '_' + str(n)
            cp_open.append(ComplexSpecies([dna]+[self.pol for i in range(n)],name=name_open))
            if n > 1:
                name_closed = self.pol.name + 'x' + dna.name + '_closed' + '_' + str(n-1)
                cp_closed.append(ComplexSpecies([dna]+[self.pol for i in range(n-1)],name=name_closed))
            else:
                name_closed = self.pol.name + 'x' + dna.name + '_closed' + '_' + str(0)
                cp_closed.append(ComplexSpecies([dna]+[self.pol for i in range(1)],name=name_closed))


        # Reactions
        # polymerase + complex(n) --> complex(n_closed)
        rxn_open_pf = [Reaction(inputs=[self.pol, cp_open[n]], outputs=[cp_closed[n+1]], k=k1) for n in range(0,max_occ-1)]
        rxn_open_pr = [Reaction(inputs=[cp_closed[n+1]], outputs=[self.pol, cp_open[n],], k=k2) for n in range(0,max_occ-1)]

        # isomerization
        rxn_iso = [Reaction(inputs=[cp_closed[n]], outputs=[cp_open[n]], k=k_iso) for n in range(0,max_occ)]

        # release/transcription from open and closed states
        rxn_release_open =  []
        rxn_release_closed = []
        for n in range(0,max_occ):
            rxn_temp1 = Reaction(inputs= [cp_open[n]], outputs=[self.pol for i in range(n+1)] +
                                 [transcript for i in range(n+1)] + [dna], k=ktx_solo)
            rxn_release_open.append(rxn_temp1)

        for n in range(1,max_occ):
            rxn_temp2 = Reaction(inputs= [cp_closed[n]], outputs=[self.pol for i in range(n)] +
                                 [transcript for i in range(n)] + [cp_closed[0]], k=ktx_solo)
            rxn_release_closed.append(rxn_temp2)

        # missing reactions (0 --> 0_closed and v.v. 0_closed --> 0)
        rxn_m1 = Reaction(inputs=[dna,self.pol], outputs=[cp_closed[0]], k=k1)
        rxn_m2 = Reaction(inputs=[cp_closed[0]], outputs=[dna,self.pol], k=k2)

        rxn_all = rxn_open_pf + rxn_open_pr + rxn_iso + rxn_release_open + rxn_release_closed + [rxn_m1, rxn_m2]

        return rxn_all

class multi_tl(Mechanism):
    '''
    Multi-RBZ Translation w/ Isomerization:
    Detailed translation mechanism accounting for each individual
    RBZ occupancy states of mRNA. Still needs some work, so use with caution,
    read all warnings and consult the example notebook.

    n ={0, max_occ}
    mRNA:RBZ_n + RBZ <--> mRNA:RBZ_n_c --> mRNA:RBZ_n+1
    mRNA:RBZ_n --> mRNA:RBZ_0 + n RBZ + n Protein
    mRNA:RBZ_n_c --> mRNA:RBZ_0_c + n RBZ + n Protein

    n --> number of open configuration RBZ on mRNA
    max_occ --> Physical maximum number of RBZ on mRNA (based on RBZ and mRNA dimensions)
    mRNA:RBZ_n --> mRNA with n open configuration RBZ on it
    mRNA:RBZ_n_c --> mRNA with n open configuration RBZ and 1 closed configuration RBZ on it

    For more details, see examples/MultiTX_Demo.ipynb
    '''

    # initialize mechanism subclass
    def __init__(self, ribosome, name='multi_tl', mechanism_type='translation', **keywords):

        if isinstance(ribosome,str):
            self.ribosome = Species(name=ribosome, material_type='protein')

        elif isinstance(ribosome,Species):
            self.ribosome = ribosome

        else:
            raise ValueError("'ribosome' must be a string or Species")

        warn('This mechanism still needs some extra validation, use at your own peril and read the warnings!')
        warn("To properly use this mechanism, set dilution for mRNA-RBZ complexes!")
        warn("I've set RBZ and mRNA-RBZ complexes as protein Species to apply dilution to them, edit if you want something else!")

        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type, **keywords)

    # species update
    def update_species(self, transcript, protein, component, part_id, **keywords):
        max_occ = int(component.get_parameter("max_occ", part_id = part_id, mechanism = self))
        cp_open = []
        cp_closed = []
        for n in range(1,max_occ + 1):
            name_open = self.ribosome.name + 'x' + transcript.name + '_' + str(n)
            cp_open.append(ComplexSpecies([transcript]+[self.ribosome for i in range(n)],name=name_open))

            if n > 1:
                name_closed = self.ribosome.name + 'x' + transcript.name + '_closed' + '_' + str(n-1)
                cp_closed.append(ComplexSpecies([transcript]+[self.ribosome for i in range(n-1)],name=name_closed))
            else:
                name_closed = self.ribosome.name + 'x' + transcript.name + '_closed' + '_' + str(0)
                cp_closed.append(ComplexSpecies([transcript]+[self.ribosome for i in range(1)],name=name_closed))


        cp_misc = [self.ribosome,transcript,protein]

        return cp_open + cp_closed + cp_misc

    def update_reactions(self, transcript, protein, component, part_id, **keywords):
        '''
    mRNA:RBZ_n + RBZ <--> mRNA:RBZ_n_c --> mRNA:RBZ_n+1
    kf1 = kbr, kr1 = kur, kf2 = k_iso_r
    mRNA:RBZ_n --> mRNA:RBZ_0 + n RBZ + n Protein
    kf = ktl_solo
    mRNA:RBZ_n_c --> mRNA:RBZ_0_c + n RBZ + n Protein
    kf = ktl_solo
        '''

        # parameter loading
        kbr = component.get_parameter("kbr", part_id = part_id, mechanism = self)
        kur = component.get_parameter("kur", part_id = part_id, mechanism = self)
        k_iso_r = component.get_parameter("k_iso_r", part_id = part_id, mechanism = self)
        ktl_solo = component.get_parameter("ktl_solo", part_id = part_id, mechanism = self)
        max_occ = int(component.get_parameter("max_occ", part_id = part_id, mechanism = self))


        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(1,max_occ + 1):
            name_open = self.ribosome.name + 'x' + transcript.name + '_' + str(n)
            cp_open.append(ComplexSpecies([transcript]+[self.ribosome for i in range(n)],name=name_open))

            if n > 1:
                name_closed = self.ribosome.name + 'x' + transcript.name + '_closed' + '_' + str(n-1)
                cp_closed.append(ComplexSpecies([transcript]+[self.ribosome for i in range(n-1)],name=name_closed))
            else:
                name_closed = self.ribosome.name + 'x' + transcript.name + '_closed' + '_' + str(0)
                cp_closed.append(ComplexSpecies([transcript]+[self.ribosome for i in range(1)],name=name_closed))

        # Reactions
        # ribosome + complex(n) --> complex(n_closed)
        rxn_open_pf = [Reaction(inputs=[self.ribosome, cp_open[n]], outputs=[cp_closed[n+1]], k=kbr) for n in range(0,max_occ-1)]
        rxn_open_pr = [Reaction(inputs=[cp_closed[n+1]], outputs=[self.ribosome, cp_open[n],], k=kur) for n in range(0,max_occ-1)]

        # isomerization
        rxn_iso = [Reaction(inputs=[cp_closed[n]], outputs=[cp_open[n]], k=k_iso_r) for n in range(0,max_occ)]

        # release/translation from open and closed states
        rxn_release_open =  []
        rxn_release_closed = []
        for n in range(0,max_occ):
            rxn_temp1 = Reaction(inputs= [cp_open[n]], outputs=[self.ribosome for i in range(n+1)] +
                                 [protein for i in range(n+1)] + [transcript], k=ktl_solo)
            rxn_release_open.append(rxn_temp1)

        for n in range(1,max_occ):
            rxn_temp2 = Reaction(inputs= [cp_closed[n]], outputs=[self.ribosome for i in range(n)] +
                                 [protein for i in range(n)] + [cp_closed[0]], k=ktl_solo)
            rxn_release_closed.append(rxn_temp2)

        # missing reactions (0 --> 0_closed and v.v. 0_closed --> 0)
        rxn_m1 = Reaction(inputs=[transcript,self.ribosome], outputs=[cp_closed[0]], k=kbr)
        rxn_m2 = Reaction(inputs=[cp_closed[0]], outputs=[transcript,self.ribosome], k=kur)

        rxn_all = rxn_open_pf + rxn_open_pr + rxn_iso + rxn_release_open + rxn_release_closed + [rxn_m1, rxn_m2]

        return rxn_all
