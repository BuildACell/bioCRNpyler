#  Copyright (c) 2019, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import copy
from typing import List, Union

from .component import Component
from .components_basic import DNA, RNA, Protein
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .mechanism import Mechanism
from .mixture import Mixture
from .parameter import ParameterDatabase
from .reaction import Reaction
from .species import Species


class DNAassembly(DNA):
    """A class that contains a Promoter, RBS, transcript, and protein."""

    def __init__(self, name: str, dna=None, promoter=None, transcript=None,
                 rbs=None, protein=None, length=None,
                 attributes=None, mechanisms=None, parameters=None, initial_conc=None, **keywords):
        """
            Note:
        If transcript is None and protein is not None,
        the DNAassembly will use its transcription mechanisms to produce the protein.
        This is used by Expression Mixtures.
        :param name: name of the DNA assembly
        :param dna:
        :param promoter:
        :param transcript:
        :param rbs:
        :param protein:
        :param length:
        :param attributes:
        :param mechanisms:
        :param parameters:
        :param initial_conc:
        :param keywords: passed into the parent object (DNA)
        """
        self.promoter = None
        self.rbs = None
        self.transcript = None
        self.initial_concentration = initial_conc
        self.name = name
        
        # This has to be called at the end so mechanisms are set for the promoter, RBS, etc.
        DNA.__init__(self, name, length=length, mechanisms=mechanisms,
                     parameters=parameters, initial_conc=initial_conc,
                     attributes=attributes, **keywords)

        self.update_dna(dna, attributes=attributes)
        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_promoter(promoter, transcript=self.transcript, protein=self.protein)
        self.update_rbs(rbs, transcript=self.transcript, protein=self.protein)

    def set_mixture(self, mixture: Mixture) -> None:
        """Set the mixture the Component is in.

        :param mixture: reference to a Mixture instance
        :return: None
        """
        self.mixture = mixture
        if self.promoter is not None:
            self.promoter.set_mixture(mixture)
        if self.rbs is not None:
            self.rbs.set_mixture(mixture)

    def update_dna(self, dna: Union[None, DNA, str], attributes=None) -> None:
        """sets up the dna attribute with a valid DNA instance.

        :param dna: name of a dna sequence or a DNA instance
        :param attributes: Species attribute
        :return: None
        """
        if dna is None:
            self.dna = self.set_species(self.name, material_type="dna", attributes=attributes)
        else:
            self.dna = self.set_species(dna, material_type="dna", attributes=attributes)

        if self.promoter is not None:
            self.promoter.dna = self.dna
        if self.rbs is not None:
            self.rbs.dna = self.dna
        
    def update_transcript(self, transcript: Union[None, RNA, str, bool], attributes=None) -> None:
        """sets up the transcript attribute with a valid RNA instance.

        :param transcript: name of a RNA transcript or RNA instance
        :param attributes: Species attribute
        :return: None
        """
        if transcript is None:
            self.transcript = self.set_species(self.name, material_type="rna", attributes=attributes)

        # this is used for expression mixtures where there is no transcript!
        elif transcript is False:
            self.transcript = None
        else:
            self.transcript = self.set_species(transcript, material_type="rna", attributes=attributes)

        if self.promoter is not None:
            self.promoter.transcript = self.transcript
        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(self, protein: Union[None, Protein, str], attributes=None) -> None:
        """sets up the protein attribute with a valid Protein instance.

        :param protein: name of a protein or Protein instance
        :param attributes: Species attribute
        :return: None
        """
        if protein is None:
            self.protein = self.set_species(self.name, material_type="protein", attributes=attributes)
        else:
            self.protein = self.set_species(protein, material_type="protein", attributes=attributes)

        if self.rbs is not None:
            self.rbs.protein = self.protein
        if self.promoter is not None:
            self.promoter.protein = self.protein

    def update_promoter(self, promoter: Union[Protein, str], transcript: RNA=None, protein: Protein=None) -> None:
        """sets up the promoter attribute with a valid Promoter instance.

        :param promoter: name of a promoter or Promoter instance
        :param transcript: reference to the RNA transcript
        :param protein:
        :return: None
        """
        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(promoter, str):
            self.promoter = Promoter(assembly=self, name=promoter,
                                     transcript=self.transcript,
                                     protein=protein)
        elif isinstance(promoter, Promoter):
            self.promoter = copy.deepcopy(promoter)
            self.promoter.assembly = self
            self.promoter.transcript = self.transcript
            self.promoter.protein = protein
        elif promoter is not None:
            raise ValueError("Improper promoter type received by DNAassembly. "
                             "Expected string or promoter object. "
                             f"Received {repr(promoter)}.")
        else:
            self.promoter = None

        if self.promoter is not None:
            self.promoter.update_parameters(parameter_database=self.parameter_database, overwrite_parameters=False)
            self.promoter.set_mixture(self.mixture)
            self.promoter.add_mechanisms(self.mechanisms, optional_mechanism=True)

    def update_rbs(self, rbs: Union[RBS, str], transcript: RNA=None, protein: Protein=None) -> None:
        """sets up the rbs attribute with a valid RBS instance.

        :param rbs: name of the ribosome binding site or RBS instance
        :param transcript: RNA that contains the ribosome binding site
        :param protein: protein that RNA contains
        :return: None
        """
        if protein is not None:
            self.update_protein(protein)

        if transcript is not None:
            self.update_transcript(transcript)

        if isinstance(rbs, str):
            self.rbs = RBS(assembly=self, name=rbs, protein=self.protein,
                           transcript=self.transcript)
        elif isinstance(rbs, RBS):
            self.rbs = copy.deepcopy(rbs)
            self.rbs.assembly = self
            self.rbs.transcript = self.transcript
            self.rbs.protein = self.protein
        elif rbs is not None:
            raise ValueError(f"Improper rbs type received by DNAassemby. Expected string or RBS object. Received {repr(rbs)}.")
        else:
            self.rbs = None

        if self.rbs is not None:
            self.rbs.update_parameters(parameter_database=self.parameter_database, overwrite_parameters=False)
            self.rbs.set_mixture(self.mixture)
            self.rbs.add_mechanisms(self.mechanisms, optional_mechanism=True)

    def update_species(self) -> List[Species]:
        """collects the list of Species that a DNAassemlby instance holds.

        :return: list of Species that a DNAassemlby instance holds
        """
        species = []
        species.append(self.dna)
        if self.promoter is not None:
            species += self.promoter.update_species()

        if self.rbs is not None:
            species += self.rbs.update_species()

        #deg_mech = self.get_mechanism("rna_degredation", optional_mechanism=True)
        #if deg_mech is not None and self.promoter is not None and self.transcript is not None:
        #    species += deg_mech.update_species(rna=self.transcript, component=self.promoter, part_id=self.transcript.name)

        return species

    def update_reactions(self) -> List[Reaction]:
        """collects the list of Reactions that a DNAassemlby instance holds.

        :return: list of Reactions that a DNAassemlby instance holds.
        """
        reactions = []
        if self.promoter is not None:
            reactions += self.promoter.update_reactions()

        if self.rbs is not None:
            reactions += self.rbs.update_reactions()

        #deg_mech = self.get_mechanism("rna_degredation", optional_mechanism=True)
        #if deg_mech is not None and self.promoter is not None and self.transcript is not None:
        #    reactions += deg_mech.update_reactions(rna=self.transcript, component=self.promoter, part_id=self.transcript.name)

        return reactions

    def update_parameters(self, parameter_file: str=None, parameters:ParameterDatabase=None, overwrite_parameters: bool=True) -> None:
        """updates the parameters stored in dna, promoter and rbs.

        :param parameter_file: valid parameter file
        :param parameters: a parameter database instance
        :param overwrite_parameters: whether to overwrite existing parameters
        :return: None
        """

        DNA.update_parameters(self=self, parameter_file=parameter_file, parameters=parameters, overwrite_parameters=overwrite_parameters)

        if self.promoter is not None:
            self.promoter.update_parameters(parameter_file=parameter_file, parameters=parameters, overwrite_parameters=overwrite_parameters)

        if self.rbs is not None:
            self.rbs.update_parameters(parameter_file=parameter_file, parameters=parameters, overwrite_parameters=overwrite_parameters)

    def add_mechanism(self, mechanism: Mechanism, mech_type: str=None, overwrite: bool=False, optional_mechanism: bool=False) -> None:
        """adds a mechanism of type mech_type to the Component Mechanism dictionary.

        DNA_assembly also adds the mechanisms to its promoter and rbs (but never overwrites them!)

        :param mechanism: reference to a Mechanism instance
        :param mech_type: type of mechanism
        :param overwrite: whether to overwrite the mechanism in Component
        :param optional_mechanism:
        :return: None
        """
        Component.add_mechanism(self, mechanism, mech_type=mech_type, overwrite=overwrite, optional_mechanism=optional_mechanism)

        if self.promoter is not None:
            self.promoter.add_mechanism(mechanism, mech_type=mech_type, optional_mechanism=True)

        if self.rbs is not None:
            self.rbs.add_mechanism(mechanism, mech_type=mech_type, optional_mechanism=True)

    def __str__(self):
        return type(self).__name__ + ": " + self.name

    def __repr__(self):
        txt = str(self)
        if self.promoter is not None:
            txt += "\n\t" + repr(self.promoter)
            txt += "\n\ttranscript = " + repr(self.transcript)
        if self.rbs is not None:
            txt += "\n\t" + repr(self.rbs)
            txt += "\n\tprotein = " + repr(self.protein)
        return txt
