#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details

from .construct import DNA_part


class Terminator(DNA_part):
    """Transcriptional terminator component for ending transcription.

    A Terminator represents a DNA sequence that signals the end of
    transcription, causing RNA polymerase to dissociate from the DNA
    template and release the newly synthesized RNA transcript. This
    component serves as a structural annotation within genetic constructs
    but does not directly generate species or reactions during CRN
    compilation.

    Parameters
    ----------
    name : str
        Name of the terminator.
    **kwargs
        Additional keyword arguments passed to the parent `DNA_part` class.

    Attributes
    ----------
    name : str
        Name of the terminator.

    See Also
    --------
    Promoter : Component that initiates transcription.
    DNAassembly : Container for terminators and other genetic parts.
    DNA_part : Base class for DNA component parts.

    Notes
    -----
    The Terminator component itself does not generate any species or
    reactions during CRN compilation. It serves primarily as a structural
    element to mark the end of transcription units in genetic constructs.

    Termination behavior, if modeled, would typically be implemented through
    termination efficiency parameters in the transcription mechanism rather
    than through the terminator component itself.

    Examples
    --------
    Create a basic terminator:

    >>> terminator = bcp.Terminator(name='T7_terminator')

    Use a terminator in a DNA construct:

    >>> promoter = bcp.Promoter('ptet')
    >>> rbs = bcp.RBS('RBS_standard')
    >>> cds = bcp.CDS('GFP')
    >>> terminator = bcp.Terminator('BBa_B0022')
    >>> construct = bcp.DNA_construct(
    ...     [promoter, rbs, cds, terminator],
    ...     name='complete_gene'
    ... )

    Create a terminator with custom attributes:

    >>> terminator = bcp.Terminator(
    ...     name='BBa_B0015',
    ...     attributes=['double_terminator']
    ... )

    """

    def __init__(self, name, **kwargs):
        DNA_part.__init__(self, name, **kwargs)
        self.name = name

    def update_species(self):
        """Generate species associated with this terminator.

        Returns
        -------
        list of Species
            Empty list. Terminator components do not generate species
            directly.

        """
        return []

    def update_reactions(self):
        """Generate reactions associated with this terminator.

        Returns
        -------
        list of Reaction
            Empty list. Terminator components do not generate reactions
            directly. Termination behavior is typically modeled through
            transcription mechanism parameters.

        """
        return []
