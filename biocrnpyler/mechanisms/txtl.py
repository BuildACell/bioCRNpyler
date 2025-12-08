from typing import List

from ..core.mechanism import Mechanism
from ..core.propensities import (
    ProportionalHillNegative,
    ProportionalHillPositive,
)
from ..core.reaction import Reaction
from ..core.species import Complex, Species
from ..utils import parameter_to_value
from .enzyme import MichaelisMentenCopy


class OneStepGeneExpression(Mechanism):
    r"""Single-step gene expression mechanism without explicit TX-TL steps.

    A 'transcription' mechanism that models gene expression as a single
    direct reaction from DNA to protein, without explicitly modeling
    transcription and translation as separate steps. This simplified
    mechanism is useful for models where the intermediate mRNA dynamics are
    not important.

    The reaction follows the schema:

        $$ 'G' --> 'G' + 'P' $$

    where G is the gene (DNA) and P is the protein product.

    Parameters
    ----------
    name : str, default='gene_expression'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transcription'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    SimpleTranscription : Explicit transcription mechanism.
    SimpleTranslation : Explicit translation mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is appropriate for modeling scenarios where:

    - mRNA dynamics are much faster than protein dynamics
    - The model focuses on protein-level behavior
    - Computational efficiency is prioritized over mechanistic detail

    The single-step abstraction combines transcription and translation into
    a single effective rate constant. This is often used in coarse-grained
    models of gene regulatory networks.

    Common applications include:

    - Simplified gene regulatory network models
    - Steady-state or quasi-steady-state analyses
    - Systems where mRNA lifetime is negligible
    - High-level circuit design and prototyping

    Required parameters for this mechanism:

    - 'kexpress' : Combined expression rate constant

    Examples
    --------
    Model simple gene expression in a minimal system:

    >>> gene = bcp.DNAassembly(
    ...     name='gfp', promoter='pconst', protein='GFP',
    ... )
    >>> mechanism = bcp.OneStepGeneExpression()
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': mechanism,
    ...     },
    ...     parameters={'kexpress': 0.1}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name='gene_expression', mechanism_type='transcription'
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, transcript=None, protein=None, **kwargs):
        """Generate species for one-step gene expression.

        Returns the DNA and protein species involved in the expression
        reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) that expresses the protein.
        transcript : Species, optional
            Transcript species (unused in this mechanism, accepted for API
            consistency).
        protein : Species
            The protein species produced by expression. If None, no species
            are returned.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [dna, protein] if protein is not None, otherwise
            [dna].

        """
        species = [dna]
        if protein is not None:
            species += [protein]

        return species

    def update_reactions(
        self,
        dna,
        component=None,
        kexpress=None,
        protein=None,
        transcript=None,
        part_id=None,
        **kwargs,
    ):
        r"""Generate reactions for one-step gene expression.

        Creates a single mass-action reaction representing direct gene
        expression from DNA to protein without intermediate transcript.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) that expresses the protein.
        component : Component, optional
            Component containing parameter values. Required if kexpress is
            not provided directly.
        kexpress : Parameter or float, optional
            Expression rate constant. If None, retrieved from component
            parameters.
        protein : Species, optional
            The protein species produced. If None, no reactions are
            generated.
        transcript : Species, optional
            Transcript species (unused in this mechanism, accepted for API
            consistency).
        part_id : str, optional
            Identifier for parameter lookup in the component's parameter
            database.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single mass-action reaction for gene
            expression if protein is not None, otherwise empty list.

        Raises
        ------
        ValueError
            If component is None and kexpress is not provided.

        Notes
        -----
        The reaction has the form:

            $$ 'DNA' --> 'DNA' + 'Protein' \quad ('rate' = 'kexpress') $$

        The DNA is catalytic and appears on both sides of the reaction.

        """
        if kexpress is None and component is not None:
            kexpress = component.get_parameter(
                'kexpress', part_id=part_id, mechanism=self
            )
        elif component is None and kexpress is None:
            raise ValueError("Must pass in component or a value for kexpress")

        if protein is not None:
            return [
                Reaction.from_massaction(
                    inputs=[dna], outputs=[dna, protein], k_forward=kexpress
                )
            ]
        else:
            return []


class SimpleTranscription(Mechanism):
    r"""Simple catalytic transcription mechanism.

    A 'transcription' mechanism that models transcription as a single
    catalytic reaction where DNA directly produces mRNA without explicitly
    modeling RNA polymerase binding and unbinding. This simplified mechanism
    is appropriate when polymerase dynamics are fast compared to other
    processes.

    The reaction follows the schema:

        $$ 'G' --> 'G' + 'T' $$

    where G is the gene (DNA) and T is the transcript (mRNA).

    Parameters
    ----------
    name : str, default='simple_transcription'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transcription'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    Transcription_MM : Explicit Michaelis-Menten transcription.
    OneStepGeneExpression : Combined transcription-translation.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is appropriate for modeling scenarios where:

    - RNA polymerase dynamics are fast and at quasi-steady-state
    - The focus is on transcript-level regulation
    - Computational efficiency is prioritized

    The catalytic abstraction treats the DNA as a template that is not
    consumed in the reaction. In mixtures without explicit transcription
    (e.g., expression mixtures), this mechanism can combine with translation
    rates to produce protein directly.

    Common applications include:

    - Basic transcription models
    - Models where RNAP is not limiting
    - Constitutive or weakly regulated promoters

    Required parameters for this mechanism:

    - 'ktx' : Transcription rate constant
    - 'ktl' : Translation rate constant (when used in expression mixtures
      without explicit transcripts)

    Examples
    --------
    Model constitutive transcription:

    >>> gene = bcp.DNAassembly(
    ...     name='constitutive_GFP',
    ...     promoter='pconst', protein='GFP')
    >>> tx_mechanism = bcp.SimpleTranscription()
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': tx_mechanism,
    ...     },
    ...     parameters={'ktx': 0.05}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name='simple_transcription', mechanism_type='transcription'
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, transcript=None, protein=None, **kwargs):
        """Generate species for simple transcription.

        Returns the DNA and transcript species (and optionally protein)
        involved in the transcription reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) that produces the transcript.
        transcript : Species, optional
            The mRNA transcript species produced. If None, only DNA is
            returned.
        protein : Species, optional
            Protein species (included for API consistency with expression
            mixtures).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing DNA and any non-None transcript/protein species.

        """
        species = [dna]
        if transcript is not None:
            species += [transcript]
        if protein is not None:
            species += [protein]

        return species

    def update_reactions(
        self,
        dna,
        component=None,
        ktx=None,
        part_id=None,
        transcript=None,
        protein=None,
        **kwargs,
    ):
        """Generate reactions for simple transcription.

        Creates a single mass-action reaction representing transcription from
        DNA to mRNA. In expression mixtures without explicit transcription,
        combines transcription and translation rates to produce protein
        directly.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) that produces the transcript.
        component : Component, optional
            Component containing parameter values. Required if ktx is not
            provided directly.
        ktx : Parameter or float, optional
            Transcription rate constant. If None, retrieved from component
            parameters.
        part_id : str, optional
            Identifier for parameter lookup in the component's parameter
            database.
        transcript : Species, optional
            The mRNA transcript species produced. If None and protein is
            not None, produces protein directly.
        protein : Species, optional
            Protein species for expression mixtures without explicit
            transcription.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single mass-action reaction for transcription.

        Raises
        ------
        ValueError
            If component is None and ktx is not provided.

        Notes
        -----
        The reaction has two possible forms:

        1. With transcript: DNA --> DNA + mRNA (rate: 'ktx')
        2. Without transcript: DNA --> DNA + Protein (rate: 'ktx' * 'ktl')

        The second form is used in expression mixtures that skip explicit
        transcript modeling.

        """
        if ktx is None and component is not None:
            ktx = component.get_parameter(
                'ktx', part_id=part_id, mechanism=self
            )
        elif component is None and ktx is None:
            raise ValueError("Must pass in component or a value for ktx")

        # First case only true in Mixtures without transcription (eg
        # Expression Mixtures)
        if transcript is None and protein is not None:
            ktl = component.get_parameter(
                'ktl', part_id=part_id, mechanism=self
            )
            rxns = [
                Reaction.from_massaction(
                    inputs=[dna], outputs=[dna, protein], k_forward=ktx * ktl
                )
            ]
        else:
            rxns = [
                Reaction.from_massaction(
                    inputs=[dna], outputs=[dna, transcript], k_forward=ktx
                )
            ]

        return rxns


class SimpleTranslation(Mechanism):
    r"""Simple catalytic translation mechanism.

    A 'translation' mechanism that models translation as a single catalytic
    reaction where mRNA directly produces protein without explicitly modeling
    ribosome binding and unbinding. This simplified mechanism is appropriate
    when ribosome dynamics are fast compared to other processes.

    The reaction follows the schema:

        $$ 'T' --> 'T' + 'P' $$

    where T is the transcript (mRNA) and P is the protein.

    Parameters
    ----------
    name : str, default='simple_translation'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='translation'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('translation').

    See Also
    --------
    Translation_MM : Explicit Michaelis-Menten translation.
    SimpleTranscription : Simple transcription mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is appropriate for modeling scenarios where:

    - Ribosome dynamics are fast and at quasi-steady-state
    - The focus is on protein-level regulation
    - Computational efficiency is prioritized

    The catalytic abstraction treats the mRNA as a template that is not
    consumed in the reaction. The mRNA persists and can produce multiple
    protein copies.

    Common applications include:

    - Basic translation models
    - Models where ribosomes are not limiting
    - Constitutive protein expression

    Required parameters for this mechanism:

    - 'ktl' : Translation rate constant

    Examples
    --------
    Model simple translation:

    >>> gene = bcp.DNAassembly(
    ...     name='dna_assembly',
    ...     promoter='pconst', rbs='RBS_medium', protein='GFP')
    >>> tl_mechanism = bcp.SimpleTranslation()
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': bcp.SimpleTranscription(),
    ...         'translation': tl_mechanism
    ...     },
    ...     parameter_file='mixtures/pure_parameters.tsv'
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self, name='simple_translation', mechanism_type='translation'
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, transcript, protein=None, **kwargs):
        """Generate species for simple translation.

        Returns the transcript and protein species involved in the
        translation reaction. If protein is None, creates a default protein
        species based on the transcript name.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species that produces protein.
        protein : Species or list of Species, optional
            The protein species produced. If None, a default protein with
            the same name as the transcript is created. Can be a single
            Species or a list.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [transcript] plus protein species (single or
            list).

        """
        if protein is None:
            protein = Species(transcript.name, material_type='protein')
        outlst = [transcript]
        if isinstance(protein, list):
            outlst += protein
        else:
            outlst += [protein]
        return outlst

    def update_reactions(
        self,
        transcript,
        component=None,
        ktl=None,
        part_id=None,
        protein=None,
        **kwargs,
    ):
        r"""Generate reactions for simple translation.

        Creates a single mass-action reaction representing translation from
        mRNA to protein.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species that produces protein.
        component : Component, optional
            Component containing parameter values. Required if ktl is not
            provided directly.
        ktl : Parameter or float, optional
            Translation rate constant. If None, retrieved from component
            parameters.
        part_id : str, optional
            Identifier for parameter lookup in the component's parameter
            database.
        protein : Species, optional
            The protein species produced.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single mass-action reaction for translation if
            transcript is not None, otherwise empty list.

        Raises
        ------
        ValueError
            If component is None and ktl is not provided.

        Notes
        -----
        The reaction has the form:

            $$ 'mRNA' --> 'mRNA' + 'Protein' \quad ('rate' = 'ktl') $$

        The mRNA is catalytic and appears on both sides of the reaction. If
        transcript is None (only occurs in mixtures without transcription),
        no reactions are generated.

        """
        if ktl is None and component is not None:
            ktl = component.get_parameter(
                'ktl', part_id=part_id, mechanism=self
            )
        elif component is None and ktl is None:
            raise ValueError("Must pass in component or a value for ktl")

        # First case only true in Mixtures without transcription (eg
        # Expression Mixtures)
        if transcript is None and protein is not None:
            rxns = []
        else:
            rxns = [
                Reaction.from_massaction(
                    inputs=[transcript],
                    outputs=[transcript, protein],
                    k_forward=ktl,
                )
            ]

        return rxns


class PositiveHillTranscription(Mechanism):
    r"""Transcription regulated by positive Hill function (activation).

    A 'transcription' mechanism that models transcriptional activation using
    a proportional positive Hill function. The transcription rate increases
    with regulator (activator) concentration according to Hill kinetics,
    capturing cooperative binding and activation.

    The reaction follows the schema:

        $$ 'G' --> 'G' + 'T' $$

    with rate:

        $$ 'rate' = k [G] \frac{[R]^n}{K + [R]^n}$$

    where R is the regulator (activator), $n$ is the Hill coefficient, and $K$
    is the activation constant. Optionally includes a basal leak reaction at
    rate 'kleak'.

    Parameters
    ----------
    name : str, default='positivehill_transcription'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transcription'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    NegativeHillTranscription : Transcriptional repression with Hill
        function.
    SimpleTranscription : Simple transcription without regulation.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models transcriptional activation by regulatory proteins
    such as transcription factors. The Hill function captures the sigmoidal
    response typical of cooperative binding, where multiple regulator
    molecules bind to the promoter to activate transcription.

    Key features:

    - Sigmoidal activation response
    - Cooperative binding through Hill coefficient
    - Optional basal transcription (leak)
    - Saturation at high regulator concentrations

    Common applications include:

    - Activatable promoters
    - Positive feedback loops
    - Transcriptional cascades
    - Genetic switches and toggles

    Required parameters for this mechanism:

    - 'k' : Maximum transcription rate constant
    - 'K' : Activation constant (regulator concentration for half-maximal
      activation)
    - 'n' : Hill coefficient (cooperativity)
    - 'kleak' : Basal transcription rate (optional, for leak reaction)

    Examples
    --------
    Model transcriptional activation by an inducer:

    >>> LacI = bcp.Protein('AraC')
    >>> plac = bcp.ActivatablePromoter('pBAD', LacI)
    >>> gene = bcp.DNAassembly(
    ...    name='activiated_GFP',
    ...    promoter=plac, rbs='RBS_medium', protein='GFP')
    >>> tx_mechanism = bcp.PositiveHillTranscription()
    >>> mixture = bcp.Mixture(
    ...     components=[LacI, gene],
    ...     mechanisms={
    ...         'transcription': tx_mechanism,
    ...         'translation': bcp.SimpleTranslation(),
    ...     },
    ...     parameters={'k': 1.0, 'K': 10.0, 'n': 2, 'kleak': 0.01},
    ...     parameter_file='mixtures/pure_parameters.tsv',
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='positivehill_transcription',
        mechanism_type='transcription',
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(
        self,
        dna,
        regulator,
        transcript=None,
        leak=False,
        protein=None,
        **kwargs,
    ):
        """Generate species for positive Hill transcription.

        Returns all species involved in the regulated transcription
        reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (promoter) being transcribed.
        regulator : Species
            The activator species that regulates transcription.
        transcript : Species, optional
            The mRNA transcript species produced. If None and protein is
            provided, protein is produced directly.
        leak : bool, default=False
            If True, includes a leak reaction for basal transcription.
        protein : Species, optional
            Protein species for expression mixtures without explicit
            transcription.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [dna, regulator] plus any non-None
            transcript/protein species.

        """
        species = [dna, regulator]
        if transcript is not None:
            species += [transcript]
        if protein is not None:
            species += [protein]

        return species  # return all species involved in the reactions

    def update_reactions(
        self,
        dna,
        regulator,
        component,
        part_id,
        transcript=None,
        leak=False,
        protein=None,
        **kwargs,
    ):
        r"""Generate reactions for positive Hill transcription.

        Creates regulated transcription reaction(s) using a proportional
        positive Hill function, with optional basal leak reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (promoter) being transcribed.
        regulator : Species
            The activator species that regulates transcription.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup in the component's parameter
            database. Required for parameter lookup.
        transcript : Species, optional
            The mRNA transcript species produced. If None and protein is
            provided, protein is produced directly.
        leak : bool, default=False
            If True, includes a basal leak reaction.
        protein : Species, optional
            Protein species for expression mixtures.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing one or two reactions:

            - Regulated transcription with ProportionalHillPositive
              propensity
            - Optional leak reaction (if leak=True)

        Raises
        ------
        AttributeError
            If component or part_id is None (required for parameter lookup).

        Notes
        -----
        The regulated reaction uses ProportionalHillPositive propensity:
        $$
          'DNA' --> 'DNA' + 'Product' \quad ('rate' = k [G] [R]^n/(K + [R]^n))
        $$
        where Product is either transcript or protein depending on mixture
        type. If `leak=True`, adds:

          $$ 'DNA' --> 'DNA' + 'Product' \quad ('rate' = 'kleak') $$

        """
        ktx = component.get_parameter('k', part_id=part_id, mechanism=self)
        n = component.get_parameter('n', part_id=part_id, mechanism=self)
        K = component.get_parameter('K', part_id=part_id, mechanism=self)
        kleak = component.get_parameter(
            'kleak', part_id=part_id, mechanism=self
        )

        prophill = ProportionalHillPositive(
            k=ktx, K=K, s1=regulator, n=n, d=dna
        )

        reactions = []

        # First case only true in Mixtures without transcription (eg
        # Expression Mixtures)
        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        reactions.append(
            Reaction(
                inputs=[dna],
                outputs=[dna, tx_output],
                propensity_type=prophill,
            )
        )

        if leak:
            reactions.append(
                Reaction.from_massaction(
                    inputs=[dna], outputs=[dna, tx_output], k_forward=kleak
                )
            )

        # In this case, we just return one reaction
        return reactions


class NegativeHillTranscription(Mechanism):
    r"""Transcription regulated by negative Hill function (repression).

    A 'transcription' mechanism that models transcriptional repression using
    a proportional negative Hill function. The transcription rate decreases
    with regulator (repressor) concentration according to Hill kinetics,
    capturing cooperative binding and repression.

    The reaction follows the schema:

        $$ 'G' --> 'G' + 'T' $$

    with rate

        $$ 'rate' = k * [G] * \frac{1}{K + [R]^n} $$

    where R is the regulator (repressor), $n$ is the Hill coefficient, and $K$
    is the repression constant. Optionally includes a basal leak reaction at
    rate 'kleak'.

    Parameters
    ----------
    name : str, default='negativehill_transcription'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transcription'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    PositiveHillTranscription : Transcriptional activation with Hill
        function.
    SimpleTranscription : Simple transcription without regulation.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models transcriptional repression by regulatory proteins
    such as repressor proteins. The Hill function captures the response
    where increasing repressor concentration decreases transcription rate,
    with cooperativity determined by the Hill coefficient.

    Key features:

    - Decreasing transcription with increasing repressor
    - Cooperative repressor binding through Hill coefficient
    - Optional basal transcription (leak)
    - Saturation at high repressor concentrations

    Common applications include:

    - Repressible promoters
    - Negative feedback loops
    - Transcriptional repression cascades
    - Genetic switches and oscillators

    Required parameters for this mechanism:

    - 'k' : Maximum transcription rate constant (when repressor is absent)
    - 'K' : Repression constant (repressor concentration for half-maximal
      repression)
    - 'n' : Hill coefficient (cooperativity)
    - 'kleak' : Basal transcription rate (optional, for leak reaction)

    Examples
    --------
    Model transcriptional repression:

    >>> LacI = bcp.Protein('LacI')
    >>> plac = bcp.RepressiblePromoter('plac', LacI)
    >>> gene = bcp.DNAassembly(
    ...    name='repressed_GFP',
    ...    promoter=plac, rbs='RBS_medium', protein='GFP')
    >>> tx_mechanism = bcp.NegativeHillTranscription()
    >>> mixture = bcp.Mixture(
    ...     components=[LacI, gene],
    ...     mechanisms={
    ...         'transcription': tx_mechanism,
    ...         'translation': bcp.SimpleTranslation(),
    ...     },
    ...     parameters={'k': 1.0, 'K': 10.0, 'n': 2, 'kleak': 0.01},
    ...     parameter_file='mixtures/pure_parameters.tsv',
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        name='negativehill_transcription',
        mechanism_type='transcription',
    ):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(
        self,
        dna,
        regulator,
        transcript=None,
        leak=False,
        protein=None,
        **kwargs,
    ):
        """Generate species for negative Hill transcription.

        Returns all species involved in the regulated transcription
        reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (promoter) being transcribed.
        regulator : Species
            The repressor species that regulates transcription.
        transcript : Species, optional
            The mRNA transcript species produced.
        leak : bool, default=False
            If True, includes a leak reaction for basal transcription.
        protein : Species, optional
            Protein species for expression mixtures.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [dna, regulator] plus any non-None
            transcript/protein species.

        """
        species = [dna, regulator]
        if transcript is not None:
            species += [transcript]
        if protein is not None:
            species += [protein]

        return species  # return all species involved in the reactions

    def update_reactions(
        self,
        dna,
        regulator,
        component,
        part_id,
        transcript=None,
        leak=False,
        protein=None,
        **kwargs,
    ):
        r"""Generate reactions for negative Hill transcription.

        Creates regulated transcription reaction(s) using a proportional
        negative Hill function, with optional basal leak reaction.

        Parameters
        ----------
        dna : Species
            The DNA species (promoter) being transcribed.
        regulator : Species
            The repressor species that regulates transcription.
        component : Component
            Component containing parameter values. Required.
        part_id : str
            Identifier for parameter lookup. Required.
        transcript : Species, optional
            The mRNA transcript species produced.
        leak : bool, default=False
            If True, includes a basal leak reaction.
        protein : Species, optional
            Protein species for expression mixtures.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing one or two reactions: regulated transcription
            with ProportionalHillNegative propensity, and optional leak
            reaction.

        Notes
        -----
        The regulated reaction uses ProportionalHillNegative propensity:
        $$
            'DNA' --> 'DNA' + 'Product' \quad ('rate' = k [G] / (K + [R]^n))
        $$
        If leak=True, adds:
        $$
            'DNA' --> 'DNA' + 'Product' \quad (rate = 'kleak')
        $$

        """
        ktx = component.get_parameter('k', part_id=part_id, mechanism=self)
        n = component.get_parameter('n', part_id=part_id, mechanism=self)
        K = component.get_parameter('K', part_id=part_id, mechanism=self)
        kleak = component.get_parameter(
            'kleak', part_id=part_id, mechanism=self
        )

        prop_hill = ProportionalHillNegative(
            k=ktx, K=K, n=n, s1=regulator, d=dna
        )

        reactions = []

        # First case only true in Mixtures without transcription (eg
        # Expression Mixtures)
        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        reactions.append(
            Reaction(
                inputs=[dna],
                outputs=[dna, tx_output],
                propensity_type=prop_hill,
            )
        )

        if leak:
            reactions.append(
                Reaction.from_massaction(
                    inputs=[dna], outputs=[dna, tx_output], k_forward=kleak
                )
            )

        # In this case, we just return one reaction
        return reactions


class Transcription_MM(MichaelisMentenCopy):
    r"""Michaelis-Menten transcription with explicit RNA polymerase.

    A 'transcription' mechanism that explicitly models RNA polymerase (RNAP)
    binding to DNA, followed by transcription and release. This mechanism
    follows Michaelis-Menten kinetics and is appropriate when RNAP is
    limiting or when explicit modeling of polymerase-DNA interactions is
    needed.

    The reaction follows the schema:

        $$ 'G' + 'RNAP' <--> 'G':'RNAP' --> 'G' + 'RNAP' + 'mRNA' $$

    Parameters
    ----------
    rnap : Species
        RNA polymerase species that catalyzes transcription.
    name : str, default='transcription_mm'
        Name identifier for this mechanism instance.

    Attributes
    ----------
    rnap : Species
        The RNA polymerase species.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    SimpleTranscription : Simple transcription without explicit RNAP.
    Translation_MM : Michaelis-Menten translation.
    MichaelisMentenCopy : Base class for MM copy mechanisms.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models transcription using standard Michaelis-Menten
    kinetics where RNA polymerase acts as an enzyme that binds DNA,
    catalyzes mRNA synthesis, and is released unchanged. This is appropriate
    when:

    - RNAP concentration affects transcription rate
    - Explicit modeling of RNAP-DNA binding is important
    - RNAP sequestration or competition effects are relevant

    The mechanism uses the MichaelisMentenCopy base class which generates:

    1. Reversible binding: DNA + RNAP <--> DNA:RNAP (rates: 'kb', 'ku')
    2. Catalysis: DNA:RNAP --> DNA + RNAP + mRNA (rate: 'ktx')

    Common applications include:

    - TX-TL systems with limited RNAP
    - Models of transcriptional resource allocation
    - Competition between promoters for RNAP
    - Detailed mechanistic models of gene expression

    Required parameters for this mechanism:

    - 'ktx' : Transcription/catalysis rate constant
    - 'kb' : Forward binding rate for RNAP to DNA
    - 'ku' : Reverse unbinding rate for RNAP from DNA

    Examples
    --------
    Model transcription with explicit RNA polymerase:

    >>> gene = bcp.DNAassembly(
    ...     name='gene_assembly',
    ...     promoter='pconst', transcript='mRNA')
    >>> mechanism = bcp.Transcription_MM(rnap=rnap)
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={'transcription': mechanism},
    ...     parameters={'ktx': 0.05, 'kb': 1.0, 'ku': 0.1}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(self, rnap: Species, name='transcription_mm', **kwargs):
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' parameter must be a Species.")

        MichaelisMentenCopy.__init__(
            self=self, name=name, mechanism_type='transcription'
        )

    def update_species(self, dna, transcript=None, protein=None, **kwargs):
        """Generate species for Michaelis-Menten transcription.

        Creates species involved in RNAP-mediated transcription including DNA,
        RNAP, the DNA:RNAP complex, and transcript or protein product.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        transcript : Species, optional
            The mRNA transcript species produced. If None and protein is
            provided, protein is produced directly (for expression mixtures).
        protein : Species, optional
            Protein species for expression mixtures without explicit
            transcription.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [dna, rnap, dna:rnap complex, product] where
            product is either transcript or protein.

        """
        species = [dna]

        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        species += MichaelisMentenCopy.update_species(
            self, enzyme=self.rnap, substrate=dna, product=tx_output
        )

        return species

    def update_reactions(
        self,
        dna,
        component,
        part_id=None,
        complex=None,
        transcript=None,
        protein=None,
        **kwargs,
    ):
        """Generate reactions for Michaelis-Menten transcription.

        Creates Michaelis-Menten transcription reactions including reversible
        RNAP-DNA binding and catalytic transcript production.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Complex, optional
            Pre-specified DNA:RNAP complex species. If None, automatically
            created.
        transcript : Species, optional
            The mRNA transcript species produced. If None and protein is
            provided, protein is produced directly.
        protein : Species, optional
            Protein species for expression mixtures.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of two reactions:

            - Reversible RNAP-DNA binding (rates: 'kb', 'ku')
            - Catalytic transcription (rate: 'ktx')

        Notes
        -----
        The reactions follow the Michaelis-Menten scheme:

        1. DNA + RNAP <--> DNA:RNAP (rates: 'kb' and 'ku')
        2. DNA:RNAP --> DNA + RNAP + Product (rate: 'ktx')

        Where Product is either transcript or protein depending on mixture
        type.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        ktx = component.get_parameter('ktx', part_id=part_id, mechanism=self)
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)

        rxns = []

        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        rxns += MichaelisMentenCopy.update_reactions(
            self,
            enzyme=self.rnap,
            substrate=dna,
            product=tx_output,
            complex=complex,
            kb=kb,
            ku=ku,
            kcat=ktx,
        )

        return rxns


class Translation_MM(MichaelisMentenCopy):
    r"""Michaelis-Menten translation with explicit ribosome.

    A 'translation' mechanism that explicitly models ribosome binding to
    mRNA, followed by translation and release. This mechanism follows
    Michaelis-Menten kinetics and is appropriate when ribosomes are limiting
    or when explicit modeling of ribosome-mRNA interactions is needed.

    The reaction follows the schema:
    $$
        'mRNA' + 'Ribo' <--> 'mRNA':'Ribo' --> 'mRNA' + 'Ribo' + 'Protein'
    $$

    Parameters
    ----------
    ribosome : Species
        Ribosome species that catalyzes translation.
    name : str, default='translation_mm'
        Name identifier for this mechanism instance.

    Attributes
    ----------
    ribosome : Species
        The ribosome species.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('translation').

    See Also
    --------
    SimpleTranslation : Simple translation without explicit ribosomes.
    Transcription_MM : Michaelis-Menten transcription.
    MichaelisMentenCopy : Base class for MM copy mechanisms.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism models translation using standard Michaelis-Menten
    kinetics where ribosomes act as enzymes that bind mRNA, catalyze
    protein synthesis, and are released unchanged. This is appropriate when:

    - Ribosome concentration affects translation rate
    - Explicit modeling of ribosome-mRNA binding is important
    - Ribosome sequestration or competition effects are relevant

    The mechanism uses the MichaelisMentenCopy base class which generates:

    1. Reversible binding: mRNA + Rib <--> mRNA:Rib (rates: 'kb', 'ku')
    2. Catalysis: mRNA:Rib --> mRNA + Rib + Protein (rate: 'ktl')

    Common applications include:

    - TX-TL systems with limited ribosomes
    - Models of translational resource allocation
    - Competition between mRNAs for ribosomes
    - Detailed mechanistic models of protein expression

    Required parameters for this mechanism:

    - 'ktl' : Translation/catalysis rate constant
    - 'kb' : Forward binding rate for ribosome to mRNA
    - 'ku' : Reverse unbinding rate for ribosome from mRNA

    Examples
    --------
    Model translation with explicit ribosomes:

    >>> gene = bcp.DNAassembly(
    ...     name='gene_assembly',
    ...     promoter='pconst', transcript='mRNA',
    ...     rbs='RBS_medium', protein='GFP')
    >>> ribosome = bcp.Species('Ribosome', material_type='protein')
    >>> tl_mechanism = bcp.Translation_MM(ribosome=ribosome)
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': bcp.SimpleTranscription(),
    ...         'translation': mechanism},
    ...     parameters={'ktx': 0.05, 'ktl': 0.1, 'kb': 1.0, 'ku': 0.1}
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(self, ribosome: Species, name='translation_mm', **kwargs):
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("ribosome must be a Species!")
        MichaelisMentenCopy.__init__(
            self=self, name=name, mechanism_type='translation'
        )

    def update_species(self, transcript, protein, **kwargs):
        """Generate species for Michaelis-Menten translation.

        Creates species involved in ribosome-mediated translation including
        mRNA, ribosome, the mRNA:ribosome complex, and protein product.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated. Can be None in
            expression mixtures without explicit transcription.
        protein : Species
            The protein species produced by translation.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing translation species. In expression mixtures
            without transcription (transcript is None), returns [protein].
            Otherwise returns [ribosome, transcript, mRNA:ribosome complex,
            protein].

        """
        species = []

        # This can only occur in expression mixtures
        if transcript is None and protein is not None:
            species += Species.flatten_list([protein])
        else:
            species += MichaelisMentenCopy.update_species(
                self,
                enzyme=self.ribosome,
                substrate=transcript,
                product=protein,
            )

        return species

    def update_reactions(
        self,
        transcript,
        protein,
        component,
        part_id=None,
        complex=None,
        **kwargs,
    ):
        """Generate reactions for Michaelis-Menten translation.

        Creates Michaelis-Menten translation reactions including reversible
        ribosome-mRNA binding and catalytic protein production.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated. Can be None in
            expression mixtures.
        protein : Species
            The protein species produced by translation.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Complex, optional
            Pre-specified mRNA:ribosome complex species. If None,
            automatically created.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of two reactions for translation if transcript is not None:

            - Reversible ribosome-mRNA binding (rates: 'kb', 'ku')
            - Catalytic translation (rate: 'ktl')

            Returns empty list if transcript is None (expression mixtures).

        Notes
        -----
        The reactions follow the Michaelis-Menten scheme:

        1. mRNA + Ribosome <--> mRNA:Ribosome (rates: 'kb' and 'ku')
        2. mRNA:Ribosome --> mRNA + Ribosome + Protein (rate: 'ktl')

        In expression mixtures without explicit transcription, no translation
        reactions are generated (empty list returned).

        """
        rxns = []

        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        ktl = component.get_parameter('ktl', part_id=part_id, mechanism=self)
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)

        # This can only occur in expression mixtures
        if transcript is None and protein is not None:
            pass
        else:
            rxns += MichaelisMentenCopy.update_reactions(
                self,
                enzyme=self.ribosome,
                substrate=transcript,
                product=protein,
                complex=complex,
                kb=kb,
                ku=ku,
                kcat=ktl,
            )
        return rxns


class Energy_Transcription_MM(Mechanism):
    r"""Michaelis-Menten transcription with explicit energy consumption.

    A 'transcription' mechanism that models transcription with explicit
    consumption of energy sources (fuel species like NTPs) and production of
    waste products. This mechanism couples RNAP-DNA binding with
    length-dependent fuel consumption to model realistic transcription
    energetics.

    The reaction follows the schema:
    $$
        & 'G' + 'RNAP' <--> 'G':'RNAP' \\
        & 'Fuel' + 'G':'RNAP' --> 'G' + 'RNAP' + 'T' \\
        & 'Fuel' + 'G':'RNAP' --> 'G':'RNAP' + 'wastes'
    $$

    Transcription occurs at rate 'ktx' / L (length-dependent), while fuel
    consumption occurs at rate 'ktx', resulting in L times more fuel
    consumption than transcripts produced.

    Parameters
    ----------
    rnap : Species
        RNA polymerase species that catalyzes transcription.
    fuels : list of Species
        List of fuel species (e.g., NTPs) consumed during transcription.
    wastes : list of Species
        List of waste species (e.g., pyrophosphate) produced during
        transcription.
    name : str, default='energy_transcription_mm'
        Name identifier for this mechanism instance.

    Attributes
    ----------
    rnap : Species
        The RNA polymerase species.
    fuels : list of Species
        Fuel species consumed during transcription.
    wastes : list of Species
        Waste species produced during transcription.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    Transcription_MM : MM transcription without explicit energy.
    MichaelisMentenCopy : Base class for MM copy mechanisms.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism provides a more realistic model of transcription by
    explicitly tracking energy consumption. Key features:

    - Length-dependent transcription rate ('ktx' / L)
    - Explicit fuel (NTP) consumption
    - Waste product generation
    - RNAP-DNA binding dynamics

    The length parameter L represents the gene length in appropriate units,
    and the mechanism automatically scales fuel consumption to match the
    energetic requirements of synthesizing an mRNA of length L.

    Common applications include:

    - Detailed TX-TL models with explicit resources
    - Models of transcriptional burden and resource depletion
    - Systems where NTP availability affects gene expression

    Required parameters for this mechanism:

    - 'ktx' : Base transcription rate constant
    - 'kb' : Forward binding rate for RNAP to DNA
    - 'ku' : Reverse unbinding rate for RNAP from DNA
    - 'length' : Gene length (for length-dependent transcription)

    Examples
    --------
    Model transcription with explicit NTP consumption:

    >>> gene = bcp.DNAassembly(
    ...     name='constitutive_promoter',
    ...     promoter='pconst', rbs='RBS_medium', protein='GFP')
    >>> rnap = bcp.Species('RNAP')
    >>> ntp = bcp.Species('NTP')
    >>> ppi = bcp.Species('PPi')
    >>> mechanism = bcp.Energy_Transcription_MM(
    ...     rnap=rnap,
    ...     fuels=[ntp],
    ...     wastes=[ppi]
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': mechanism,
    ...         'translation': bcp.SimpleTranslation()
    ...     },
    ...     parameters={'ktx': 0.05, 'kb': 1.0, 'ku': 0.1, 'length': 1000},
    ...     parameter_file='mixtures/pure_parameters.tsv'
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        rnap: Species,
        fuels: List[Species],
        wastes: List[Species],
        name='energy_transcription_mm',
        **kwargs,
    ):
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' parameter must be a Species.")

        if all([isinstance(s, Species) for s in fuels]):
            self.fuels = fuels
        else:
            raise ValueError(
                "recieved a non-Species object in the fuels list."
            )

        if all([isinstance(s, Species) for s in fuels]):
            self.wastes = wastes
        else:
            raise ValueError("wastes must be a list of Species!")

        Mechanism.__init__(
            self=self, name=name, mechanism_type='transcription'
        )

    def update_species(self, dna, transcript=None, protein=None, **kwargs):
        """Generate species for energy-consuming transcription.

        Creates species involved in transcription with explicit energy
        consumption including DNA, RNAP, the DNA:RNAP complex, transcript,
        and fuel species.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        transcript : Species, optional
            The mRNA transcript species produced.
        protein : Species, optional
            Protein species (unused in this mechanism, accepted for API
            consistency).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [dna, rnap, transcript, fuels..., dna:rnap
            complex].

        """
        species = [dna, self.rnap, transcript] + self.fuels
        bound_complex = Complex([dna, self.rnap])
        species += [bound_complex]

        return species

    def update_reactions(
        self,
        dna,
        component,
        part_id=None,
        complex=None,
        transcript=None,
        protein=None,
        **kwargs,
    ):
        """Generate reactions for energy-consuming transcription.

        Creates three reactions modeling transcription with explicit fuel
        consumption and waste production: RNAP-DNA binding, length-dependent
        transcription, and fuel consumption.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Complex, optional
            Pre-specified DNA:RNAP complex species (unused, complex is
            created internally).
        transcript : Species, optional
            The mRNA transcript species produced.
        protein : Species, optional
            Protein species (unused, accepted for API consistency).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of three reactions:

            - RNAP-DNA binding (rates: 'kb', 'ku')
            - Transcription with fuel (rate: 'ktx' / 'length')
            - Fuel consumption producing wastes (rate: 'ktx')

        Notes
        -----
        The reactions model transcription energetics:

        1. DNA + RNAP <--> DNA:RNAP (rates: 'kb' and 'ku')
        2. Fuel + DNA:RNAP --> Fuel + DNA + RNAP + mRNA (rate: 'ktx' / L)
        3. Fuel + DNA:RNAP --> DNA:RNAP + wastes (rate: 'ktx')

        The length-dependent transcription rate ('ktx' / L) ensures that L
        times more fuel is consumed than transcripts produced, reflecting the
        energetic cost of synthesizing a transcript of length L.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        ktx = component.get_parameter('ktx', part_id=part_id, mechanism=self)
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)
        L = component.get_parameter('length', part_id=part_id, mechanism=self)

        bound_complex = Complex([dna, self.rnap])

        # RNAP DNA Binding
        r1 = Reaction.from_massaction(
            [dna, self.rnap], [bound_complex], k_forward=kb, k_reverse=ku
        )
        # Transcription
        r2 = Reaction.from_massaction(
            self.fuels + [bound_complex],
            self.fuels + [dna, self.rnap, transcript],
            k_forward=parameter_to_value(ktx.value) / parameter_to_value(L),
        )
        # Fuel consumption
        r3 = Reaction.from_massaction(
            self.fuels + [bound_complex],
            [bound_complex] + self.wastes,
            k_forward=ktx,
        )

        return [r1, r2, r3]


class Energy_Translation_MM(Mechanism):
    r"""Michaelis-Menten translation with explicit energy consumption.

    A 'translation' mechanism that models translation with explicit
    consumption of energy sources (fuel species like amino acids/NTPs) and
    production of waste products. This mechanism couples ribosome-mRNA binding
    with length-dependent fuel consumption to model realistic translation
    energetics.

    The reaction follows the schema:
    $$
        & 'mRNA' + 'Ribo' <--> 'mRNA':'Ribo' \\
        & 'Fuel' + 'mRNA':'Ribo' --> 'mRNA' + 'Ribo' + 'Protein' + 'Fuel' \\
        & 'Fuel' + 'mRNA':'Ribo' --> 'mRNA':'Ribo' + 'wastes'
    $$
    Translation occurs at rate 'ktl' / L (length-dependent), while fuel
    consumption occurs at rate 'ktl', resulting in L times more fuel
    consumption than proteins produced.

    Parameters
    ----------
    ribosome : Species
        Ribosome species that catalyzes translation.
    fuels : list of Species
        List of fuel species (e.g., amino acids, GTP) consumed during
        translation.
    wastes : list of Species
        List of waste species produced during translation.
    name : str, default='energy_translation_mm'
        Name identifier for this mechanism instance.

    Attributes
    ----------
    ribosome : Species
        The ribosome species.
    fuels : list of Species
        Fuel species consumed during translation.
    wastes : list of Species
        Waste species produced during translation.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('translation').

    See Also
    --------
    Translation_MM : MM translation without explicit energy.
    Energy_Transcription_MM : MM transcription with explicit energy.
    MichaelisMentenCopy : Base class for MM copy mechanisms.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism provides a more realistic model of translation by
    explicitly tracking energy consumption. Key features:

    - Length-dependent translation rate ('ktl' / L)
    - Explicit fuel (amino acid/NTP) consumption
    - Waste product generation
    - Ribosome-mRNA binding dynamics

    The length parameter L represents the gene/protein length in appropriate
    units, and the mechanism automatically scales fuel consumption to match
    the energetic requirements of synthesizing a protein of length L.

    Common applications include:

    - Detailed TX-TL models with explicit resources
    - Models of translational burden and resource depletion
    - Systems where amino acid availability affects protein expression

    Required parameters for this mechanism:

    - 'ktl' : Base translation rate constant
    - 'kb' : Forward binding rate for ribosome to mRNA
    - 'ku' : Reverse unbinding rate for ribosome from mRNA
    - 'length' : Gene/protein length (for length-dependent translation)

    Examples
    --------
    Model translation with explicit amino acid consumption:

    >>> gene = bcp.DNAassembly(
    ...    name='constitutive_promoter',
    ...    promoter='pconst', rbs='RBS_medium', protein='GFP')
    >>> ribosome = bcp.Species('Ribosome')
    >>> amino_acids = bcp.Species('AA')
    >>> waste = bcp.Species('waste')
    >>> mechanism = bcp.Energy_Translation_MM(
    ...     ribosome=ribosome,
    ...     fuels=[amino_acids],
    ...     wastes=[waste]
    ... )
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': bcp.SimpleTranscription(),
    ...         'translation': mechanism,
    ...     },
    ...     parameters={'ktl': 0.1, 'kb': 1.0, 'ku': 0.1, 'length': 300},
    ...     parameter_file='mixtures/pure_parameters.tsv'
    ... )
    >>> mixture.compile_crn()

    """

    def __init__(
        self,
        ribosome: Species,
        fuels: List[Species],
        wastes=List[Species],
        name='energy_translation_mm',
        **kwargs,
    ):
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("ribosome must be a Species!")

        if all([isinstance(s, Species) for s in fuels]):
            self.fuels = fuels
        else:
            raise ValueError("Fuels must be a list of Species!")

        if all([isinstance(s, Species) for s in fuels]):
            self.wastes = wastes
        else:
            raise ValueError("wastes must be a list of Species!")

        Mechanism.__init__(self=self, name=name, mechanism_type='translation')

    def update_species(self, transcript, protein, **kwargs):
        """Generate species for energy-consuming translation.

        Creates species involved in translation with explicit energy
        consumption including mRNA, ribosome, the mRNA:ribosome complex,
        protein, and fuel species.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated.
        protein : Species
            The protein species produced by translation.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [fuels..., ribosome, protein, mRNA:ribosome
            complex].

        """
        species = self.fuels + [self.ribosome, protein]
        bound_complex = Complex([transcript, self.ribosome])
        species += [bound_complex]

        return species

    def update_reactions(
        self,
        transcript,
        protein,
        component,
        part_id=None,
        complex=None,
        **kwargs,
    ):
        """Generate reactions for energy-consuming translation.

        Creates three reactions modeling translation with explicit fuel
        consumption and waste production: ribosome-mRNA binding,
        length-dependent translation, and fuel consumption.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated.
        protein : Species
            The protein species produced by translation.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            component.name.
        complex : Complex, optional
            Pre-specified mRNA:ribosome complex species (unused, complex is
            created internally).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of three reactions:

            - Ribosome-mRNA binding (rates: 'kb', 'ku')
            - Translation with fuel (rate: 'ktl' / 'length')
            - Fuel consumption producing wastes (rate: 'ktl')

        Notes
        -----
        The reactions model translation energetics:

        1. mRNA + Ribosome <--> mRNA:Ribosome (rates: 'kb' and 'ku')
        2. Fuel + mRNA:Ribosome --> Fuel + mRNA + Ribosome + Protein (rate:
           'ktl' / L)
        3. Fuel + mRNA:Ribosome --> mRNA:Ribosome + Wastes (rate: 'ktl')

        The length-dependent translation rate ('ktl' / L) ensures that L
        times more fuel is consumed than proteins produced, reflecting the
        energetic cost of synthesizing a protein of length L.

        """
        # Get Parameters
        if part_id is None and component is not None:
            part_id = component.name

        ktl = component.get_parameter('ktl', part_id=part_id, mechanism=self)
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)
        L = component.get_parameter('length', part_id=part_id, mechanism=self)

        bound_complex = Complex([transcript, self.ribosome])

        # RNAP DNA Binding
        r1 = Reaction.from_massaction(
            [transcript, self.ribosome],
            [bound_complex],
            k_forward=kb,
            k_reverse=ku,
        )
        # Transcription
        r2 = Reaction.from_massaction(
            self.fuels + [bound_complex],
            self.fuels + [transcript, self.ribosome, protein],
            k_forward=parameter_to_value(ktl.value) / parameter_to_value(L),
        )
        # Fuel consumption
        r3 = Reaction.from_massaction(
            self.fuels + [bound_complex],
            [bound_complex] + self.wastes,
            k_forward=ktl,
        )

        return [r1, r2, r3]


class multi_tx(Mechanism):
    r"""Multi-polymerase transcription with isomerization and occupancy.

    A 'transcription' mechanism that explicitly models multiple RNA
    polymerases (RNAPs) binding to a single gene simultaneously, accounting
    for polymerase spacing, isomerization between open and closed
    configurations, and competitive binding. This detailed mechanism captures
    transcriptional queueing and interference effects.

    The reaction scheme follows:
    $$
        & 'DNA':'RNAP'_n + 'RNAP' <--> 'DNA':'RNAP'_{n_c}
            --> 'DNA':'RNAP'_{n+1} \\
        & 'DNA':'RNAP'_n --> 'DNA':'RNAP'_0 + n 'RNAP' + n 'mRNA' \\
        & 'DNA':'RNAP'_{n_c} --> 'DNA':'RNAP'_{0_c} + n 'RNAP' + n 'mRNA'
    $$
    where:

    - $n$ = {0, 1, ..., max_occ} is the number of polymerases
    - max_occ is the physical maximum based on RNAP and DNA dimensions
    - $'DNA':'RNAP'_n$ represents DNA with $n$ open configuration RNAPs
    - $'DNA':'RNAP'_{n_c}$ represents DNA with $n$ open RNAPs and 1 closed
      RNAP

    Parameters
    ----------
    pol : Species
        RNA polymerase species.
    name : str, default='multi_tx'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='transcription'
        Type classification of this mechanism.

    Attributes
    ----------
    pol : Species
        The RNA polymerase species.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('transcription').

    See Also
    --------
    Transcription_MM : Simple Michaelis-Menten transcription.
    multi_tl : Multi-ribosome translation mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism provides a detailed, spatially-aware model of
    transcription that captures:

    - Multiple polymerases on a single gene
    - Polymerase queueing and spacing constraints
    - Open/closed conformational states (isomerization)
    - Coordinated transcript release from multiple polymerases

    The model is appropriate for detailed mechanistic studies where:

    - Polymerase density and spacing matter
    - Transcriptional queueing affects expression
    - Multiple concurrent transcription events are important

    The mechanism generates many species (2 * max_occ complexes) and
    reactions (O(max_occ)), so it should be used with caution for
    computational efficiency.

    Required parameters for this mechanism:

    - 'ktx' : Transcription/release rate constant
    - 'kb' : Forward binding rate for polymerase to DNA
    - 'ku' : Reverse unbinding rate for polymerase from DNA
    - 'k_iso' : Isomerization rate from closed to open configuration
    - 'max_occ' : Maximum polymerase occupancy (integer)

    Examples
    --------
    Model multi-polymerase transcription:

    >>> gene = bcp.DNAassembly(
    ...     name='dna_assembly',
    ...     promoter='pconst', rbs='RBS_medium', protein='GFP')
    >>> rnap = bcp.Species('RNAP')
    >>> tx_mechanism = bcp.multi_tx(pol=rnap)
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': tx_mechanism,
    ...         'translation': bcp.SimpleTranslation()
    ...     },
    ...     parameters={
    ...         'ktx': 0.05, 'ktl': 0.1, 'kb': 1.0, 'ku': 0.1,
    ...         'k_iso': 0.5, 'max_occ': 5
    ...     }
    ... )
    >>> mixture.compile_crn()

    For more details, see examples/MultiTX_Demo.ipynb.

    """

    def __init__(
        self,
        pol: Species,
        name: str = 'multi_tx',
        mechanism_type: str = 'transcription',
        **kwargs,
    ):
        if isinstance(pol, Species):
            self.pol = pol
        else:
            raise ValueError("'pol' must be a Species")

        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    # species update
    def update_species(
        self, dna, transcript, component, part_id, protein=None, **kwargs
    ):
        """Generate species for multi-polymerase transcription.

        Creates all species and complexes involved in multi-polymerase
        transcription including DNA with varying numbers of bound polymerases
        in both open and closed configurations.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        transcript : Species
            The mRNA transcript species produced.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup. Required to retrieve 'max_occ'
            parameter.
        protein : Species, optional
            Protein species (unused, accepted for API consistency).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing:

            - DNA:RNAP_n complexes in open configuration (n = 0 to max_occ-1)
            - DNA:RNAP_n complexes in closed configuration (n = 0 to
              max_occ-1)
            - Individual species [polymerase, dna, transcript]

        Notes
        -----
        For each occupancy level $n$ from 0 to max_occ-1, two complex species
        are created:

        - Open configuration: DNA with $n+1$ polymerases, all in open state
        - Closed configuration: DNA with $n+1$ polymerases (n open, 1 closed)

        The 'max_occ' parameter determines the maximum number of polymerases
        that can occupy the gene simultaneously, typically based on physical
        spacing constraints.

        """
        max_occ = int(
            component.get_parameter(
                'max_occ',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        )
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex(
                    [dna] + [self.pol for i in range(n + 1)],
                    attributes=['open'],
                )
            )
            cp_closed.append(
                Complex(
                    [dna] + [self.pol for i in range(n + 1)],
                    attributes=['closed'],
                )
            )

        cp_misc = [self.pol, dna, transcript]

        return cp_open + cp_closed + cp_misc

    def update_reactions(
        self, dna, transcript, component, part_id, protein=None, **kwargs
    ):
        """Generate reactions for multi-polymerase transcription.

        Creates reactions modeling multiple polymerases binding to DNA,
        isomerization between configurations, and transcript release with
        coordination among multiple polymerases.

        Parameters
        ----------
        dna : Species
            The DNA species (gene) being transcribed.
        transcript : Species
            The mRNA transcript species produced.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup. Required to retrieve parameters.
        protein : Species, optional
            Protein species (unused, accepted for API consistency).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of reactions including:

            - Polymerase binding reactions
              (DNA:RNAP$_n$ + RNAP <--> DNA:RNAP$_{n_c}$)
            - Isomerization reactions (DNA:RNAP$_{n_c}$ --> $DNA:RNAP$_n$)
            - Release reactions from open states
              (DNA:RNAP$_n$ --> DNA + $n$ RNAP + $n$ mRNA)
            - Release reactions from closed states
              (DNA:RNAP$_{n_c}$ --> DNA:RNAP$_{0_c}$ + $n$ RNAP + $n$ mRNA)
            - Base binding reaction (DNA + RNAP <--> DNA:RNAP$_{0_c}$)

        Notes
        -----
        The reaction scheme captures:

        1. DNA:RNAP$_n$ + RNAP <--> DNA:RNAP$_{(n+1)_c}$ (rates: 'kb', 'ku')
        2. DNA:RNAP$_{n_c}$ --> DNA:RNAP$_n$ (rate: 'k_iso')
        3. DNA:RNAP$_n$ --> DNA + $n$ RNAP + $n$ mRNA (rate: 'ktx')
        4. DNA:RNAP$_{n_c}$ --> DNA:RNAP$_{0_c}$ + $n$ RNAP + $n$ mRNA
           (rate: 'ktx')

        Where $n$ ranges from 0 to max_occ-1. The 'max_occ' parameter
        represents the physical maximum occupancy of polymerases on the gene.

        """
        # parameter loading
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)
        k_iso = component.get_parameter(
            'k_iso', part_id=part_id, mechanism=self
        )
        ktx = component.get_parameter('ktx', part_id=part_id, mechanism=self)
        max_occ = int(
            component.get_parameter(
                'max_occ',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        )

        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            # has n polymerases all open
            cp_open.append(
                Complex(
                    [dna] + [self.pol for i in range(n + 1)],
                    attributes=['open'],
                )
            )
            # has n-1 open polymerases and 1 closed polymerase
            cp_closed.append(
                Complex(
                    [dna] + [self.pol for i in range(n + 1)],
                    attributes=['closed'],
                )
            )

        # Reactions
        # polymerase + complex(n) <--> complex(n+1)_closed
        rxn_open_p = [
            Reaction.from_massaction(
                inputs=[self.pol, cp_open[n]],
                outputs=[cp_closed[n + 1]],
                k_forward=kb,
                k_reverse=ku,
            )
            for n in range(0, max_occ - 1)
        ]
        # rxn_open_pr = [
        #    Reaction.from_massaction(
        #        inputs=[cp_closed[n + 1]], outputs=[self.pol, cp_open[n], ],
        #        k_forward=k2)
        #    for n in range(0, max_occ - 1)
        #    ]

        # isomerization
        # complex(n)_closes --> complex(n)
        rxn_iso = [
            Reaction.from_massaction(
                inputs=[cp_closed[n]], outputs=[cp_open[n]], k_forward=k_iso
            )
            for n in range(0, max_occ)
        ]

        # release/transcription from open and closed states
        rxn_release_open = []
        rxn_release_closed = []
        for n in range(0, max_occ):
            rxn_temp1 = Reaction.from_massaction(
                inputs=[cp_open[n]],
                outputs=[self.pol for i in range(n + 1)]
                + [transcript for i in range(n + 1)]
                + [dna],
                k_forward=ktx,
            )
            rxn_release_open.append(rxn_temp1)

        for n in range(1, max_occ):
            rxn_temp2 = Reaction.from_massaction(
                inputs=[cp_closed[n]],
                outputs=[self.pol for i in range(n)]
                + [transcript for i in range(n)]
                + [cp_closed[0]],
                k_forward=ktx,
            )
            rxn_release_closed.append(rxn_temp2)

        # base case pol + dna <--> complex(n=1)_open
        rxn_m1 = Reaction.from_massaction(
            inputs=[dna, self.pol],
            outputs=[cp_closed[0]],
            k_forward=kb,
            k_reverse=ku,
        )

        rxn_all = (
            rxn_open_p
            + rxn_iso
            + rxn_release_open
            + rxn_release_closed
            + [rxn_m1]
        )

        return rxn_all


class multi_tl(Mechanism):
    r"""Multi-ribosome translation with isomerization and occupancy.

    A 'translation' mechanism that explicitly models multiple ribosomes
    binding to a single mRNA simultaneously, accounting for ribosome spacing,
    isomerization between open and closed configurations, and competitive
    binding. This detailed mechanism captures translational queueing and
    ribosome traffic effects.

    The reaction scheme follows:
    $$
        & 'mRNA':'RBZ'_n + 'RBZ' <--> 'mRNA':'RBZ'_{n_c}
            --> 'mRNA':'RBZ'_{n+1} \\
        & 'mRNA':'RBZ'_n --> 'mRNA':'RBZ'_0 + n 'RBZ' + n 'Protein' \\
        & 'mRNA':'RBZ'_{n_c} --> 'mRNA':'RBZ'_{0_c} + n 'RBZ' + n 'Protein'
    $$
    where:

    - $n$ = {0, 1, ..., max_occ} is the number of ribosomes
    - max_occ is the physical maximum based on ribosome and mRNA dimensions
    - mRNA:RBZ$_n$ represents mRNA with $n$ open configuration ribosomes
    - mRNA:RBZ$_{n_c}$ represents mRNA with $n$ open and 1 closed ribosome

    Parameters
    ----------
    ribosome : Species
        Ribosome species.
    name : str, default='multi_tl'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='translation'
        Type classification of this mechanism.

    Attributes
    ----------
    ribosome : Species
        The ribosome species.
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('translation').

    See Also
    --------
    Translation_MM : Simple Michaelis-Menten translation.
    multi_tx : Multi-polymerase transcription mechanism.
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism provides a detailed, spatially-aware model of translation
    that captures:

    - Multiple ribosomes on a single mRNA (polysome formation)
    - Ribosome queueing and spacing constraints
    - Open/closed conformational states (isomerization)
    - Coordinated protein release from multiple ribosomes

    The model is appropriate for detailed mechanistic studies where:

    - Ribosome density and spacing matter
    - Translational queueing affects expression
    - Multiple concurrent translation events are important
    - Polysome dynamics are relevant

    The mechanism generates many species (2 * max_occ complexes) and
    reactions (O(max_occ)), so it should be used with caution for
    computational efficiency.

    CAUTION: This mechanism is still under development. Use with care, read
    all warnings, and consult the example notebook before use.

    Required parameters for this mechanism:

    - 'ktl' : Translation/release rate constant
    - 'kb' : Forward binding rate for ribosome to mRNA
    - 'ku' : Reverse unbinding rate for ribosome from mRNA
    - 'k_iso' : Isomerization rate from closed to open configuration
    - 'max_occ' : Maximum ribosome occupancy (integer)

    Examples
    --------
    Model multi-ribosome translation:

    >>> gene = bcp.DNAassembly(
    ...     name='gene_assembly',
    ...     promoter='pconst', transcript='mRNA',
    ...     rbs='RBS_medium', protein='GFP')
    >>> ribosome = bcp.Species('Ribosome')
    >>> tl_mechanism = bcp.multi_tl(ribosome=ribosome)
    >>> mixture = bcp.Mixture(
    ...     components=[gene],
    ...     mechanisms={
    ...         'transcription': bcp.SimpleTranscription(),
    ...         'translation': tl_mechanism},
    ...     parameters={
    ...         'ktx': 0.05, 'ktl': 0.1, 'kb': 1.0, 'ku': 0.1,
    ...         'k_iso': 0.5, 'max_occ': 5
    ...     }
    ... )
    >>> mixture.compile_crn()

    For more details, see examples/MultiTX_Demo.ipynb.

    """

    def __init__(
        self,
        ribosome: Species,
        name: str = 'multi_tl',
        mechanism_type: str = 'translation',
        **kwargs,
    ):
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("'ribosome' must be a Species.")

        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    # species update
    def update_species(
        self, transcript, protein, component, part_id, **kwargs
    ):
        """Generate species for multi-ribosome translation.

        Creates all species and complexes involved in multi-ribosome
        translation including mRNA with varying numbers of bound ribosomes in
        both open and closed configurations.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated.
        protein : Species
            The protein species produced by translation.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup. Required to retrieve 'max_occ'
            parameter.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing:

            - mRNA:Ribosome_n complexes in open configuration (n = 0 to
              max_occ-1)
            - mRNA:Ribosome_n complexes in closed configuration (n = 0 to
              max_occ-1)
            - Individual species [ribosome, transcript, protein]

        Notes
        -----
        For each occupancy level $n$ from 0 to max_occ-1, two complex species
        are created:

        - Open configuration: mRNA with n+1 ribosomes, all in open state
        - Closed configuration: mRNA with n+1 ribosomes (n open, 1 closed)

        The 'max_occ' parameter determines the maximum number of ribosomes
        that can occupy the mRNA simultaneously, typically based on physical
        spacing constraints and mRNA length.

        """
        max_occ = int(
            component.get_parameter(
                'max_occ',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        )
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=['open'],
                )
            )
            cp_closed.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=['closed'],
                )
            )

        cp_misc = [self.ribosome, transcript, protein]

        return cp_open + cp_closed + cp_misc

    def update_reactions(
        self, transcript, protein, component, part_id, **kwargs
    ):
        """Generate reactions for multi-ribosome translation.

        Creates reactions modeling multiple ribosomes binding to mRNA,
        isomerization between configurations, and protein release with
        coordination among multiple ribosomes.

        Parameters
        ----------
        transcript : Species
            The mRNA transcript species being translated.
        protein : Species
            The protein species produced by translation.
        component : Component
            Component containing parameter values. Required for parameter
            lookup.
        part_id : str
            Identifier for parameter lookup. Required to retrieve parameters.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List of reactions including:

            - Ribosome binding reactions (mRNA:RBZ$_n$ + RBZ <-->
              mRNA:RBZ$_{n_c}$)
            - Isomerization reactions (mRNA:RBZ$_{n_c}$ --> mRNA:RBZ$_n$)
            - Release reactions from open states (mRNA:RBZ$_n$ --> mRNA +
              $n$ RBZ + $n$ Protein)
            - Release reactions from closed states (mRNA:RBZ$_{n_c}$ -->
              mRNA:RBZ$_{0_c}$ + $n$ RBZ + $n$ Protein)
            - Base binding reaction (mRNA + RBZ <--> mRNA:RBZ$_{0_c}$)

        Notes
        -----
        The reaction scheme captures:

        1. mRNA:RBZ$_n$ + RBZ <--> mRNA:RBZ_${(n+1)_c}$ (rates: 'kb', 'ku')
        2. mRNA:RBZ$_{n_c}$ --> mRNA:RBZ$_n$ (rate: 'k_iso')
        3. mRNA:RBZ$_n$ --> mRNA + $n$ RBZ + $n$ Protein (rate: 'ktl')
        4. mRNA:RBZ$_{n_c}$ --> mRNA:RBZ$_{0_c}$ + $n$ RBZ + $n$ Protein
           (rate: 'ktl')

        Where $n$ ranges from 0 to max_occ-1. The 'max_occ' parameter
        represents the physical maximum occupancy of ribosomes on the mRNA,
        determined by ribosome spacing and mRNA length.

        """
        # parameter loading
        kb = component.get_parameter('kb', part_id=part_id, mechanism=self)
        ku = component.get_parameter('ku', part_id=part_id, mechanism=self)
        k_iso = component.get_parameter(
            'k_iso', part_id=part_id, mechanism=self
        )
        ktl = component.get_parameter('ktl', part_id=part_id, mechanism=self)
        max_occ = int(
            component.get_parameter(
                'max_occ',
                part_id=part_id,
                mechanism=self,
                return_numerical=True,
            )
        )

        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=['open'],
                )
            )
            cp_closed.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=['closed'],
                )
            )

        # Reactions
        # ribosome + complex(n) <--> complex(n+1)_closed
        rxn_open_p = [
            Reaction.from_massaction(
                inputs=[self.ribosome, cp_open[n]],
                outputs=[cp_closed[n + 1]],
                k_forward=kb,
                k_reverse=ku,
            )
            for n in range(0, max_occ - 1)
        ]

        # isomerization
        # complex(n)_closed --> complex(n)
        rxn_iso = [
            Reaction.from_massaction(
                inputs=[cp_closed[n]], outputs=[cp_open[n]], k_forward=k_iso
            )
            for n in range(0, max_occ)
        ]

        # release/translation from open and closed states
        rxn_release_open = []
        rxn_release_closed = []
        for n in range(0, max_occ):
            rxn_temp1 = Reaction.from_massaction(
                inputs=[cp_open[n]],
                outputs=[self.ribosome for i in range(n + 1)]
                + [protein for i in range(n + 1)]
                + [transcript],
                k_forward=ktl,
            )
            rxn_release_open.append(rxn_temp1)

        for n in range(1, max_occ):
            rxn_temp2 = Reaction.from_massaction(
                inputs=[cp_closed[n]],
                outputs=[self.ribosome for i in range(n)]
                + [protein for i in range(n)]
                + [cp_closed[0]],
                k_forward=ktl,
            )
            rxn_release_closed.append(rxn_temp2)

        # missing reactions (0 --> 0_closed and v.v. 0_closed --> 0)
        rxn_m1 = Reaction.from_massaction(
            inputs=[transcript, self.ribosome],
            outputs=[cp_closed[0]],
            k_forward=kb,
            k_reverse=ku,
        )
        # rxn_m2 = Reaction.from_massaction(
        #     inputs=[cp_closed[0]], outputs=[transcript, self.ribosome],
        #     k_forward=kur)

        rxn_all = (
            [rxn_m1]
            + rxn_iso
            + rxn_open_p
            + rxn_release_open
            + rxn_release_closed
        )

        return rxn_all
