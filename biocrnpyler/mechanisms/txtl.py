from typing import List
from ..core.mechanism import Mechanism
from ..core.propensities import ProportionalHillNegative, ProportionalHillPositive
from ..core.reaction import Reaction
from ..core.species import Complex, Species
from ..utils import parameter_to_value
from .enzyme import MichaelisMentenCopy


class OneStepGeneExpression(Mechanism):
    """A mechanism to model gene expression without transcription or translation.

    G --> G + P
    """

    def __init__(self, name="gene_expression", mechanism_type="transcription"):
        """Initializes a OneStepGeneExpression instance.

        :param name: name of the Mechanism, default: gene_expression
        :param mechanism_type: type of the Mechanism, default: transcription

        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, protein, transcript=None, **keywords):
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
        **keywords,
    ):

        if kexpress is None and component is not None:
            kexpress = component.get_parameter(
                "kexpress", part_id=part_id, mechanism=self
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
    """A Mechanism to model simple catalytic transcription.

    G --> G + T
    """

    def __init__(self, name="simple_transcription", mechanism_type="transcription"):
        """Initializes a SimpleTranscription instance.

        :param name: name of the Mechanism, default: simple_transcription
        :param mechanism_type: type of the Mechanism, default: transcription

        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, dna, transcript=None, protein=None, **keywords):

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
        **keywords,
    ):

        if ktx == None and component != None:
            ktx = component.get_parameter("ktx", part_id=part_id, mechanism=self)
        elif component == None and ktx == None:
            raise ValueError("Must pass in component or a value for ktx")

        # First case only true in Mixtures without transcription (eg Expression Mixtures)
        if transcript is None and protein is not None:
            ktl = component.get_parameter("ktl", part_id=part_id, mechanism=self)
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
    """A mechanism to model simple catalytic translation.

    T --> T + P
    """

    def __init__(self, name="simple_translation", mechanism_type="translation"):
        """Initializes a SimpleTranslation instance.

        :param name: name of the Mechanism, default: simple_translation
        :param mechanism_type: type of the Mechanism, default: translation

        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, transcript, protein=None, **keywords):
        if protein is None:
            protein = Species(transcript.name, material_type="protein")
        outlst = [transcript]
        if type(protein) == list:
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
        **keywords,
    ):

        if ktl is None and component is not None:
            ktl = component.get_parameter("ktl", part_id=part_id, mechanism=self)
        elif component is None and ktl is None:
            raise ValueError("Must pass in component or a value for ktl")

        # First case only true in Mixtures without transcription (eg Expression Mixtures)
        if transcript is None and protein is not None:
            rxns = []
        else:
            rxns = [
                Reaction.from_massaction(
                    inputs=[transcript], outputs=[transcript, protein], k_forward=ktl
                )
            ]

        return rxns


class PositiveHillTranscription(Mechanism):
    """A mechanism to model transcription as a proprotional positive hill function:
    G --> G + P
    rate = k*G*(R^n)/(K+R^n)
    where R is a regulator (activator).
    Optionally includes a leak reaction
    G --> G + P @ rate kleak.
    """

    def __init__(
        self, name="positivehill_transcription", mechanism_type="transcription"
    ):
        """Initializes a PositiveHillTranscription instance.

        :param name: name of the Mechanism, default: positivehill_transcription
        :param mechanism_type: type of the Mechanism, default: transcription

        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(
        self, dna, regulator, transcript=None, leak=False, protein=None, **keywords
    ):

        species = [dna, regulator]
        if transcript is not None:
            species += [transcript]
        if protein is not None:
            species += [protein]

        return species  # it is best to return all species that will be involved in the reactions

    def update_reactions(
        self,
        dna,
        regulator,
        component,
        part_id,
        transcript=None,
        leak=False,
        protein=None,
        **keywords,
    ):
        """This always requires the inputs component and part_id to find the relevant parameters

        :param dna:
        :param regulator:
        :param component:
        :param part_id:
        :param transcript:
        :param leak:
        :param protein:
        :param keywords:
        :return:
        """

        ktx = component.get_parameter("k", part_id=part_id, mechanism=self)
        n = component.get_parameter("n", part_id=part_id, mechanism=self)
        K = component.get_parameter("K", part_id=part_id, mechanism=self)
        kleak = component.get_parameter("kleak", part_id=part_id, mechanism=self)

        prophill = ProportionalHillPositive(k=ktx, K=K, s1=regulator, n=n, d=dna)

        reactions = []

        # First case only true in Mixtures without transcription (eg Expression Mixtures)
        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        reactions.append(
            Reaction(inputs=[dna], outputs=[dna, tx_output], propensity_type=prophill)
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
    """A mechanism to model transcription as a proprotional negative hill function:
    G --> G + P
    rate = k*G*(1)/(K+R^n)
    where R is a regulator (repressor).
    Optionally includes a leak reaction
    G --> G + P @ rate kleak.
    """

    def __init__(
        self, name="negativehill_transcription", mechanism_type="transcription"
    ):
        """Initializes a NegativeHillTranscription instance.

        :param name: name of the Mechanism, default: negativehill_transcription
        :param mechanism_type: type of the Mechanism, default: transcription

        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(
        self, dna, regulator, transcript=None, leak=False, protein=None, **keywords
    ):

        species = [dna, regulator]
        if transcript is not None:
            species += [transcript]
        if protein is not None:
            species += [protein]

        return species  # it is best to return all species that will be involved in the reactions

    def update_reactions(
        self,
        dna,
        regulator,
        component,
        part_id,
        transcript=None,
        leak=False,
        protein=None,
        **keywords,
    ):
        """This always requires the inputs component and part_id to find the relevant parameters

        :param dna:
        :param regulator:
        :param component:
        :param part_id:
        :param transcript:
        :param leak:
        :param protein:
        :param keywords:
        :return:
        """

        ktx = component.get_parameter("k", part_id=part_id, mechanism=self)
        n = component.get_parameter("n", part_id=part_id, mechanism=self)
        K = component.get_parameter("K", part_id=part_id, mechanism=self)
        kleak = component.get_parameter("kleak", part_id=part_id, mechanism=self)

        prop_hill = ProportionalHillNegative(k=ktx, K=K, n=n, s1=regulator, d=dna)

        reactions = []

        # First case only true in Mixtures without transcription (eg Expression Mixtures)
        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        reactions.append(
            Reaction(inputs=[dna], outputs=[dna, tx_output], propensity_type=prop_hill)
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
    """Michaelis Menten Transcription.

    G + RNAP <--> G:RNAP --> G+RNAP+mRNA
    """

    def __init__(self, rnap: Species, name="transcription_mm", **keywords):
        """Initializes a Transcription_MM instance.

        :param rnap: Species instance that is representing an RNA polymerase
        :param name: name of the Mechanism, default: transcription_mm
        """
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' parameter must be a Species.")

        MichaelisMentenCopy.__init__(
            self=self, name=name, mechanism_type="transcription"
        )

    def update_species(self, dna, transcript=None, protein=None, **keywords):
        species = [dna]

        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        species += MichaelisMentenCopy.update_species(
            self, Enzyme=self.rnap, Sub=dna, Prod=tx_output
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
        **keywords,
    ):

        # Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktx = component.get_parameter("ktx", part_id=part_id, mechanism=self)
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)

        rxns = []

        if transcript is None and protein is not None:
            tx_output = protein
        else:
            tx_output = transcript

        rxns += MichaelisMentenCopy.update_reactions(
            self,
            Enzyme=self.rnap,
            Sub=dna,
            Prod=tx_output,
            complex=complex,
            kb=kb,
            ku=ku,
            kcat=ktx,
        )

        return rxns


class Translation_MM(MichaelisMentenCopy):
    """Michaelis Menten Translation.

    mRNA + Rib <--> mRNA:Rib --> mRNA + Rib + Protein
    """

    def __init__(self, ribosome: Species, name="translation_mm", **keywords):
        """Initializes a Translation_MM instance.

        :param ribosome: Species instance that is representing a ribosome
        :param name: name of the Mechanism, default: translation_mm
        """
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("ribosome must be a Species!")
        MichaelisMentenCopy.__init__(self=self, name=name, mechanism_type="translation")

    def update_species(self, transcript, protein, **keywords):
        species = []

        # This can only occur in expression mixtures
        if transcript is None and protein is not None:
            species += Species.flatten_list([protein])
        else:
            species += MichaelisMentenCopy.update_species(
                self, Enzyme=self.ribosome, Sub=transcript, Prod=protein
            )

        return species

    def update_reactions(
        self, transcript, protein, component, part_id=None, complex=None, **keywords
    ):
        rxns = []

        # Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktl = component.get_parameter("ktl", part_id=part_id, mechanism=self)
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)

        # This can only occur in expression mixtures
        if transcript is None and protein is not None:
            pass
        else:
            rxns += MichaelisMentenCopy.update_reactions(
                self,
                Enzyme=self.ribosome,
                Sub=transcript,
                Prod=protein,
                complex=complex,
                kb=kb,
                ku=ku,
                kcat=ktl,
            )
        return rxns


class Energy_Transcription_MM(Mechanism):
    """Michaelis Menten Transcription that consumed energy.

    G + RNAP <--> G:RNAP
    Fuel + G:RNAP --> G + RNAP + T + Fuel (Transcription can only happen when there is fuel)
        at rate ktx/L (length dependent transcription rate)

    Fuel + G:RNAP --> G:RNAP +wastes (Fuel consumption treated faster)
        at rate ktx (This occurs L times faster than the above, resulting in the correct fuel use)
    """

    def __init__(
        self,
        rnap: Species,
        fuels: List[Species],
        wastes=List[Species],
        name="energy_transcription_mm",
        **keywords,
    ):
        """Initializes a Transcription_MM instance.
        :param fuels: List of Species consumed during transcription
        :param wastes: List of Species consumed during transcription
        :param rnap: Species instance that is representing an RNA polymerase
        :param name: name of the Mechanism, default: transcription_mm
        """
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' parameter must be a Species.")

        if all([isinstance(s, Species) for s in fuels]):
            self.fuels = fuels
        else:
            raise ValueError("recieved a non-Species object in the fuels list.")

        if all([isinstance(s, Species) for s in fuels]):
            self.wastes = wastes
        else:
            raise ValueError("wastes must be a list of Species!")

        Mechanism.__init__(self=self, name=name, mechanism_type="transcription")

    def update_species(self, dna, transcript=None, protein=None, **keywords):
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
        **keywords,
    ):

        # Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktx = component.get_parameter("ktx", part_id=part_id, mechanism=self)
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        L = component.get_parameter("length", part_id=part_id, mechanism=self)

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
            self.fuels + [bound_complex], [bound_complex] + self.wastes, k_forward=ktx
        )

        return [r1, r2, r3]


class Energy_Translation_MM(Mechanism):
    """Michaelis Menten Translation that consumes energy species.

    mRNA + Rib  <--> mRNA:Rib  (binding)
    fuels + mRNA:Rib --> mRNA + Rib + Protein + fuels    (translation)
    fuels + mRNA:Rib --> mRNA:Rib  +wastes (fuel consumption)
    """

    def __init__(
        self,
        ribosome: Species,
        fuels: List[Species],
        wastes=List[Species],
        name="energy_translation_mm",
        **keywords,
    ):
        """Initializes a Translation_MM instance.

        :param ribosome: Species instance that is representing a ribosome
        :param fuels: List of fuel Species that are consumed during translation
        :param wastes: List of Species consumed during translation
        :param name: name of the Mechanism, default: energy_translation_mm
        """
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

        Mechanism.__init__(self=self, name=name, mechanism_type="translation")

    def update_species(self, transcript, protein, **keywords):
        species = self.fuels + [self.ribosome, protein]
        bound_complex = Complex([transcript, self.ribosome])
        species += [bound_complex]

        return species

    def update_reactions(
        self, transcript, protein, component, part_id=None, complex=None, **keywords
    ):
        rxns = []

        # Get Parameters
        if part_id == None and component != None:
            part_id = component.name

        ktl = component.get_parameter("ktl", part_id=part_id, mechanism=self)
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        L = component.get_parameter("length", part_id=part_id, mechanism=self)

        bound_complex = Complex([transcript, self.ribosome])

        # RNAP DNA Binding
        r1 = Reaction.from_massaction(
            [transcript, self.ribosome], [bound_complex], k_forward=kb, k_reverse=ku
        )
        # Transcription
        r2 = Reaction.from_massaction(
            self.fuels + [bound_complex],
            self.fuels + [transcript, self.ribosome, protein],
            k_forward=parameter_to_value(ktl.value) / parameter_to_value(L),
        )
        # Fuel consumption
        r3 = Reaction.from_massaction(
            self.fuels + [bound_complex], [bound_complex] + self.wastes, k_forward=ktl
        )

        return [r1, r2, r3]


class multi_tx(Mechanism):
    """Multi-RNAp Transcription w/ Isomerization.

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
    """

    def __init__(
        self,
        pol: Species,
        name: str = "multi_tx",
        mechanism_type: str = "transcription",
        **keywords,
    ):
        """Initializes a multi_tx instance.

        :param pol: reference to a species instance that represents a polymerase
        :param name: name of the Mechanism, default: multi_tx
        :param mechanism_type: type of the mechanism, default: transcription
        :param keywords:
        """
        if isinstance(pol, Species):
            self.pol = pol
        else:
            raise ValueError("'pol' must be a Species")

        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    # species update
    def update_species(
        self, dna, transcript, component, part_id, protein=None, **keywords
    ):
        max_occ = int(
            component.get_parameter(
                "max_occ", part_id=part_id, mechanism=self, return_numerical=True
            )
        )
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex([dna] + [self.pol for i in range(n + 1)], attributes=["open"])
            )
            cp_closed.append(
                Complex([dna] + [self.pol for i in range(n + 1)], attributes=["closed"])
            )

        cp_misc = [self.pol, dna, transcript]

        return cp_open + cp_closed + cp_misc

    def update_reactions(
        self, dna, transcript, component, part_id, protein=None, **keywords
    ):
        """It sets up the following reactions.

        DNA:RNAp_n + RNAp <--> DNA:RNAp_n_c --> DNA:RNAp_n+1
        kf1 = k1, kr1 = k2, kf2 = k_iso
        DNA:RNAp_n --> DNA:RNAp_0 + n RNAp + n mRNA
        kf = ktx_solo
        DNA:RNAp_n_c --> DNA:RNAp_0_c + n RNAp + n mRNA
        kf = ktx_solo

        max_occ =  maximum occupancy of gene (physical limit)
        """

        # parameter loading
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        k_iso = component.get_parameter("k_iso", part_id=part_id, mechanism=self)
        ktx = component.get_parameter("ktx", part_id=part_id, mechanism=self)
        max_occ = int(
            component.get_parameter(
                "max_occ", part_id=part_id, mechanism=self, return_numerical=True
            )
        )

        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            # has n polymerases all open
            cp_open.append(
                Complex([dna] + [self.pol for i in range(n + 1)], attributes=["open"])
            )
            # has n-1 open polymerases and 1 closed polymerase
            cp_closed.append(
                Complex([dna] + [self.pol for i in range(n + 1)], attributes=["closed"])
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
        # rxn_open_pr = [Reaction.from_massaction(inputs=[cp_closed[n + 1]], outputs=[self.pol, cp_open[n], ], k_forward=k2) for n in range(0, max_occ - 1)]

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
            inputs=[dna, self.pol], outputs=[cp_closed[0]], k_forward=kb, k_reverse=ku
        )

        rxn_all = (
            rxn_open_p + rxn_iso + rxn_release_open + rxn_release_closed + [rxn_m1]
        )

        return rxn_all


class multi_tl(Mechanism):
    """Multi-RBZ Translation w/ Isomerization.

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
    """

    def __init__(
        self,
        ribosome: Species,
        name: str = "multi_tl",
        mechanism_type: str = "translation",
        **keywords,
    ):
        """Initializes a multi_tl instance.

        :param ribosome: a Species instance that represents a ribosome
        :param name: name of the Mechanism, default: multi_tl
        :param mechanism_type: type of the Mechanism, default: translation

        """
        if isinstance(ribosome, Species):
            self.ribosome = ribosome
        else:
            raise ValueError("'ribosome' must be a Species.")

        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    # species update
    def update_species(self, transcript, protein, component, part_id, **keywords):
        max_occ = int(
            component.get_parameter(
                "max_occ", part_id=part_id, mechanism=self, return_numerical=True
            )
        )
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=["open"],
                )
            )
            cp_closed.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=["closed"],
                )
            )

        cp_misc = [self.ribosome, transcript, protein]

        return cp_open + cp_closed + cp_misc

    def update_reactions(self, transcript, protein, component, part_id, **keywords):
        """It sets up the following reactions.

        mRNA:RBZ_n + RBZ <--> mRNA:RBZ_n_c --> mRNA:RBZ_n+1
        kf1 = kbr, kr1 = kur, kf2 = k_iso_r
        mRNA:RBZ_n --> mRNA:RBZ_0 + n RBZ + n Protein
        kf = ktl_solo
        mRNA:RBZ_n_c --> mRNA:RBZ_0_c + n RBZ + n Protein
        kf = ktl_solo
        """

        # parameter loading
        kb = component.get_parameter("kb", part_id=part_id, mechanism=self)
        ku = component.get_parameter("ku", part_id=part_id, mechanism=self)
        k_iso = component.get_parameter("k_iso", part_id=part_id, mechanism=self)
        ktl = component.get_parameter("ktl", part_id=part_id, mechanism=self)
        max_occ = int(
            component.get_parameter(
                "max_occ", part_id=part_id, mechanism=self, return_numerical=True
            )
        )

        # complex species instantiation
        cp_open = []
        cp_closed = []
        for n in range(0, max_occ):
            cp_open.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=["open"],
                )
            )
            cp_closed.append(
                Complex(
                    [transcript] + [self.ribosome for i in range(n + 1)],
                    attributes=["closed"],
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
        # rxn_m2 = Reaction.from_massaction(inputs=[cp_closed[0]], outputs=[transcript, self.ribosome], k_forward=kur)

        rxn_all = (
            [rxn_m1] + rxn_iso + rxn_open_p + rxn_release_open + rxn_release_closed
        )

        return rxn_all
