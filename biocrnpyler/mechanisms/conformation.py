# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from ..core.mechanism import Mechanism
from ..core.reaction import Reaction


class One_Step_Reversible_Conformation_Change(Mechanism):
    r"""Reversible conformational change mechanism.

    A mechanism that models the reversible conformational change of a species
    from one state to another. This can represent protein folding/unfolding,
    DNA structural changes, allosteric transitions, or any other molecular
    conformational switch. Additional species (cofactors, ions, etc.) can be
    required for the conformational change.

    The reaction follows:

    $$s_0 + \text{additional species} \leftrightarrow s_f + \text{additional species}$$

    where s0 is the initial conformation and sf is the final conformation.

    Parameters
    ----------
    name : str, default='one_step_conformation_change'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='conformation_change'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('conformation_change').

    See Also
    --------
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    This mechanism is used to model conformational changes including:

    - Protein folding: Native <--> denatured states
    - Allosteric transitions: Inactive <--> active enzyme forms
    - DNA structural changes: B-DNA <--> Z-DNA transitions
    - Receptor activation: Closed <--> open channel states
    - Molecular switches: Any two-state molecular system

    The mechanism requires two rate constants:

    - 'kf': Forward rate constant (s0 -> sf)
    - 'kr': Reverse rate constant (sf -> s0)

    Additional species can participate in the conformational change without
    being modified (e.g., ATP required for a conformational switch but not
    consumed).

    """

    def __init__(
        self,
        name='one_step_conformation_change',
        mechanism_type='conformation_change',
    ):
        """Initialize a One_Step_Reversible_Conformation_Change mechanism.

        See class docstring for parameter descriptions.

        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(
        self,
        s0,
        sf,
        additional_species=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate species for conformational change.

        Creates the list of species involved in the conformational change,
        including initial and final conformations plus any additional species
        required for the transition.

        Parameters
        ----------
        s0 : Species
            The initial conformation state of the molecule.
        sf : Species
            The final conformation state of the molecule.
        additional_species : list of Species, optional
            Additional species required for the conformational change (e.g.,
            cofactors, ions) but not consumed. Default is empty list.
        component : Component, optional
            Component containing this mechanism (for context, not used here).
        part_id : str, optional
            Identifier for parameter lookup (not used in species generation).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [s0, sf] plus any additional species.

        Notes
        -----
        The additional species participate in the conformational change but
        are not modified. They appear on both sides of the reaction. This is
        useful for modeling:

        - Cofactor-dependent conformational changes
        - Ion-induced structural transitions
        - Ligand-stabilized conformations

        """
        if additional_species is None:
            additional_species = []

        return [s0, sf] + additional_species

    def update_reactions(
        self,
        s0,
        sf,
        additional_species=None,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for conformational change.

        Creates a reversible mass-action reaction for the conformational
        transition between two states of a molecule.

        Parameters
        ----------
        s0 : Species
            The initial conformation state of the molecule.
        sf : Species
            The final conformation state of the molecule.
        additional_species : list of Species, optional
            Additional species required for the conformational change (e.g.,
            cofactors, ions) but not consumed. These appear on both sides of
            the reaction. Default is empty list.
        component : Component
            Component containing parameter values. Required for retrieving
            'kf' and 'kr' rate constants.
        part_id : str, optional
            Identifier for parameter lookup. If None, defaults to
            'repr(s0)-repr(sf)'.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing a single reversible mass-action reaction for the
            conformational change.

        Raises
        ------
        ValueError
            If component is None (implicitly, when parameters cannot be
            retrieved).

        Notes
        -----
        The reaction equation depends on additional_species:

        - Without additional species: s0 <--> sf
        - With additional species: s0 + additionals <--> sf + additionals

        The mechanism requires two parameters:

        - 'kf': Forward rate constant (s0 -> sf transition)
        - 'kr': Reverse rate constant (sf -> s0 transition)

        Additional species are not consumed or produced; they facilitate the
        conformational change. This models situations where cofactors or ions
        must be present for the transition but are not modified.

        """
        if part_id is None:
            part_id = repr(s0) + '-' + repr(sf)

        if additional_species is None:
            additional_species = []

        kf = component.get_parameter('kf', part_id=part_id, mechanism=self)
        kr = component.get_parameter('kr', part_id=part_id, mechanism=self)

        return [
            Reaction.from_massaction(
                inputs=[s0] + additional_species,
                outputs=[sf],
                k_forward=kf,
                k_reverse=kr,
            )
        ]
