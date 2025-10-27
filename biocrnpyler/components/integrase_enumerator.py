#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.
import copy
import itertools as it
from typing import List

from ..core.polymer import NamedPolymer, OrderedMonomer
from ..core.species import ComplexSpecies, Species
from ..utils import combine_dictionaries
from .component_enumerator import GlobalComponentEnumerator
from .dna.construct import Construct, DNAconstruct
from .dna.misc import IntegraseSite


class PolymerTransformation:
    """Template for transforming polymer sequences through recombination.

    A `PolymerTransformation` defines a template for creating new polymers
    from existing ones through recombination reactions. The template
    specifies a parts list containing placeholders (monomers from input
    polymers) and new parts (with no parent). This enables complex DNA
    rearrangements like integration, deletion, and inversion.

    Parameters
    ----------
    partslist : list
        List of parts defining the output polymer. Can contain:

        - OrderedMonomers from existing polymers (placeholders)
        - Parts with parent=None (inserted into new polymer)
        - Tuples of (part, direction)

    circular : bool, default=False
        Whether the output polymer should be circular.
    parentsdict : dict, optional
        Dictionary mapping parent polymers to input names ('input1',
        'input2', etc.). If None, automatically generated.
    material_type : str, default='dna'
        Material type for the created polymer.

    Attributes
    ----------
    number_of_inputs : int
        Number of distinct input polymers required.
    parentsdict : dict
        Mapping from parent polymers to generic input names.
    partslist : list
        Template parts list with dummy placeholders.
    circular : bool
        Whether output is circular.
    material_type : str
        Material type of output polymer.

    See Also
    --------
    IntegraseRule : Defines integrase recombination rules.
    Integrase_Enumerator : Enumerates integrase products.
    OrderedMonomer : Monomer in ordered polymer.

    Notes
    -----
    **Template Mechanism:**

    The transformation works by:

    1. Analyzing partslist to identify parent polymers
    2. Creating generic 'input#' placeholders for each parent
    3. Storing template with dummy placeholders
    4. When applied via `create_polymer`, replacing placeholders with
       actual parts from input polymers

    **Placeholder System:**

    - Parts from polymers become placeholders referencing position in
      'input#'
    - Parts with parent=None are copied directly into output
    - Complexes bound to parts are transferred to new positions

    Examples
    --------
    Create a simple transformation template:

    >>> # Define template: take element 0 from input1, element 2 from
    >>> # input2 (reversed), and insert a promoter
    >>> template = PolymerTransformation(
    ...     partslist=[
    ...         polymer1[0],  # Placeholder for position 0
    ...         [polymer2[2], 'reverse'],  # Position 2, reversed
    ...         promoter  # New part (parent=None)
    ...     ],
    ...     circular=False
    ... )
    >>> # Apply template to create new polymer
    >>> new_polymer = template.create_polymer([polymer1, polymer2])

    Integration reaction template:

    >>> # Combine two plasmids at cut sites
    >>> template = PolymerTransformation(
    ...     partslist=(
    ...         plasmid1[:cut1] +
    ...         [[prod_site1, 'forward']] +
    ...         plasmid2[cut2+1:] +
    ...         [[prod_site2, 'forward']] +
    ...         plasmid1[cut1+1:]
    ...     ),
    ...     circular=True
    ... )

    """

    def __init__(
        self, partslist, circular=False, parentsdict=None, material_type='dna'
    ):
        if parentsdict is None:
            # the input to this function is a list of monomers that belong to
            # various parents.  Each different parent that is represented in
            # these monomers is converted into "input#" arbitrarily to keep it
            # consistant, a link of that parent to the appropriate "input#" is
            # kept in this dictionary.  You can also pass a pre-populated
            # dictionary into this function in order to control which parent
            # gets which "input#". This is essential if we want to properly
            # "reverse" a transformation (input1 becomes input2 and input2
            # becomes input1)
            parentsdict = {}
        inputcount = 1
        actual_partslist = []
        partdir = (
            -1
        )  # -1 = not set. Valid values are "forward" "reverse" and None
        # go through the parts
        inputname = None
        new_parentsdict = {}
        for part in partslist:
            if isinstance(part, (list, tuple)):
                # if the part is a list that means it looks like
                # [OrderedMonomer,"direction"]
                part_ref = part[0]
                partdir = part[1]
            else:
                part_ref = part
                partdir = part_ref.direction
            # if the parent is populated, it means this part should be a
            # placeholder
            if part_ref.parent is not None:
                # if we haven't tracked it already
                if (part_ref.parent not in new_parentsdict) and (
                    part_ref.parent not in parentsdict
                ):
                    # create a 'blank' polymer based on the input, and give it
                    # a generic name
                    inputname = 'input' + str(inputcount)
                    while inputname in list(new_parentsdict.values()):
                        inputcount += 1
                        inputname = 'input' + str(inputcount)
                    dummyPolymer = self.dummify(part_ref.parent, inputname)
                    # set the value in the dictionary
                    new_parentsdict[part_ref.parent] = dummyPolymer.name
                    # increment for the next time this happens
                    inputcount += 1

                elif (part_ref.parent not in new_parentsdict) and (
                    part_ref.parent in parentsdict
                ):
                    inputname = parentsdict[part_ref.parent]
                    dummyPolymer = self.dummify(part_ref.parent, inputname)
                    new_parentsdict[part_ref.parent] = dummyPolymer.name
                elif (
                    inputname is None
                    or inputname != new_parentsdict[part_ref.parent]
                ):
                    inputname = new_parentsdict[part_ref.parent]
                    dummyPolymer = self.dummify(part_ref.parent, inputname)
                # this variable is the actual parts list that will be
                # stored. It has mostly dummy parts
                if part_ref.parent[part_ref.position] != part_ref:
                    # this will be true only for situations where a component
                    # needs to be transformed into another component but still
                    # needs to keep a reference to where it was in the parent
                    copied_part = part_ref.get_orphan()
                    copied_part.parent = dummyPolymer
                    actual_partslist += [[copied_part, partdir]]
                else:
                    actual_partslist += [
                        [dummyPolymer[part_ref.position], partdir]
                    ]
            else:
                # if the part has no parent, copy it and put it in
                actual_partslist += [[copy.deepcopy(part_ref), partdir]]
        self.number_of_inputs = len(new_parentsdict)
        self.parentsdict = new_parentsdict
        self.partslist = actual_partslist
        self.circular = circular
        self.material_type = material_type

    def renumber_output(self, output_renumbering_function):
        """Change the ordering of the output parts list.

        Applies a renumbering function to rearrange the parts in the
        template, useful for handling circular permutations or reversals.

        Parameters
        ----------
        output_renumbering_function : callable
            Function that takes an index (int) and returns tuple of
            (new_index, direction), where direction is 'f' (forward) or
            'r' (reverse).

        Notes
        -----
        Modifies self.partslist in place. Used to adjust transformations
        when matching against existing constructs.

        """
        new_partslist = []
        for i in range(len(self.partslist)):
            new_i, direc = output_renumbering_function(i)
            part = self.partslist[new_i]
            if direc == 'r':
                new_partslist += [
                    [part[0], ['forward', 'reverse'][part[1] == 'forward']]
                ]
            else:
                new_partslist += [part]
        self.partslist = new_partslist

    def get_renumbered(self, output_renumbering_function):
        """Return copy of transformation with output indexes renumbered.

        Parameters
        ----------
        output_renumbering_function : callable
            Function mapping old indexes to (new_index, direction) tuples.

        Returns
        -------
        PolymerTransformation
            Copy of this transformation with renumbered output.

        """
        rxn_copied = copy.copy(self)
        rxn_copied.renumber_output(output_renumbering_function)
        return rxn_copied

    def reversed(self):
        """Return circularly permuted version with rotated inputs.

        Creates a new transformation where inputs are shuffled: input1
        becomes input2, input2 becomes input3, ..., last input becomes
        input1. Used for handling symmetric integrase reactions.

        Returns
        -------
        PolymerTransformation
            New transformation with rotated input assignments.

        Notes
        -----
        For single-input transformations, returns self unchanged. Essential
        for bidirectional integrase reactions where site1/site2 roles can
        be swapped.

        """
        new_parentsdict = {}

        new_name_dict = {}
        for inputnum in range(self.number_of_inputs):
            if inputnum == self.number_of_inputs - 1:
                new_name_dict['input' + str(inputnum + 1)] = 'input1'
            else:
                new_name_dict['input' + str(inputnum + 1)] = 'input' + str(
                    inputnum + 2
                )

        if len(self.parentsdict) == 1:
            # in this case we have only "input1"
            return self
        else:
            for part in self.partslist:
                if (
                    part[0].parent is not None
                    and part[0].parent not in new_parentsdict
                ):
                    # this next line outputs input2 when given input1

                    newname = new_name_dict[part[0].parent.name]
                    new_parentsdict[part[0].parent] = newname
            return PolymerTransformation(
                self.partslist, self.circular, parentsdict=new_parentsdict
            )

    def create_polymer(self, polymer_list, **kwargs):
        """Create a new polymer from the template using input polymers.

        Applies the transformation template to concrete input polymers,
        replacing placeholders with actual parts and inserting new parts.

        Parameters
        ----------
        polymer_list : list of Polymer
            List of input polymers. Must have at least number_of_inputs
            polymers. Order matters: first polymer is 'input1', second is
            'input2', etc.
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        Polymer
            New polymer created by applying the transformation. Type
            matches the input polymers' class.

        Notes
        -----
        **Transformation Process:**

        1. Map polymer_list to 'input#' names
        2. For each template part:

           - If placeholder: grab part from appropriate input polymer at
             specified position
           - If new part: insert the part or its dna_species
           - Handle complex species bound to parts

        3. Create output polymer with specified circularity

        **Complex Handling:**

        When transforming parts with bound complexes, the method attempts
        to preserve bindings by replacing core parts within complexes.

        """
        polymer_dict = {
            'input' + str(a + 1): b for a, b in enumerate(polymer_list)
        }
        assert len(polymer_list) >= self.number_of_inputs
        outlst = []
        for part_list in self.partslist:
            part = part_list[0]
            partdir = part_list[1]
            outpart = None

            if (part.parent is not None) and not hasattr(part, 'name'):
                outpart = polymer_dict[part.parent.name][
                    part.position
                ]  # grab the part from the proper input
            elif (part.parent is not None) and hasattr(part, 'name'):
                # in this case we are transforming a part into a different
                # part, taking the complexes from the previous position

                if isinstance(polymer_dict['input1'], Construct):
                    outpart = part.get_orphan()
                else:
                    template_part = polymer_dict[part.parent.name][
                        part.position
                    ]
                    if isinstance(template_part, ComplexSpecies):
                        # in this case we have to replace the thing inside the
                        # complex with the new part, but leave everything else
                        # the same.
                        old_species = template_part.get_species(
                            recursive=True
                        )
                        core_parts = []
                        for spec in old_species:
                            # if you have a material type of part then you are
                            # part of a polymer and thus will be replaced
                            if (
                                spec.material_type == 'part'
                                and not spec == template_part
                            ):
                                core_parts += [spec]

                        # Now we should have found only one core part. If
                        # there are multiple core parts then we might have
                        # to get crazy.
                        new_part = part.dna_species
                        if len(core_parts) == 1:
                            new_part = template_part.replace_species(
                                core_parts[0], part.dna_species
                            )
                        else:
                            raise KeyError(
                                f"{template_part} contained more than one "
                                "species with material_type='part', they "
                                f"were {core_parts}"
                            )
                        outpart = new_part
                    else:
                        outpart = part.dna_species
            else:
                # Parts that don't come from input1 or input2.  They need
                # to be either DNApart or species objects, depending on
                # what kind of object we are making in the end.
                if isinstance(polymer_dict['input1'], Construct):
                    outpart = part
                else:
                    outpart = part.dna_species
            # assuming the stored parts have a valid direction
            if hasattr(outpart, 'linked_sites'):
                outpart = copy.copy(outpart)
                # make sure that any integrase sites we copy this way have no
                # linked sites, as those would not be links created by the
                # integrate() function
                outpart.linked_sites = {}
            outlst += [[outpart, partdir]]
        if hasattr(outlst[0][0], 'material_type') and any(
            ['complex' in a[0].material_type for a in outlst]
        ):
            outpolymer = polymer_dict['input1'].__class__(
                outlst,
                circular=self.circular,
                material_type='ordered_polymer',
            )
        else:
            outpolymer = polymer_dict['input1'].__class__(
                outlst,
                circular=self.circular,
                material_type=self.material_type,
            )
        return outpolymer

    @classmethod
    def dummify(cls, in_polymer, name):
        """Create a simplified placeholder polymer.

        Generates a generic polymer with the same structure (length,
        directions, circularity) as the input but without specific parts.
        Used for creating template placeholders.

        Parameters
        ----------
        in_polymer : Polymer
            The polymer to simplify.
        name : str
            Name for the dummy polymer (e.g., 'input1').

        Returns
        -------
        NamedPolymer
            Simplified polymer with generic OrderedMonomers at each
            position.

        Notes
        -----
        The dummy polymer preserves only structural information (number of
        monomers, their directions, and circularity) while removing
        specific part identities.

        """
        # this is used specifically with polymerTransformation. Dummified
        # version of polymers are stored in polymerTransformation as generic
        # "slots" for monomers to be properly placed into
        out_list = []
        for element in in_polymer:
            out_list += [OrderedMonomer(direction=element.direction)]
        circular = False
        if hasattr(in_polymer, 'circular'):
            circular = in_polymer.circular
        return NamedPolymer(out_list, name, circular=circular)

    def __repr__(self):
        """Return string representation of the transformation.

        Returns
        -------
        str
            Human-readable string showing the transformation template with
            part names, positions, and directions.

        """
        part_texts = []
        for plist in self.partslist:
            part = plist[0]
            part_dir = plist[1]
            if part.parent is not None:
                if part.parent[part.position] != part:
                    part_texts += [[part.name, part.position, part_dir]]
                else:
                    part_texts += [
                        [part.parent.name, part.position, part_dir]
                    ]
            else:
                part_texts += [[part.name, part_dir]]
        out_txt = 'Polymer transformation = '
        for part_text in part_texts:
            out_txt += '(' + str(part_text[0])
            if len(part_text) == 3:
                out_txt += '-' + str(part_text[1])
            if part_text[-1] != 'forward':
                out_txt += '-r'
            out_txt += ')'
        if self.circular:
            out_txt += '(circular)'
        return out_txt


class IntegraseRule:
    """Rules defining integrase recombination reactions and products.

    An `IntegraseRule` specifies how an integrase enzyme acts on
    attachment sites to generate recombined DNA products. Defines which
    site pairs can react, what products they form, and which reaction
    types (deletion, integration, inversion) are allowed.

    Parameters
    ----------
    name : str, optional
        Name of the integrase (default='int1').
    reactions : dict, optional
        Dictionary mapping (site1_type, site2_type) tuples to product
        site type. Default: {('attB', 'attP'): 'attL',
        ('attP', 'attB'): 'attR'}
    allow_deletion : bool, default=True
        Whether to allow deletion reactions (intramolecular, same
        direction).
    allow_integration : bool, default=True
        Whether to allow integration reactions (intermolecular).
    allow_inversion : bool, default=True
        Whether to allow inversion reactions (intramolecular, opposite
        directions).

    Attributes
    ----------
    name : str
        Integrase name.
    integrase_species : Species
        The integrase protein species.
    reactions : dict
        Reaction rules mapping site pairs to products.
    attsites : list
        All attachment site types involved in reactions.
    allow_deletion : bool
        Whether deletions are allowed.
    allow_integration : bool
        Whether integrations are allowed.
    allow_inversion : bool
        Whether inversions are allowed.
    integrations_to_do : list
        List of integrations to perform during compilation.

    See Also
    --------
    IntegraseSite : DNA part representing attachment sites.
    Integrase_Enumerator : Enumerator using integrase rules.
    PolymerTransformation : Template for DNA transformations.

    Notes
    -----
    **Integrase Mechanism Types:**

    1. **Serine Integrases:**

       - Recombine attB + attP --> attL + attR
       - Require matching dinucleotides
       - With directionality factors: attL + attR --> attB + attP

    2. **Tyrosine Recombinases (Cre, Flp):**

       - Homotypic sites: loxP + loxP --> loxP + loxP
       - Can be palindromic (bidirectional)

    3. **Invertases:**

       - Only perform inversion reactions
       - Set allow_deletion=False, allow_integration=False

    4. **Resolvases:**

       - Only perform deletion reactions
       - Set allow_inversion=False, allow_integration=False

    **Reaction Types:**

    - **Inversion:** Two sites on same DNA, opposite directions -->
      region between sites flips
    - **Deletion:** Two sites on same DNA, same direction --> region
      between sites excised (forms circular product)
    - **Integration:** Sites on different DNAs --> DNAs join
    - **Recombination:** Two linear DNAs --> two recombinant linear DNAs

    Examples
    --------
    Define a standard serine integrase:

    >>> int_rule = bcp.IntegraseRule(
    ...     name='PhiC31',
    ...     reactions={
    ...         ('attB', 'attP'): 'attL',
    ...         ('attP', 'attB'): 'attR'
    ...     }
    ... )

    Define a Cre recombinase (homotypic):

    >>> cre_rule = bcp.IntegraseRule(
    ...     name='Cre',
    ...     reactions={
    ...         ('loxP', 'loxP'): 'loxP',
    ...         ('loxP', 'loxP'): 'loxP'  # Symmetric
    ...     }
    ... )

    Define an invertase (inversion only):

    >>> inv_rule = bcp.IntegraseRule(
    ...     name='Hin',
    ...     reactions={('hixL', 'hixR'): 'hixL', ('hixR', 'hixL'): 'hixR'},
    ...     allow_deletion=False,
    ...     allow_integration=False
    ... )

    """

    def __init__(
        self,
        name=None,
        reactions=None,
        allow_deletion=True,
        allow_integration=True,
        allow_inversion=True,
    ):
        if reactions is None:
            reactions = {('attB', 'attP'): 'attL', ('attP', 'attB'): 'attR'}
        if name is None:
            self.name = 'int1'
        else:
            self.name = name
        self.allow_deletion = allow_deletion
        self.allow_integration = allow_integration
        self.allow_inversion = allow_inversion
        self.integrase_species = Species(self.name, material_type='protein')
        self.reactions = reactions
        self.attsites = []
        for reaction in reactions:
            self.attsites += list(reaction) + [reactions[reaction]]
        self.attsites = list(set(self.attsites))
        # these are the reactions that will be performed at compile time
        self.integrations_to_do = []

    def binds_to(self):
        """Get all attachment site types this integrase binds to.

        Returns
        -------
        list of str
            List of all site types involved in integrase reactions.

        """
        return self.attsites

    def reaction_allowed(self, site1, site2):
        """Check if two sites can undergo integrase recombination.

        Parameters
        ----------
        site1 : IntegraseSite
            First attachment site.
        site2 : IntegraseSite
            Second attachment site.

        Returns
        -------
        bool
            True if sites can react according to reaction rules.

        Raises
        ------
        AssertionError
            If sites have different integrases or don't match this
            integrase.

        """
        assert isinstance(site1, IntegraseSite)
        assert isinstance(site2, IntegraseSite)
        assert site1.integrase == site2.integrase
        assert site1.integrase == self.integrase_species
        if tuple([site1.site_type, site2.site_type]) in self.reactions:
            return True
        return False

    def reactive_sites(self):
        """Get attachment site types that participate in reactions.

        Returns
        -------
        list of str
            List of site types that can be reactants (not just products).

        """
        attsites = []
        for reaction in self.reactions:
            attsites += list(reaction)
        attsites = list(set(attsites))
        return attsites

    def generate_products(self, site1, site2, site2_parent=None):
        """Generate product sites from recombination of two sites.

        Creates IntegraseSite objects for the products of site1 + site2
        recombination according to the reaction rules.

        Parameters
        ----------
        site1 : IntegraseSite
            First attachment site (determines product ordering).
        site2 : IntegraseSite
            Second attachment site.
        site2_parent : Polymer, optional
            Parent polymer for site2 (used in intermolecular reactions).

        Returns
        -------
        tuple of (IntegraseSite, IntegraseSite)
            Product sites at positions corresponding to site1 and site2.

        Raises
        ------
        AssertionError
            If sites have mismatched integrases or dinucleotides.
        KeyError
            If site pair is not in reaction rules.

        Notes
        -----
        Product sites inherit dinucleotides and integrase from reactants.
        Product order depends on site1 direction:

        - site1 forward: return (prod1, prod2)
        - site1 reverse: return (prod2, prod1) with swapped directions

        """
        # the sites should have the same integrase and dinucleotide, otherwise
        # it won't work
        assert site1.integrase == site2.integrase
        assert site1.integrase == self.integrase_species, (
            f"sites have integrase {site1.integrase} but should "
            f"be {self.integrase_species}"
        )
        assert site1.dinucleotide == site2.dinucleotide
        integrase = site1.integrase
        dinucleotide = site1.dinucleotide
        rxn1 = (site1.site_type, site2.site_type)
        rxn2 = (site2.site_type, site1.site_type)
        # this next part checks if these parts can even react
        if rxn1 in self.reactions:
            prod1 = self.reactions[rxn1]
        else:
            raise KeyError(
                f"{site1} not found to react with {site2} in {self.reactions}"
            )

        if rxn2 in self.reactions:
            prod2 = self.reactions[rxn2]
        else:
            raise KeyError(
                f"{site2} not found to react with {site1} in {self.reactions}"
            )

        part_prod1 = IntegraseSite(
            prod1,
            prod1,
            dinucleotide=dinucleotide,
            integrase=integrase,
            direction=site1.direction,
            integrase_binding=site1.integrase_binding,
            material_type=site1.material_type,
        )

        part_prod2 = IntegraseSite(
            prod2,
            prod2,
            dinucleotide=dinucleotide,
            integrase=integrase,
            direction=site2.direction,
            integrase_binding=site1.integrase_binding,
            material_type=site2.material_type,
        )
        if site1.direction == 'forward':
            part_prod1.position = site1.position
            part_prod1.parent = site1.parent
            part_prod2.position = site2.position
            if site2_parent is not None:
                part_prod2.parent = site2_parent
            else:
                part_prod2.parent = site2.parent
            return (part_prod1, part_prod2)
        else:
            part_prod2.direction = site1.direction
            part_prod2.position = site1.position
            part_prod2.parent = site1.parent
            part_prod1.direction = site2.direction
            part_prod1.position = site2.position
            if site2_parent is not None:
                part_prod1.parent = site2_parent
            else:
                part_prod1.parent = site2.parent

            return (part_prod2, part_prod1)

    def integrate(
        self,
        site1,
        site2,
        also_inter=True,
        force_inter=False,
        existing_dna_constructs=None,
    ):
        """Perform integrase recombination between two attachment sites.

        Executes an integration reaction between site1 and site2, creating
        new DNA constructs based on the reaction type (inversion, deletion,
        integration, or recombination). Stores transformation templates in
        the sites' linked_sites attribute.

        Parameters
        ----------
        site1 : IntegraseSite
            First attachment site (must have Construct parent).
        site2 : IntegraseSite
            Second attachment site (must have Construct parent).
        also_inter : bool, default=True
            If True and reaction is intramolecular, also generate
            intermolecular version (between two copies of same plasmid).
        force_inter : bool, default=False
            Force reaction to be treated as intermolecular even if sites
            are on same construct.
        existing_dna_constructs : list of Construct, optional
            List of previously generated constructs to check for
            duplicates.

        Returns
        -------
        list of Construct
            List of newly created DNA constructs from the integration.

        Raises
        ------
        ValueError
            If either site is not part of a Construct.

        Notes
        -----
        **Four Reaction Types:**

        1. **Inversion (intramolecular, opposite directions):**

           - Same construct, sites point opposite directions
           - Result: Region between sites is flipped
           - Circularity preserved

        2. **Deletion (intramolecular, same direction):**

           - Same construct, sites point same direction
           - Result: Two constructs - one with deleted region, one
             circular excised fragment

        3. **Integration (intermolecular, one circular):**

           - Sites on different constructs, one circular
           - Result: Single construct (circular if both were circular)

        4. **Recombination (intermolecular, both linear):**

           - Sites on two linear constructs
           - Result: Two recombinant linear constructs

        **Transformation Storage:**

        PolymerTransformation templates are stored in:

        - site1.linked_sites[(site2, intermolecular)]
        - site2.linked_sites[(site1, intermolecular)]

        These templates are used during CRN compilation to generate
        reactions and species.

        **Duplicate Detection:**

        Checks existing_dna_constructs for matches (including circular
        permutations and reversals) to avoid creating duplicates.

        Examples
        --------
        Inversion reaction:

        >>> # Two sites on same plasmid, opposite directions
        >>> int_rule.integrate(attB_site, attP_site_reversed)
        # Creates inverted plasmid

        Integration reaction:

        >>> # Sites on different plasmids
        >>> int_rule.integrate(
        ...     plasmid1_attB,
        ...     plasmid2_attP,
        ...     existing_dna_constructs=prev_constructs
        ... )
        # Creates integrated plasmid

        """
        intermolecular = True  # by default, the reaction is intermolecular
        # if one of the sites is not part of a construct then raise an error!
        integ_funcs = []
        new_dna_constructs = []  # new dna constructs made by this function!
        if existing_dna_constructs is None:
            existing_dna_constructs = []
        if not (isinstance(site1.parent, Construct)):
            raise ValueError(f"{site1} not part of a construct")
        elif not (isinstance(site2.parent, Construct)):
            raise ValueError(f"{site2} not part of a construct")
        cutpos1 = site1.position
        cutpos2 = site2.position
        # below are the references to the sites in the products
        dna_inputs = []
        if site1.parent == site2.parent and not force_inter:
            intermolecular = False
            # these sites are part of the same piece of DNA, so they are going
            # to do an intramolecular reaction
            contains_no_inter = any(
                ['no_inter' in a.attributes for a in site1.parent]
            )
            if also_inter and not contains_no_inter:
                # we should generate the intermolecular reaction also!  in
                # this case we are generating multiple integration reactions
                # at the same time i think the right thing to do is NOT return
                # it?
                new_dna_constructs += self.integrate(
                    site1,
                    site2,
                    force_inter=True,
                    existing_dna_constructs=existing_dna_constructs,
                )
            if site1.position > site2.position:
                # we reverse which position is where
                cutpos2 = site1.position
                cutpos1 = site2.position
                prod1, prod2 = self.generate_products(site2, site1)
                # integrase sites are converted into different sites according
                # to this function
            else:
                prod1, prod2 = self.generate_products(site1, site2)
                # integrase sites are converted into different sites according
                # to this function

            dna = site1.parent
            dna_inputs = [dna]
            circularity = dna.circular

            if site1.direction == site2.direction:
                if self.allow_deletion:
                    # case 2: deletion if the sites point in the same
                    # direction, then we are doing a deletion reaction
                    # direction doesn't matter so we don't need to flip
                    # anything
                    cutdna_list_parts = (
                        list(dna[:cutpos1])
                        + [[prod1, dna[cutpos1].direction]]
                        + list(dna[cutpos2 + 1 :])
                    )  # delete
                    newdna_list_parts = [
                        [prod2, dna[cutpos2].direction]
                    ] + list(dna[1 + cutpos1 : cutpos2])

                    integ_funcs += [
                        PolymerTransformation(
                            cutdna_list_parts, circular=circularity
                        ),
                        PolymerTransformation(
                            newdna_list_parts, circular=True
                        ),
                    ]
            else:
                if self.allow_inversion:
                    # case 1: inversion
                    # this means we are dealing with an inversion
                    inv_segment = []

                    [
                        [
                            a,
                        ]
                        for a in dna[cutpos1 + 1 : cutpos2][::-1]
                    ]
                    for a in dna[cutpos1 + 1 : cutpos2][::-1]:
                        inv_segment += [
                            [
                                a,
                                ['forward', 'reverse'][
                                    a.direction == 'forward'
                                ],
                            ]
                        ]
                    # the inverted segment is reversed

                    invertdna_list = (
                        list(dna[:cutpos1])
                        + [[prod1, dna[cutpos1].direction]]
                        + inv_segment
                        + [[prod2, dna[cutpos2].direction]]
                        + list(dna[cutpos2 + 1 :])
                    )

                    integ_funcs += [
                        PolymerTransformation(
                            invertdna_list, circular=circularity
                        )
                    ]
        else:
            if self.allow_integration:
                # otherwise these sites are on different pieces of DNA, so
                # they are going to combine
                dna1 = site1.parent
                dna2 = site2.parent
                circ1 = dna1.circular
                circ2 = dna2.circular
                if dna1 == dna2:
                    # this will happen if we trying to do an intermolecular
                    # reaction between two copies of the same thing
                    dna2 = dna1.__class__(
                        dna1.parts_list,
                        dna1.name + '_duplicate',
                        dna1.circular,
                    )
                dna_inputs = [dna1, dna2]
                pdict = {
                    a[1]: 'input' + str(a[0] + 1)
                    for a in enumerate(dna_inputs)
                }
                # make sure everyone is forwards
                sites = [site1, site2]
                site_halves = []
                for dna_num, site_num in zip(dna_inputs, sites):
                    if site_num.direction == 'reverse':
                        dnanum_beginning = [
                            [
                                a,
                                ['forward', 'reverse'][
                                    a.direction == 'forward'
                                ],
                            ]
                            for a in dna_num[site_num.position + 1 :][::-1]
                        ]
                        dnanum_end = [
                            [
                                a,
                                ['forward', 'reverse'][
                                    a.direction == 'forward'
                                ],
                            ]
                            for a in dna_num[: site_num.position][::-1]
                        ]
                    else:
                        dnanum_beginning = dna_num[: site_num.position]
                        dnanum_end = dna_num[site_num.position + 1 :]
                    site_halves += [
                        [list(dnanum_beginning), list(dnanum_end)]
                    ]
                dna1_halves = site_halves[0]
                dna2_halves = site_halves[1]

                if (
                    site1.direction == 'reverse'
                    and site2.direction == 'reverse'
                ):
                    prod1, prod2 = self.generate_products(
                        site2, site1, site2_parent=dna2
                    )
                else:
                    prod1, prod2 = self.generate_products(
                        site1, site2, site2_parent=dna2
                    )
                # direction of everything should be forward

                if circ2:
                    # case 3: integration
                    # in this case we are combining a circular plasmid with a
                    # circular or linear plasmid either way the result is
                    # basically the same, except the result is either linear
                    # or circular result is ONE PIECE OF DNA
                    result = (
                        dna1_halves[0]
                        + [[prod1, 'forward']]
                        + dna2_halves[1]
                        + dna2_halves[0]
                        + [[prod2, 'forward']]
                        + dna1_halves[1]
                    )
                    integ_funcs += [
                        PolymerTransformation(
                            result, circ1, parentsdict=pdict
                        )
                    ]
                elif not circ2 and circ1:
                    # if the sites are backwards just reverse everything
                    new_dna_constructs += self.integrate(
                        site2,
                        site1,
                        force_inter=force_inter,
                        existing_dna_constructs=existing_dna_constructs,
                    )
                    # the above already populates the sites, so then we don't
                    # need to
                elif not circ1 and not circ2:
                    # case 4: recombination
                    # here we are recombining two linear dnas, so two linear
                    # dnas are produced
                    result1 = (
                        dna1_halves[0] + [[prod1, 'forward']] + dna2_halves[1]
                    )
                    result2 = (
                        dna2_halves[0] + [[prod2, 'forward']] + dna1_halves[1]
                    )
                    integ_funcs += [
                        PolymerTransformation(result1, parentsdict=pdict),
                        PolymerTransformation(result2, parentsdict=pdict),
                    ]
        if len(integ_funcs) > 0:
            for integ_func in integ_funcs:
                # generate new dna constructs and check if we already made
                # them before
                created_dna = integ_func.create_polymer(
                    [site1.parent, site2.parent]
                )
                output = Integrase_Enumerator.find_dna_construct(
                    created_dna, existing_dna_constructs + new_dna_constructs
                )
                if output is not None:
                    # if the construct has already been made, then fix the
                    # integ_func so it outputs the right thing
                    integ_func.renumber_output(output[1])
                else:
                    # otherwise, add it to the list
                    new_dna_constructs += [created_dna]
            # link the two sites and give them the adjusted integ_funcs
            site1.linked_sites[(site2, intermolecular)] = [integ_funcs, []]
            site2.linked_sites[(site1, intermolecular)] = [
                [a.reversed() for a in integ_funcs],
                [],
            ]
        # the return value of this function is used mostly only for generating
        # constructs
        return new_dna_constructs


class IntegraseEnumerator(GlobalComponentEnumerator):
    """Global enumerator for integrase-mediated DNA recombination products.

    An `IntegraseEnumerator` systematically enumerates all possible DNA
    constructs that can result from integrase-mediated recombination
    reactions. Examines all components for integrase attachment sites and
    generates products for all allowed site pairs.

    Parameters
    ----------
    name : str
        Name identifier for the enumerator.
    int_mechanisms : dict, optional
        Dictionary mapping integrase names (str) to IntegraseRule objects.
        Default: {'int1': IntegraseRule()}

    Attributes
    ----------
    int_mechanisms : dict
        Dictionary of integrase rules.

    See Also
    --------
    GlobalComponentEnumerator : Base class for global enumeration.
    IntegraseRule : Defines integrase recombination rules.
    IntegraseSite : DNA part for attachment sites.
    PolymerTransformation : Template for DNA rearrangements.

    Notes
    -----
    **Enumeration Process:**

    1. Identify all integrase attachment sites in components
    2. Group sites by integrase type
    3. For each integrase:

       a. Find all valid site pairs (from reactive_sites)
       b. Check if pair can react (reaction_allowed)
       c. Perform integration to generate products
       d. Store transformation templates in sites

    4. Return list of new DNA constructs

    **Global Enumeration:**

    This is a global enumerator because integrase reactions can occur
    between sites on different constructs (intermolecular reactions).
    Access to all components is necessary.

    **Integrase Types Supported:**

    - Serine integrases (attB/attP --> attL/attR)
    - Tyrosine recombinases (Cre, Flp with homotypic sites)
    - Invertases (inversion only)
    - Resolvases (deletion only)
    - Custom integrase rules

    **Duplicate Handling:**

    Uses `find_dna_construct` to detect duplicates including circular
    permutations and reversals, preventing redundant construct
    generation.

    Examples
    --------
    Create an integrase enumerator:

    >>> phi_c31 = bcp.IntegraseRule(
    ...     name='PhiC31',
    ...     reactions={
    ...         ('attB', 'attP'): 'attL',
    ...         ('attP', 'attB'): 'attR'
    ...     }
    ... )
    >>> enumerator = bcp.IntegraseEnumerator(
    ...     name='integrase_enum',
    ...     int_mechanisms={'PhiC31': phi_c31}
    ... )

    Use in a mixture:

    >>> mixture = bcp.Mixture(
    ...     components=[plasmid1, plasmid2],
    ...     global_component_enumerators=[enumerator]
    ... )
    >>> # Enumerator automatically called during compilation
    >>> crn = mixture.compile_crn()

    Manual enumeration:

    >>> constructs = [plasmid_with_attB, plasmid_with_attP]
    >>> new_constructs = enumerator.enumerate_components(constructs)
    >>> # new_constructs contains integrated plasmids

    """

    def __init__(self, name: str, int_mechanisms=None):
        if int_mechanisms is None:
            int_mechanisms = {'int1': IntegraseRule()}
        self.int_mechanisms = int_mechanisms
        GlobalComponentEnumerator.__init__(self, name=name)

    def list_integrase(self, construct):
        """List all integrase attachment sites in a construct.

        Parameters
        ----------
        construct : Construct
            DNA construct to examine.

        Returns
        -------
        dict
            Dictionary mapping integrase names (str) to lists of
            IntegraseSite objects.

        """
        int_dict = {}
        for part in construct.parts_list:
            if isinstance(part, IntegraseSite) and part.integrase is not None:
                part_integrase = part.integrase.name
                if part_integrase in int_dict:
                    int_dict.update(
                        {part_integrase: int_dict[part_integrase] + [part]}
                    )
                else:
                    int_dict[part_integrase] = [part]
        return int_dict

    def reset(self, components=None, **kwargs):
        """Reset linked_sites attribute in all attachment sites.

        Clears stored integration reactions from all integrase sites in
        components, preparing for fresh enumeration.

        Parameters
        ----------
        components : list of Component
            Components containing integrase sites to reset.
        **kwargs
            Additional keyword arguments (unused).

        Notes
        -----
        Called at the start of enumeration to clear previous integration
        data.

        """
        for component in components:
            if hasattr(component, 'parts_list'):
                for part in component:
                    if hasattr(part, 'linked_sites'):
                        part.linked_sites = {}

    @classmethod
    def find_dna_construct(
        cls, construct: Construct, conlist: List[Construct]
    ):
        """Find matching construct in list (handles permutations/reversals).

        Searches for a construct equivalent to the input, accounting for
        circular permutations and reversals.

        Parameters
        ----------
        construct : Construct
            Construct to find.
        conlist : list of Construct
            List of constructs to search.

        Returns
        -------
        tuple of (Construct, callable) or None
            If found: (matched_construct, index_function), where
            index_function maps old indexes to (new_index, direction).
            If not found: None.

        Raises
        ------
        KeyError
            If construct matches multiple constructs in list (should not
            happen with proper generation order).

        Notes
        -----
        **Matching Logic:**

        For circular constructs:

        - Try all circular permutations
        - For each permutation, try forward and reverse orientations

        For linear constructs:

        - Try forward orientation
        - Try reverse orientation

        Uses directionless_hash for fast initial filtering before
        detailed species comparison.

        """
        matched_construct = None
        for other_construct in conlist:
            if not isinstance(other_construct, Construct):
                continue
            if (
                construct.directionless_hash
                == other_construct.directionless_hash
            ):
                if matched_construct is not None:
                    # a construct must not match two constructs, since they
                    # are generated and checked in order
                    raise KeyError(
                        f"{construct} matches with {matched_construct} "
                        f"but also {other_construct}"
                    )
                matched_construct = other_construct
        if matched_construct is None:
            return None
        other_construct = matched_construct
        other_indexes = list(range(len(other_construct)))
        if construct.circular:
            for pos, _ in enumerate(construct):
                cp_construct = construct.get_circularly_permuted(pos)

                if (
                    cp_construct.get_species()
                    == other_construct.get_species()
                ):
                    return other_construct, (
                        lambda a: (
                            (other_indexes[pos:] + other_indexes[:pos])[a],
                            'f',
                        )
                    )
                elif (
                    cp_construct.get_reversed().get_species()
                    == other_construct.get_species()
                ):
                    return other_construct, (
                        lambda a: (
                            (other_indexes[pos:] + other_indexes[:pos])[::-1][
                                a
                            ],
                            'r',
                        )
                    )
        else:
            if construct.get_species() == other_construct.get_species():
                return other_construct, (lambda a: (a, 'f'))
            elif (
                construct.get_reversed().get_species()
                == other_construct.get_species()
            ):
                return other_construct, (
                    lambda a: (other_indexes[::-1][a], 'r')
                )

    def enumerate_components(
        self, components=None, previously_enumerated=None, **kwargs
    ):
        """Enumerate all possible integrase-mediated DNA configurations.

        Systematically generates all DNA constructs that can result from
        integrase recombination between attachment sites in the input
        components.

        Parameters
        ----------
        components : list of Component, optional
            List of components to enumerate. Only DNAconstruct objects
            are processed.
        previously_enumerated : list of Component, optional
            List of components already enumerated (used for duplicate
            detection).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Construct
            List of newly created DNA constructs from all allowed
            integrase reactions.

        Notes
        -----
        **Enumeration Algorithm:**

        1. Extract all DNAconstruct components
        2. List all integrase sites by integrase type
        3. For each integrase in int_mechanisms:

           a. Get all attachment sites for that integrase
           b. Find reactive site types from IntegraseRule
           c. Generate all site pairs (combinations)
           d. For each valid pair:

              - Check if reaction_allowed
              - Perform `integrate` to generate products
              - Add new constructs to output list

        4. Return all newly generated constructs

        **Reaction Types Generated:**

        Depending on IntegraseRule settings:

        - Inversions (same construct, opposite directions)
        - Deletions (same construct, same direction)
        - Integrations (different constructs, at least one circular)
        - Recombinations (two linear constructs)

        **Duplicate Prevention:**

        The `integrate` method checks existing_dna_constructs (includes
        both previously_enumerated and newly created constructs) to avoid
        generating duplicates.

        Examples
        --------
        Enumerate integration products:

        >>> enumerator = bcp.IntegraseEnumerator(
        ...     name='enum',
        ...     int_mechanisms={'PhiC31': phi_c31_rule}
        ... )
        >>> plasmids = [donor_plasmid, target_plasmid]
        >>> products = enumerator.enumerate_components(plasmids)
        >>> # products contains integrated plasmids

        """
        if previously_enumerated is None:
            previously_enumerated = []
        construct_list = []
        if previously_enumerated is None:
            previously_enumerated = []
        for component in components:
            if isinstance(component, DNAconstruct):
                construct_list += [component]
        int_dict = {}
        for construct in construct_list:
            # list each integrase that exists and which sites they react with

            con_dict = self.list_integrase(construct)

            int_dict = combine_dictionaries(int_dict, con_dict)
        constructlist = []
        for integrase in int_dict:
            if integrase in self.int_mechanisms:
                int_mech = self.int_mechanisms[integrase]
                # now, going through each one, generate the reactions and
                # species that arise
                attsites = int_dict[integrase]
                # but now we need to know what kind of integrase reactions are
                # possible
                reactive_sites = int_mech.reactive_sites()
                attcombos = [
                    a
                    for a in it.combinations(attsites, 2)
                    if (
                        (a[0].site_type in reactive_sites)
                        and (a[1].site_type in reactive_sites)
                    )
                ]
                for combo in attcombos:
                    # first question: is this combo legal?
                    if int_mech.reaction_allowed(combo[0], combo[1]):
                        # this means the reaction can exist
                        # integrate now
                        new_dnas = int_mech.integrate(
                            combo[0],
                            combo[1],
                            existing_dna_constructs=(
                                previously_enumerated + constructlist
                            ),
                        )
                        constructlist += new_dnas

        return constructlist


# Legacy names
Polymer_transformation = PolymerTransformation
Integrase_Enumerator = IntegraseEnumerator
