#  Copyright (c) 2026, BioCRNpyler Developers. All rights reserved.
#  See LICENSE file in the project root directory for details.

"""Consistency checks across the default parameter files.

The same parameter is often declared in several files: a mechanism default
in `mechanisms/*_parameters.tsv` and again in each `mixtures/*_parameters.tsv`
that wants a different value.  Nothing keeps those copies in step, so a value
updated in one file can silently go stale in another.  The tests here pin down
which divergences are deliberate and check that everything else agrees.
"""

import glob
import io
import os
from collections import defaultdict

import pytest

import biocrnpyler as bcp

# Parameters that are deliberately different between files, with the reason.
# Anything else that differs is treated as drift and fails the first test.
INTENTIONAL_DIVERGENCES = {
    # Cells lose RNA to growth dilution as well as degradation, so the
    # one step expression rate is lower than in a cell free extract.
    ('gene_expression', '', 'kexpress'),
    # Degradation is faster in vivo than in extract, so both the
    # Michaelis-Menten constants and the first order rate differ.
    ('rna_degradation', '', 'kdil'),
    ('rna_degradation_mm', '', 'kdeg'),
    ('rna_degradation_mm', '', 'ku'),
    # Machinery and metabolite levels are properties of the system being
    # modelled: whole cells, TX-TL extract, and reconstituted PURE each
    # have their own measurements.
    ('initial concentration', '', 'ATP'),
    ('initial concentration', '', 'NTPs'),
    ('initial concentration', '', 'RNAP'),
    ('initial concentration', '', 'RNase'),
    ('initial concentration', '', 'Ribo'),
}

# Files holding mechanism defaults rather than the settings of one mixture.
MECHANISM_FILES = {
    'binding_parameters.tsv',
    'enzyme_parameters.tsv',
    'transport_parameters.tsv',
    'txtl_parameters.tsv',
}

# Mixtures are exercised with a single constitutively expressed gene.
MIXTURES = [
    'ExpressionExtract',
    'SimpleTxTlExtract',
    'TxTlExtract',
    'EnergyTxTlExtract',
    'PURE',
    'ExpressionDilutionMixture',
    'SimpleTxTlDilutionMixture',
    'TxTlDilutionMixture',
]


def parameter_files():
    root = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    )
    pattern = os.path.join(root, 'biocrnpyler', '**', '*.tsv')
    return sorted(glob.glob(pattern, recursive=True))


def parameters_by_key():
    """Map (mechanism, part_id, param_name) to {file name: value}."""
    table = defaultdict(dict)
    for path in parameter_files():
        for line in io.open(path, encoding='utf-8'):
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 4 or fields[0] in ('mechanism', 'mechanism_id'):
                continue
            key = (fields[0], fields[1], fields[2])
            table[key][os.path.basename(path)] = fields[3]
    return table


def make_gene():
    return bcp.DNAassembly(
        name='gfp', promoter='pconst', rbs='rbs_strong', protein='GFP'
    )


def concentration(entry):
    """Initial concentrations come back as Parameters or as plain numbers."""
    return float(getattr(entry, 'value', entry))


def test_parameter_files_are_readable():
    # Guards the file based tests, which quietly pass if nothing is parsed
    table = parameters_by_key()
    assert len(parameter_files()) >= 5
    assert len(table) > 50


def test_no_unintended_divergence():
    """Duplicated parameters agree unless the difference is deliberate."""
    drifted = {}
    for key, values in parameters_by_key().items():
        if (
            len(set(values.values())) > 1
            and key not in INTENTIONAL_DIVERGENCES
        ):
            drifted[key] = values
    assert not drifted, (
        "These parameters are declared in several files with different "
        "values. Either bring them back into agreement, or add the key to "
        f"INTENTIONAL_DIVERGENCES with a reason: {drifted}"
    )


def test_mechanism_defaults_are_not_stale():
    """A mechanism default must be a value some mixture actually uses.

    Mixtures on the allowlist are expected to disagree with each other, so
    that test alone cannot tell a deliberate difference from a mechanism
    default left behind when the mixtures moved on.  Requiring the default
    to match at least one mixture catches the stale copy without
    constraining the mixtures.
    """
    stale = {}
    for key, values in parameters_by_key().items():
        defaults = {f: v for f, v in values.items() if f in MECHANISM_FILES}
        used = {v for f, v in values.items() if f not in MECHANISM_FILES}
        if not defaults or not used:
            continue
        for name, value in defaults.items():
            if value not in used:
                stale[(key, name)] = (value, sorted(used))
    assert not stale, (
        "These mechanism defaults do not match any value a mixture uses, "
        "which usually means the default was missed when the mixture files "
        f"were updated: {stale}"
    )


def test_allowlisted_parameters_are_set_explicitly():
    """A mixture must set any allowlisted parameter its mechanisms use.

    Values on the allowlist are expected to differ per mixture, so relying
    on a default from another file would silently give the wrong one.
    """
    missing = []
    for name in MIXTURES:
        mixture = getattr(bcp, name)(name='m', components=[make_gene()])
        mechanisms = list(mixture.mechanisms.values()) + list(
            mixture.global_mechanisms.values()
        )
        used = {getattr(m, 'name', None) for m in mechanisms}
        declared = {
            (key.mechanism, key.part_id or '', key.name)
            for key in mixture.parameter_database.parameters
        }
        for key in sorted(INTENTIONAL_DIVERGENCES):
            if key[0] in used and key not in declared:
                missing.append((name, key))
    assert not missing, (
        "These mixtures use a mechanism whose parameter is expected to "
        "differ per mixture, but do not set it in their own parameter "
        f"file: {missing}"
    )


@pytest.mark.parametrize('name', MIXTURES)
def test_mixture_compiles_and_expresses(name):
    """Every mixture builds a CRN that actually makes protein.

    A missing parameter can leave a species at zero rather than raising,
    so compiling is not on its own evidence that a mixture works.
    """
    crn = getattr(bcp, name)(name='m', components=[make_gene()]).compile_crn()
    assert crn.reactions, f'{name} compiled with no reactions'

    species = {str(s) for s in crn.species}
    assert 'protein_GFP' in species

    initial = crn.initial_concentration_dict
    machinery = {
        str(s): concentration(initial.get(s, 0))
        for s in crn.species
        if ('RNAP' in str(s) or 'Ribo' in str(s)) and 'complex' not in str(s)
    }
    assert all(v > 0 for v in machinery.values()), (
        f'{name} leaves machinery at zero concentration: {machinery}'
    )
