#  Copyright (c) 2025, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

import pytest

from biocrnpyler import PURE, BasicPURE, DNAassembly


def make_gene():
    return DNAassembly('GFP', promoter='pconst', rbs='RBS', protein='GFP')


def species_names(mixture):
    return {str(species) for species in mixture.compile_crn().species}


def mechanism_names(mixture):
    return {
        key: type(mixture.mechanisms[key]).__name__
        for key in ('transcription', 'translation')
    }


def test_pure_without_machinery():
    # With everything turned off, expression is a simple two step process
    mixture = PURE(
        name='simple',
        components=[make_gene()],
        include_machinery=False,
        include_resources=False,
        include_fuel=False,
    )
    assert mechanism_names(mixture) == {
        'transcription': 'SimpleTranscription',
        'translation': 'SimpleTranslation',
    }

    # No machinery, resources, or energy carriers are created
    assert mixture.atp is None
    assert mixture.adp is None
    assert mixture.fuel is None
    assert species_names(mixture) == {'dna_GFP', 'rna_GFP', 'protein_GFP'}


def test_pure_with_machinery():
    # Machinery adds RNAP and ribosome, and switches to MM kinetics
    mixture = PURE(
        name='machinery',
        components=[make_gene()],
        include_machinery=True,
        include_resources=False,
        include_fuel=False,
    )
    assert mechanism_names(mixture) == {
        'transcription': 'Transcription_MM',
        'translation': 'Translation_MM',
    }

    names = species_names(mixture)
    assert 'protein_RNAP' in names
    assert 'protein_Ribo' in names

    # Resources and energy carriers are still absent
    assert 'metabolite_NTPs' not in names
    assert 'metabolite_ATP' not in names
    assert mixture.atp is None


def test_pure_with_resources():
    # Resources add NTPs and amino acids, and switch to energy mechanisms
    mixture = PURE(
        name='resources',
        components=[make_gene()],
        include_machinery=True,
        include_resources=True,
        include_fuel=False,
    )
    assert mechanism_names(mixture) == {
        'transcription': 'Energy_Transcription_MM',
        'translation': 'Energy_Translation_MM',
    }

    names = species_names(mixture)
    assert 'metabolite_NTPs' in names
    assert 'metabolite_AAs' in names

    # Energy carriers require include_fuel
    assert 'metabolite_ATP' not in names
    assert mixture.atp is None


def test_pure_with_fuel():
    # Fuel adds the energy carriers, which are regenerated from the fuel
    mixture = PURE(
        name='fuel',
        components=[make_gene()],
        include_machinery=True,
        include_resources=True,
        include_fuel=True,
    )
    names = species_names(mixture)
    assert 'metabolite_ATP' in names
    assert 'metabolite_ADP' in names
    assert 'metabolite_Fuel_CP' in names

    assert mixture.atp is not None
    assert mixture.adp is not None
    assert mixture.fuel is not None


def test_pure_switches_are_nested():
    # Each switch requires the ones before it
    with pytest.raises(ValueError, match='include_fuel requires'):
        PURE(
            name='bad',
            components=[make_gene()],
            include_machinery=True,
            include_resources=False,
            include_fuel=True,
        )

    with pytest.raises(ValueError, match='include_resources requires'):
        PURE(
            name='bad',
            components=[make_gene()],
            include_machinery=False,
            include_resources=True,
            include_fuel=False,
        )


def test_pure_species_names_can_be_changed():
    # Each species name argument is used to name the species that is created
    mixture = PURE(
        name='renamed',
        components=[make_gene()],
        rnap='T7',
        ribosome='Ribosome',
        ntps='NTP',
        amino_acids='AminoAcids',
        adp='SpentEnergy',
        fuel='CreatinePhosphate',
    )
    names = species_names(mixture)
    for expected in [
        'protein_T7',
        'protein_Ribosome',
        'metabolite_NTP',
        'metabolite_AminoAcids',
        'metabolite_SpentEnergy',
        'metabolite_CreatinePhosphate',
    ]:
        assert expected in names

    # The defaults should no longer be present
    assert 'protein_RNAP' not in names
    assert 'metabolite_NTPs' not in names


@pytest.mark.xfail(
    reason="ATP production and degradation parameters in "
    "pure_parameters.tsv are keyed to the part_ids ATP_production and "
    "ATP_degradation, so renaming atp orphans them",
    raises=ValueError,
    strict=True,
)
def test_pure_atp_name_can_be_changed():
    PURE(
        name='renamed', components=[make_gene()], atp='Energy'
    ).compile_crn()


def test_basic_pure_is_deprecated():
    # BasicPURE warns, but still builds the same model as PURE
    with pytest.warns(DeprecationWarning, match='BasicPURE is deprecated'):
        legacy = BasicPURE(name='legacy', components=[make_gene()])

    current = PURE(
        name='legacy',
        components=[make_gene()],
        include_machinery=True,
        include_resources=True,
        include_fuel=True,
    )
    assert species_names(legacy) == species_names(current)


def test_basic_pure_forwards_species_names():
    # Species names given to BasicPURE are passed through to PURE
    with pytest.warns(DeprecationWarning):
        mixture = BasicPURE(
            name='legacy', components=[make_gene()], ribosome='Ribosome'
        )
    names = species_names(mixture)
    assert 'protein_Ribosome' in names
    assert 'protein_Ribo' not in names
