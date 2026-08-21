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
        include_energy=False,
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
        include_energy=False,
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
        include_energy=False,
    )
    assert mechanism_names(mixture) == {
        'transcription': 'Energy_Transcription_MM',
        'translation': 'Energy_Translation_MM',
    }

    names = species_names(mixture)
    assert 'metabolite_NTPs' in names
    assert 'metabolite_AAs' in names

    # Energy carriers require include_energy
    assert 'metabolite_ATP' not in names
    assert mixture.atp is None


def test_pure_with_fuel():
    # Fuel adds the energy carriers, which are regenerated from the fuel
    mixture = PURE(
        name='fuel',
        components=[make_gene()],
        include_machinery=True,
        include_resources=True,
        include_energy=True,
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
    with pytest.raises(ValueError, match='include_energy requires'):
        PURE(
            name='bad',
            components=[make_gene()],
            include_machinery=True,
            include_resources=False,
            include_energy=True,
        )

    with pytest.raises(ValueError, match='include_resources requires'):
        PURE(
            name='bad',
            components=[make_gene()],
            include_machinery=False,
            include_resources=True,
            include_energy=False,
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


def test_pure_atp_name_requires_its_own_parameters():
    # Parameters are keyed by species name, so renaming the energy carrier
    # orphans the ATP_production and ATP_degradation entries shipped with
    # the mixture.  That should be an error rather than a model that is
    # quietly missing its energy pathway.
    with pytest.raises(ValueError, match='Energy_production'):
        PURE(
            name='renamed', components=[make_gene()], atp='Energy'
        ).compile_crn()

    # Supplying replacements for the new name resolves it
    mixture = PURE(
        name='renamed',
        components=[make_gene()],
        atp='Energy',
        parameters={
            ('one_step_pathway', 'Energy_production', 'k'): 0.02,
            ('one_step_pathway', 'Energy_degradation', 'k'): 0.0000177,
            ('initial concentration', None, 'Energy'): 1000,
        },
    )
    names = species_names(mixture)
    assert 'metabolite_Energy' in names
    assert 'metabolite_ATP' not in names

    # The initial concentration has to be given for the new name too
    initial = {
        str(key): getattr(value, 'value', value)
        for key, value in (
            mixture.compile_crn().initial_concentration_dict or {}
        ).items()
    }
    assert initial['metabolite_Energy'] == 1000


def test_pure_without_regeneration():
    # fuel=None keeps ATP consumption but drops the regeneration pathway
    mixture = PURE(name='noregen', components=[make_gene()], fuel=None)

    names = species_names(mixture)
    assert 'metabolite_ATP' in names
    assert 'metabolite_ADP' in names
    assert not any('Fuel' in name for name in names)

    assert mixture.atp is not None
    assert mixture.adp is not None
    assert mixture.fuel is None

    # ATP keeps its initial concentration, and gains no pathway reactions
    # of its own, so the ATP_production parameters are never needed
    crn = mixture.compile_crn()
    initial = {
        str(key): getattr(value, 'value', value)
        for key, value in (crn.initial_concentration_dict or {}).items()
    }
    assert initial['metabolite_ATP'] == 1000


def test_basic_pure_is_deprecated():
    # BasicPURE warns, and builds PURE without the regeneration pathway
    with pytest.warns(DeprecationWarning, match='BasicPURE is deprecated'):
        legacy = BasicPURE(name='legacy', components=[make_gene()])

    current = PURE(
        name='legacy',
        components=[make_gene()],
        include_machinery=True,
        include_resources=True,
        include_energy=True,
        fuel=None,
    )
    assert species_names(legacy) == species_names(current)

    # and so has no fuel species, unlike PURE's default
    assert not any('Fuel' in name for name in species_names(legacy))


def test_basic_pure_fuel_names_the_carrier():
    # In BasicPURE, fuel names the energy carrier, as it did before PURE
    # was introduced; it is an alias for PURE's atp argument
    with pytest.warns(DeprecationWarning):
        legacy = BasicPURE(
            name='legacy', components=[make_gene()], fuel='ATP'
        )
    assert 'metabolite_ATP' in species_names(legacy)

    with pytest.warns(DeprecationWarning):
        by_fuel = BasicPURE(
            name='a', components=[make_gene()], fuel='GTP',
            parameters={('initial concentration', None, 'GTP'): 1000},
        )
    with pytest.warns(DeprecationWarning):
        by_atp = BasicPURE(
            name='a', components=[make_gene()], atp='GTP',
            parameters={('initial concentration', None, 'GTP'): 1000},
        )
    assert species_names(by_fuel) == species_names(by_atp)

    # Giving both is an error rather than one silently winning
    with pytest.warns(DeprecationWarning):
        with pytest.raises(ValueError, match='not both'):
            BasicPURE(
                name='a', components=[make_gene()], fuel='X', atp='Y'
            )


def test_basic_pure_forwards_species_names():
    # Species names given to BasicPURE are passed through to PURE
    with pytest.warns(DeprecationWarning):
        mixture = BasicPURE(
            name='legacy', components=[make_gene()], ribosome='Ribosome'
        )
    names = species_names(mixture)
    assert 'protein_Ribosome' in names
    assert 'protein_Ribo' not in names
