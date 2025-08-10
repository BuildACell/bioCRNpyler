from biocrnpyler.components import Protein, DNAassembly
from biocrnpyler.components import RepressiblePromoter
from biocrnpyler.core import Compartment
from biocrnpyler.mixtures import SimpleTxTlExtract


def test_compartment_dna_assembly():
    """Test compartment assignment for DNA assembly in a CRN."""

    # Create compartments for the inside and outside of the detector cell
    detector_internal = Compartment(name="DetectorInternal")

    # Create the protein for TetR
    tetr = Protein("TetR", compartment=detector_internal)

    # Create a DNA construct for a GFP reporter
    ptet = RepressiblePromoter(
        name="ptet", repressor=tetr, compartment=detector_internal
    )
    
    dna_gfp = DNAassembly(
        name="GFP",
        promoter=ptet,
        rbs="RBS_standard",
        protein="GFP",
        compartment=detector_internal,
    )

    parameter_dict = {
        "ktx": 1.0,  # Transcription rate
        "ktl": 1.0,  # Translation rate
        "kdil": 0.1,  # Dilution rate
        "kleak": 0.01, # Leak rate
        "K": 10, # Half-activation value
        "k": 1, # activation rate
        "n": 2, # Hill coefficient
    }

    # Create mixture
    mixture = SimpleTxTlExtract(
        "DetectorCell",
        components=[tetr, dna_gfp],
        parameters=parameter_dict,
    )
    # Create the chemical reaction network
    crn = mixture.compile_crn()
    # assert that compartments are correctly set
    assert tetr.compartment == detector_internal
    assert ptet.compartment == detector_internal
    assert dna_gfp.compartment == detector_internal
    # assert that the compartment names match
    assert tetr.compartment.name == detector_internal.name
    assert ptet.compartment.name == detector_internal.name
    assert dna_gfp.compartment.name == detector_internal.name
    # assert that compartments are correct after compilation
    for species in crn.species:
        assert species.compartment.name == detector_internal.name
    
