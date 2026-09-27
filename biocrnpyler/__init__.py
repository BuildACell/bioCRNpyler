# All core classes
# Library of all components
from .components import *
from .core import *

# Library of all mechanisms
from .mechanisms import *

# Library of all mixtures
from .mixtures import *

# All utilities
from .utils import *
from .utils.plotting import (  # noqa: F401
    plot_all_species_containing,
    plot_gene_expression_data,
)

try:
    from ._version import version as __version__
except Exception:
    __version__ = "0+unknown"

# Mapping of deprecated legacy mechanism names to their updated class names
_LEGACY_NAMES = {
    'Simple_Diffusion': 'Diffusion_Simple',
    'Simple_Transport': 'Diffusion_Facilitated_Channel',
    'Facilitated_Transport_MM': 'Diffusion_Facilitated_Carrier',
    'Primary_Active_Transport_MM': 'Transport_PrimaryActive_ABCexporter',
    'Membrane_Signaling_Pathway_MM': 'Sensor_TwoComponentSystem',
    'Membrane_Protein_Integration': 'Integration_MembraneProtein',
}


def __getattr__(name):
    if name in _LEGACY_NAMES:
        new = _LEGACY_NAMES[name]

        if name == 'Facilitated_Transport_MM':
            warn(
                f"'{name}' is deprecated; use '{new}'. Note: New transport" \
                "mechanisms are also available, including" \
                "'Transport_SecondaryActive_Symporter' and" \
                "'Transport_SecondaryActive_Antiporter'.",
                DeprecationWarning,
                stacklevel=2,
            )
        else:
            warn(
                f"'{name}' is deprecated; use '{new}'",
                DeprecationWarning,
                stacklevel=2,
            )

        return globals()[new]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")