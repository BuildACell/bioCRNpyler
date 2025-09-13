"""BioCRNpyler mixture library.

A Mixture in BioCRNpyler defines the *context* in which Components are
compiled into a chemical reaction network (CRN). A mixture ties
together components, mechanisms, and parameters by specifying *which*
Mechanisms are available, *which* Components are present, and *what*
parameters to use.

"""

from .cell import *
from .extract import *
