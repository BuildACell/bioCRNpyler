"""BioCRNpyler mixture library.

A mixture in BioCRNpyler defines the *context* in which components are
compiled into a chemical reaction network (CRN). A mixture ties
together components, mechanisms, and parameters by specifying *which*
Mechanisms are available, *which* components are present, and *what*
parameters to use.

"""

from .cell import *
from .extract import *
from .pure import *
