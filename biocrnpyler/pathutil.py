# pathutil.py - path utilities
# RMM, 16 Aug 2018
#
# This file contains some utility functions for manipulating and using
# paths for finding models and configuration files.
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from importlib import import_module

# Load a model from a file
def load_model(prefix, name, length):
    # Look to see if we have a model for this component
    #! Expand this to look in other locations
    model = None
    try:
        module = import_module("txtl.components.%s_%s" %
                               (prefix.lower(), name.lower()))
        model = eval("module.%s_%s('name=%s', length=%d)" %
                     (prefix.lower(), name.lower(), name, length))
        print("Warning: Eval Statements for Class Construction are being "
              "Depricated.")

    except ModuleNotFoundError as error:
        print(error)
        print("couldn't find component %s_%s" % (prefix, name))
    return model
