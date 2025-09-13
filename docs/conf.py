# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

project = 'BioCRNPyler'
copyright = '2025, Build-a-Cell'
author = 'William Poole, Ayush Pandey, Andrey Shur, Zoltan Tuza, Richard M. Murray'

# Import the package
import biocrnpyler

# The short X.Y version
version = '1.3'
# The full version, including alpha/beta/rc tags
release = '1.3.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc', 'sphinx.ext.linkcode', 'sphinx.ext.doctest',
    'sphinx.ext.autosummary', 'sphinx_copybutton', 'sphinx_toggleprompt',
    'nbsphinx', 'nbsphinx_link',
    'recommonmark'
]

source_suffix = ['.rst']

# scan documents for autosummary directives and generate stub pages for each.
autosummary_generate = True

# list of autodoc directive flags that should be automatically applied
# to all autodoc directives.
autodoc_default_options = {
    'members': True,
    'inherited-members': True,
    'exclude-members': '__init__, __weakref__, __repr__, __str__'
}

# For classes, include both the class docstring and the init docstring
autoclass_content = 'both'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Set the default role to render items in backticks as code
default_role = 'py:obj'

# Skip prompts when using copy button
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']
html_css_files = ['css/custom.css']

# -----------------------------------------------------------------------------
# Source code links (from numpy)
# -----------------------------------------------------------------------------

import inspect
from os.path import relpath, dirname

def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    # print(f"{domain=}, {info=}")
    if domain != 'py':
        # print("  domain != 'py'")
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        # print("  submod is None")
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except Exception:
            # print("  getattr Exception")
            return None

    # strip decorators, which would resolve to the source of the decorator
    # possibly an upstream bug in getsourcefile, bpo-1764286
    try:
        unwrap = inspect.unwrap
    except AttributeError:
        pass
    else:
        obj = unwrap(obj)

    # Get the filename for the function
    try:
        fn = inspect.getsourcefile(obj)
    except Exception:
        fn = None
    if not fn:
        # print("  not fn")
        return None

    # Ignore re-exports as their source files are not within the numpy repo
    module = inspect.getmodule(obj)
    if module is not None and not module.__name__.startswith("biocrnpyler"):
        # print("module is not None but doesnt start with biocrnpyler")
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except Exception:
        lineno = None

    fn = relpath(fn, start=dirname(biocrnpyler.__file__))

    if lineno:
        linespec = "#L%d-L%d" % (lineno, lineno + len(source) - 1)
    else:
        linespec = ""

    base_url = "https://github.com/BuildACell/BioCRNPyler/blob/"
    if release != version:      # development release
        # TODO: replace 'refactor-modules' with 'master'
        # print("  --> ", base_url + "refactor-modules/control/%s%s" % (fn, linespec))
        return base_url + "refactor-modules/biocrnpyler/%s%s" % (fn, linespec)
    else:                       # specific version
        return base_url + "%s/biocrnpyler/%s%s" % (version, fn, linespec)

# -- Options for doctest ----------------------------------------------

# Import biocrnpyler as bcp
doctest_global_setup = """
import numpy as np
import biocrnpyler as bcp
"""
    
