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
import inspect
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# Use the readthedocs.org theme if installed
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    try:
        import sphinx_rtd_theme  # noqa: F401

        html_theme = 'sphinx_rtd_theme'
    except ImportError:
        html_theme = 'default'

# -- Project information -----------------------------------------------------

project = 'BioCRNPyler'
copyright = '2025, Build-a-Cell'
author = """Zoila Jurado, William Poole, Ayush Pandey, Andrey Shur, Zoltan
    Tuza, Richard M. Murray"""

# Import the package
import biocrnpyler
from setuptools_scm import get_version

release = get_version(root="..", relative_to=__file__)

# Short X.Y
version = ".".join(release.split(".", 2)[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.linkcode',
    'sphinx.ext.doctest',
    'sphinx.ext.imgmath',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'sphinx_toggleprompt',
    'nbsphinx',
    'nbsphinx_link',
    'recommonmark',
    'numpydoc',
]

source_suffix = ['.rst']

# scan documents for autosummary directives and generate stub pages for each.
autosummary_generate = True

# list of autodoc directive flags that should be automatically applied
# to all autodoc directives.
autodoc_default_options = {
    #    'members': True,
    #    'inherited-members': True,
    #    'special-members': True,
    'exclude-members': '__init__, __weakref__, __repr__, __str__, __hash__',
}

# For classes, include both the class docstring and the init docstring
autoclass_content = 'class'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# This config value contains the locations and names of other projects that
# should be linked to in this documentation.
intersphinx_mapping = {
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'python': ('https://docs.python.org/3/', None),
}

# Don't generate external links to (local) keywords
intersphinx_disabled_reftypes = ["py:keyword"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Set the default role to render items in backticks as code
default_role = 'py:obj'

# Align inline math with text
imgmath_use_preview = True

# Skip prompts when using copy button
copybutton_prompt_text = r'>>> |\.\.\. '
copybutton_prompt_is_regexp = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".

html_static_path = ['_static']
html_css_files = ['css/custom.css']

# Don't automatically show all members of class in Methods & Attributes section
numpydoc_show_class_members = False

# Don't create a Sphinx TOC for the lists of class methods and attributes
numpydoc_class_members_toctree = False

# Leave Attributes documentation right after Parameters
napoleon_use_ivar = False
napoleon_custom_sections = [
    ('Attributes', 'params_style'),
]

# Aliases to allow objects to avoid including module names
napoleon_use_param = True
napoleon_preprocess_types = True  # convert refs in types to std form
napoleon_type_aliases = dict()
for name, obj in inspect.getmembers(biocrnpyler):
    if inspect.isclass(obj):
        parent_path = os.path.dirname(obj.__module__.replace('.', os.sep))
        parent_subpkg = parent_path.replace(os.sep, '.')
        napoleon_type_aliases[name] = f":class:`~{parent_subpkg}.{name}`"

# Set autodoc aliases here to avoid recompiling everything every time
autodoc_type_aliases = {
    k: napoleon_type_aliases[k] for k in sorted(napoleon_type_aliases)
}


# -----------------------------------------------------------------------------
# Source code links (from numpy)
# -----------------------------------------------------------------------------

import inspect
from os.path import dirname, relpath


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
    if module is not None and not module.__name__.startswith('biocrnpyler'):
        # print("module is not None but doesnt start with biocrnpyler")
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except Exception:
        lineno = None

    fn = relpath(fn, start=dirname(biocrnpyler.__file__))

    if lineno:
        linespec = '#L%d-L%d' % (lineno, lineno + len(source) - 1)
    else:
        linespec = ''

    base_url = "https://github.com/BuildACell/BioCRNPyler/blob/"
    if release != version:  # development release
        # TODO: replace 'refactor-modules' with 'master' -> replaced with main
        # print("  --> ", base_url + "refactor-modules/control/%s%s" % (fn, linespec))
        return base_url + 'main/biocrnpyler/%s%s' % (fn, linespec)
    else:  # specific version
        return base_url + '%s/biocrnpyler/%s%s' % (version, fn, linespec)


# -- Options for doctest ----------------------------------------------

# Import biocrnpyler as bcp
doctest_global_setup = """
import numpy as np
import biocrnpyler as bcp
"""
