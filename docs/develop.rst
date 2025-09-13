.. currentmodule:: biocrnpyler

***************
Developer Notes
***************

This chapter contains notes for developers who wish to contribute to
the BioCRNpyler package.

Package Structure
=================

The bioCRNpyler package is maintained on GitHub, with documentation
hosted by ReadTheDocs and a mailing list on SourceForge:

  * Project home page: 
  * Source code repository: https://github.com/buildacell/BioCRNPyler
  * Documentation: https://biocrnpyler.readthedocs.io/
  * Issue tracker: https://github.com/buildacell/BioCRNPyler/issues

GitHub repository file and directory layout:
  - **biocrnpyler/** - main repository

    * LICENSE, Manifest, pyproject.toml, README.rst - package information

    * **biocrnpyler/** - primary package source code

      + **core/** -

        - __init__.py

    * **docs/** - user guide and reference manual

      + index.rst - main documentation index

      + conf.py, Makefile - sphinx configuration files

      + intro.rst

      + functions.rst, classes.rst, develop.rst -
        Reference Manual

      + **examples/**

        - \*.py, \*.rst - Python scripts (linked to ../examples/\*.py)

        - \*.ipynb - Jupyter notebooks (linked to ../examples.ipynb)

    * **examples/**

      + \*.py - Python scripts

      + \*.ipynb - Jupyter notebooks


Naming Conventions
==================

Generally speaking, standard Python and NumPy naming conventions are
used throughout the package.

* Python PEP 8 (code style): https://peps.python.org/pep-0008/


Filenames
---------

* Source files are lower case, usually less than 10 characters (and 8
  or less is better).

* Unit tests (in `Tests/*/`) are of the form `test_functinality.py`.


Class names
-----------

* Most class names are in camel case, with long form descriptions of
  the object purpose/contents (`EnzymaticReaction`).


Function names
--------------

* Function names are lower case with words separated by underscores.

* Function names usually describe what they do
  (`create_statefbk_iosystem`, `find_operating_points`) or what they
  generate (`input_output_response`, `find_operating_point`).


Parameter names
---------------

Function parameter names are not (yet) very uniform across the package.  A few
general patterns are emerging:

* Use longer description parameter names that describe the action or
  role (e.g., `trajectory_constraints` and `print_summary` in
  `optimal.solve_optimal_trajectory`.


Documentation Guidelines
========================

The bioCRNpyler package is documented using docstrings and Sphinx.
Reference documentation (class and function descriptions, with details
on parameters) should all go in docstrings.  User documentation in
more narrative form should be in the `.rst` files in `docs/`, where it
can be incorporated into the User Guide.  All significant
functionality should have a narrative description in the User Guide in
addition to docstrings.

Generally speaking, standard Python and NumPy documentation
conventions are used throughout the package:

* Python PEP 257 (docstrings): https://peps.python.org/pep-0257/
* Numpydoc Style guide: https://numpydoc.readthedocs.io/en/latest/format.html


General docstring info
----------------------

The guiding principle used to guide how docstrings are written is
similar to NumPy (as articulated in the `numpydoc style guide
<https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_):

   A guiding principle is that human readers of the text are given
   precedence over contorting docstrings so our tools produce nice
   output. Rather than sacrificing the readability of the docstrings,
   we have written pre-processors to assist Sphinx in its task.

To that end, docstrings in `python-control` should use the following
guidelines:

* Use single backticks around all Python objects. The Sphinx
  configuration file (`doc/conf.py`) defines `default_role` to be
  `py:obj`, so everything in a single backtick will be rendered in
  code form and linked to the appropriate documentation if it exists.

  - Note: consistent with numpydoc recommendations, parameters names
    for functions should be in single backticks, even though they
    don't generate a link (but the font will still be OK).

  - The `doc/_static/custom.css` file defines the style for Python
    objects and is configured so that linked objects will appear in a
    bolder type, so that it is easier to see what things you can click
    on to get more information.

  - By default, the string \`sys\` in docstrings would normally
    generate a link to the :mod:`sys` Python module.  To avoid this,
    `conf.py` includes code that converts \`sys\` in docstrings to
    \:code\:\`sys`, which renders as :code:`sys` (code style, with no
    link).  In ``.rst`` files this construction should be done
    manually, since ``.rst`` files are not pre-processed as a
    docstring.

* Use double backticks for inline code, such as a Python code fragments.

  - In principle single backticks might actually work OK given the way
    that the `py:obj` processing works in Sphinx, but the inclusion of
    code is somewhat rare and the extra two backticks seem like a
    small sacrifice (and far from a "contortion").

* Avoid the use of backticks and \:math\: for simple formulas where
  the additional annotation or formatting does not add anything.  For
  example "-c <= x <= c" (without the double quotes) in
  `relay_hysteresis_nonlinearity`.

  - Some of these formulas might be interpreted as Python code
    fragments, but they only need to be in double quotes if that makes
    the documentation easier to understand.

  - Examples:

      * \`dt\` > 0 not \`\`dt > 0\`\` (`dt` is a parameter)
      * \`squeeze\` = True not \`\`squeeze = True\`\` nor squeeze = True.
      * -c <= x <= c not \`\`-c <= x <= c\`\` nor \:math\:\`-c \\leq x
        \\leq c`.
      * \:math\:\`|x| < \\epsilon\` (becomes :math:`|x| < \epsilon`)

* Built-in Python objects (True, False, None) should be written with no
  backticks and should be properly capitalized.

  - Another possibility here is to use a single backtick around
    built-in objects, and the `py:obj` processing will then generate a
    link back to the primary Python documentation.  That seems
    distracting for built-ins like `True`, `False` and `None` (written
    here in single backticks) and using double backticks looks fine in
    Sphinx (``True``, ``False``, ``None``), but seemed to cross the
    "contortions" threshold.

* Strings used as arguments to parameters should be in single
  (forward) ticks ('eval', 'rows', etc) and don't need to be rendered
  as code if just listed as part of a docstring.

  - The rationale here is similar to built-ins: adding 4 backticks
    just to get them in a code font seems unnecessary.

  - Note that if a string is included in Python assignment statement
    (e.g., ``method='slycot'``) it looks quite ugly in text form to
    have it enclosed in double backticks (\`\`method='slycot'\`\`), so
    OK to use method='slycot' (no backticks) or `method` = 'slycot'
    (backticks with extra spaces).


Function docstrings
-------------------

Follow numpydoc format with the following additional details:

* All functions should have a short (< 64 character) summary line that
  starts with a capital letter and ends with a period.

* All parameter descriptions should start with a capital letter and
  end with a period.  An exception is parameters that have a list of
  possible values, in which case a phrase sending in a colon (:)
  followed by a list (without punctuation) is OK.

* All parameters and keywords must be documented.  The
  `docstrings_test.py` unit test tries to flag as many of these as
  possible.

* Include an "Examples" section for all non-trivial functions, in a
  form that can be checked by running `make doctest` in the `doc`
  directory.  This is also part of the CI checks.


Class docstrings
----------------

Follow numpydoc format with the follow additional details:

* Parameters used in creating an object go in the class docstring and
  not in the `__init__` docstring (which is not included in the
  Sphinx-based documentation).  OK for the `__init__` function to have
  no docstring.

* Parameters that are also attributes only need to be documented once
  (in the "Parameters" or "Additional Parameters" section of the class
  docstring).

* Attributes that are created within a class and that might be of
  interest to the user should be documented in the "Attributes"
  section of the class docstring.

* Classes should not include a "Returns" section (since they always
  return an instance of the class).

* Functions and attributes that are not intended to be accessed by
  users should start with an underscore.


User Guide
----------

The purpose of the User Guide is provide a *narrative* description of
the key functions of the package.  It is not expected to cover every
command, but should allow someone who knows about biological circuit
design to get up and running quickly.

The User Guide consists of chapters that are each their own separate
`.rst` file and each of them generates a separate page.  Chapters are
divided into sections whose names appear in the index on the left of
the web page when that chapter is being viewed.  In some cases a
section may be in its own file, included in the chapter page by using
the `include` directive.

Sphinx files guidelines:

* Each file should declare the `currentmodule` at or near the top of
  the file.

* When possible, sample code in the User Guide should use Sphinx
  doctest directives so that the code is executed by `make doctest`.
  Two styles are possible: doctest-style blocks (showing code with a
  prompt and the expected response) and code blocks (using the
  `testcode` directive).

* Unlike docstrings, the documentation in the User Guide should use
  backticks and \:math\: more liberally when it is appropriate to
  highlight/format code properly.  However, Python built-ins should
  still just be written as True, False, and None (no backticks), for
  formatting consistency.

  - The Sphinx documentation is not read in "raw" form, so OK to add
    the additional annotations.

  - The Python built-ins occur frequently and are capitalized, and so
    the additional formatting doesn't add much and would be
    inconsistent if you jump from the User Guide to the Reference
    Manual (e.g., to look at a function more closely via a link in the
    User Guide).


Reference Manual
----------------

The Reference Manual should provide a fairly comprehensive description
of every class and function.  All primary functions and classes bust
be included here, since the Reference Manual generates the stub files
used by Sphinx.
