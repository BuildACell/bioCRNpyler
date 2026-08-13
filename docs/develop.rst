.. currentmodule:: biocrnpyler

***************
Developer Notes
***************

This chapter contains notes for developers who wish to contribute to
the BioCRNpyler package.

Package Structure
=================

The BioCRNpyler package is maintained on GitHub, with documentation
hosted by ReadTheDocs:

  * Documentation: https://biocrnpyler.readthedocs.io/
  * Source code repository: https://github.com/buildacell/BioCRNpyler
  * Issue tracker: https://github.com/buildacell/BioCRNpyler/issues

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


Naming Conventions and Code Style
=================================

Generally speaking, standard Python and NumPy naming conventions are
used throughout the package.

* Python PEP 8 (code style): https://peps.python.org/pep-0008/

The `ruff <https://docs.astral.sh/ruff>`_ linter and code formatter
should be used to maintain a consistent code style, using the
following commands (from the main package directory):

* `ruff format`: apply standard code style
* `ruff check`: check against standard coding and formating conventions

The `pyproject.toml` file contains the appropriate
directives to implement the desired style, which follows the `Black
<https://pypi.org/project/black>`_ formatting conventions.


Filenames
---------

* Source files are lower case, usually less than 10 characters (and 8
  or less is better).

* Unit tests (in `Tests/*/`) are of the form `test_functionality.py`.


Class names
-----------

* Most class names are in camel case, with long form descriptions of
  the object purpose/contents (`~core.ChemicalReactionNetwork`).


Function and method names
-------------------------

* Function and method names are lower case with words separated by underscores.

* Function and method names usually describe what they do
  (`~core.Mixture.compile_crn`,
  `~core.ChemicalReactionNetwork.simulate_with_bioscrape_via_sbml`).


Parameter names
---------------

Function parameter names are not (yet) very uniform across the package.  A
few general patterns are emerging:

* Use longer description parameter names that describe the action or role
  (e.g., `overwrite_parameters` and `parameter_file` in `~core.Mixture`.


Documentation Guidelines
========================

The BioCRNpyler package is documented using docstrings and Sphinx.
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

To that end, docstrings in `biocrnpyler` should use the following
guidelines:

* Use single backticks around all Python objects. The Sphinx
  configuration file (`doc/conf.py`) defines `default_role` to be
  `py:obj`, so everything in a single backtick will be rendered in
  code form and linked to the appropriate documentation if it exists.

  - Note: consistent with numpydoc recommendations, parameters names
    for functions should be in single backticks, even though they
    do not generate a link (but the font will still be OK).

  - The `doc/_static/custom.css` file defines the style for Python
    objects and is configured so that linked objects will appear in a
    bolder type, so that it is easier to see what things you can click
    on to get more information.

* Use double backticks for inline code, such as a Python code fragments.

  - In principle single backticks might actually work OK given the way
    that the `py:obj` processing works in Sphinx, but the inclusion of
    code is somewhat rare and the extra two backticks seem like a
    small sacrifice (and far from a "contortion").

* Mathematical equations in docstrings can be written using LaTeX
  syntax.  Inline equations should be enclosed in a single set of
  dollar signs ($), with no space before or afer the equation (so
  `$A + B --> C$` to produce :math:`A + B \rightarrow C`.  Displayed
  equations can be written in one of two ways::

    text

        $$ A + B --> C $$

    text

  or

  .. code::

     text
     $$
         A + B --> C
     $$
     text

  both of which will produce:

  .. math::
    A + B \rightarrow C

  The first form is prefered when the contents fit on a single line.
  The second form can be used for equations that require more than one
  line in the docstring.  Multi-line formulas are also possible, using
  standard `aligned` syntax from LaTex.

* Several string substitutions are implemented on equation to allow
  docstrings that are easier to read in plain text while being
  formatted in LaTex for sphinx-processed documentation:

  - --> becomes :math:`\rightarrow`
  - <--> becomes :math:`\rightleftharpoons`
  - ... becomes :math:`\dots`
  - >> becomes :math:`\gg`
  - << becomes :math:`\ll`
  - A:B becomes :math:`A \mathord{:} B` (no space)

  In addition, any strings contained in single quotes will be
  converted to text::

    $$ 'substrate' + 'enzyme' --> 'substrate':'enzyme' $$

  becomes:

  .. math::

     \text{substrate} + \text{enzyme} \rightarrow
     \text{substrate} \mathord{:} \text{enzyme}

  Note: docstrings are processed line by line, so any inline equations
  (single $) must appear in the same line.

* As a general rule, chemical species in reactions should be in regular
  (non-math) font.  For chemical species that appear in regular text, no
  annotations are required.  For chemical species in equations, use single
  quotes to allow the preprocessing to convert to text.  When refering to
  the concentration of a chemical species, use square brackets ([A]; the
  chemical species will automatically be converted to regular text).

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
  (forward) ticks ('eval', 'rows', etc) and do not need to be rendered
  as code if just listed as part of a docstring.

  - The rationale here is similar to built-ins: adding 4 backticks
    just to get them in a code font seems unnecessary.

  - Note that if a string is included in Python assignment statement (e.g.,
    ``material_type='protein'``) it looks quite ugly in text form to have
    it enclosed in double backticks (\`\`material_type='protein'\`\`), so
    OK to use material_type='protein' (no backticks) or `material_type` =
    'protein' (backticks with extra spaces).


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
  possible. [not yet implemented]

* Include an "Examples" section for all non-trivial functions, in a
  form that can be checked by running `make doctest` in the `doc`
  directory.  This is also part of the CI checks. [not yet implemented]


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

Contributing and Releases
=========================

AI policy
---------

Follow the `AI Policy of NumPy`_ when writing issues and pull requests,
with one exception: AI coding tools may be used to draft the text of
issues and pull requests, not just to translate or edit it.  The
contributor is responsible for reviewing and editing any drafted text
before submitting it, and must disclose the AI assistance as described in
that policy.

.. _AI Policy of NumPy: https://numpy.org/doc/stable/dev/ai_policy.html

Releasing new versions
----------------------

BioCRNpyler uses semantic versioning (`MAJOR.MINOR.PATCH`) and Git tags for
automatic versioning of package distributions. The package wheels/sdists
are published to PyPI via GitHub Actions using the `pypi-release.yml`
action that can be triggered manually (with appropriate
permissions). Follow this checklist when planning a new release.

Release checklist for a new stable version:

1. Ensure your working tree is clean and CI is green.

2. Update the docs as needed.

3. Commit and push all changes.

4. Create an annotated tag on the release commit::

     git tag -a vX.Y.Z -m "Release X.Y.Z"
     git push --tags

5. Start the `pypi-release` workflow in GitHub Actions:

   * Prefer running the workflow from the tag, or
   * Provide the tag name in the `ref` input on the workflow.

6. Once the job has finished, you can verify from a fresh environment::

     python -m pip install -U biocrnpyler
     python -c "import biocrnpyler; print(biocrnpyler.__version__)"

7. Create a GitHub Release for the tag and paste the CHANGELOG entry.


Reference Manual
----------------

The Reference Manual should provide a fairly comprehensive description
of every class and function. All primary functions and classes must
be included here, since the Reference Manual generates the stub files
used by Sphinx.
