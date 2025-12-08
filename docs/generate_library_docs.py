# generate_library_docs.py - create docs for components, mechanisms, mixtures
#
# This script creates files for the component, mechanism, and mixture
# libraries within BioCRNpyler, as a meanso of making sure that all such
# elements are document.  The files are of the form _autogen_<type>.rst and
# can be included in the main documentation to generate sections for each
# module in the package, with a list of the relevant classes.
#
# For each module, the module description is pulled from the module
# docstring (at the top of the file for that module), followed by a
# autosummary table of all classes defined in the file.

import ast
from pathlib import Path

# Text to generate a summary table of the methods available in the module
autosummary = """
.. autosummary::
   :toctree: generated/
   :nosignatures:
   :recursive:
"""

# Text to extract the module description at the start of the session
automodule = """
.. automodule:: biocrnpyler.{0}.{1}
   :no-index:

.. currentmodule:: biocrnpyler.{2}
"""


def get_classes_in_file(py_path):
    """
    Parse a Python file and return all top-level class names.
    """
    with open(py_path, 'r', encoding='utf-8') as f:
        tree = ast.parse(f.read(), filename=str(py_path))
    return [node.name for node in tree.body if isinstance(node, ast.ClassDef)]


def write_rst_file(out_path, header, sections):
    """
    Write an RST file with a main header and multiple sections.
    sections: list of tuples (section_title, list of lines)
    """
    with open(out_path, 'w', encoding='utf-8') as f:
        # Main title
        for title, lines in sections:
            f.write(title + '\n')
            f.write('-' * len(title) + '\n')
            for line in lines:
                f.write(line + '\n')
            f.write('\n')


# Generate sections for all mechanisms
def generate_mechanisms_rst(src_root, out_file):
    mech_dir = Path(src_root) / 'biocrnpyler' / 'mechanisms'
    sections = []
    for py in sorted(mech_dir.glob('*.py')):
        if py.name == '__init__.py':
            continue
        module = py.stem
        classes = get_classes_in_file(py)
        if not classes:
            continue
        section_title = module.replace('_', ' ').title()
        lines = (
            [automodule.format('mechanisms', module, 'mechanisms')]
            + [autosummary]
            + [f"    {cls}" for cls in classes]
        )
        sections.append((section_title, lines))
    write_rst_file(out_file, 'Mechanisms', sections)


# Generate sections for all mixtures
def generate_mixtures_rst(src_root, out_file):
    mix_dir = Path(src_root) / 'biocrnpyler' / 'mixtures'
    sections = []
    for py in sorted(mix_dir.glob('*.py')):
        if py.name == '__init__.py':
            continue
        module = py.stem
        classes = get_classes_in_file(py)
        if not classes:
            continue
        section_title = module.replace('_', ' ').title()
        lines = (
            [automodule.format('mixtures', module, 'mixtures')]
            + [autosummary]
            + [f"    {cls}" for cls in classes]
        )
        sections.append((section_title, lines))
    write_rst_file(out_file, 'Mixtures', sections)


# Generate sections for all mixtures
def generate_components_rst(src_root, out_file):
    comp_dir = Path(src_root) / 'biocrnpyler' / 'components'
    sections = []

    # top-level components
    top_py = [p for p in comp_dir.glob('*.py') if p.name != '__init__.py']
    if top_py:
        for py in sorted(top_py):
            module = py.stem
            classes = get_classes_in_file(py)
            section_title = module.replace('_', ' ').title()
            lines = (
                [automodule.format('components', module, 'components')]
                + [autosummary]
                + [f"    {cls}" for cls in classes]
            )
            sections.append((section_title, lines))

    # DNA subpackage - include as a single subsection
    dna_dir = comp_dir / 'dna'
    if dna_dir.exists():
        title = "DNA Assemblies"
        lines = [automodule.format('components', 'dna', 'components.dna')] + [
            autosummary
        ]
        for py in sorted(dna_dir.glob('*.py')):
            if py.name == '__init__.py':
                continue
            module = py.stem
            classes = get_classes_in_file(py)
            if not classes:
                continue
            lines += [f'   {cls}' for cls in classes]
        sections.append((title, lines))
    write_rst_file(out_file, 'Components', sections)


if __name__ == '__main__':
    # 1) The directory *this script* lives in → that's your docs/ folder
    docs_dir = Path(__file__).resolve().parent

    # 2) The project root is the parent of docs/
    project_root = docs_dir.parent

    # 3) Source code root is that project root
    src_root = project_root

    # 4) Ensure docs/ exists (it does, since script is there, but no harm)
    docs_dir.mkdir(exist_ok=True)

    # 5) Generate into docs/
    generate_components_rst(src_root, docs_dir / '_autogen_components.rst')
    generate_mechanisms_rst(src_root, docs_dir / '_autogen_mechanisms.rst')
    generate_mixtures_rst(src_root, docs_dir / '_autogen_mixtures.rst')

    print('RST auto-generation complete! Files written to:', docs_dir)
