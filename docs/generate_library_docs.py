import ast
from pathlib import Path


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
        lines = [f'- :class:`biocrnpyler.mechanisms.{module}.{cls}`' for cls in classes]
        sections.append((section_title, lines))
    write_rst_file(out_file, 'Mechanisms', sections)


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
        lines = [f'- :class:`biocrnpyler.mixtures.{module}.{cls}`' for cls in classes]
        sections.append((section_title, lines))
    write_rst_file(out_file, 'Mixtures', sections)


def generate_components_rst(src_root, out_file):
    comp_dir = Path(src_root) / 'biocrnpyler' / 'components'
    sections = []
    # top-level components
    top_py = [p for p in comp_dir.glob('*.py') if p.name != '__init__.py']
    if top_py:
        lines = []
        for py in sorted(top_py):
            module = py.stem
            classes = get_classes_in_file(py)
            for cls in classes:
                lines.append(f'- :class:`biocrnpyler.components.{module}.{cls}`')
        sections.append(('Components', lines))
    # DNA subpackage
    dna_dir = comp_dir / 'dna'
    if dna_dir.exists():
        for py in sorted(dna_dir.glob('*.py')):
            if py.name == '__init__.py':
                continue
            module = py.stem
            classes = get_classes_in_file(py)
            if not classes:
                continue
            title = f'DNA: {module.replace("_", " ").title()}'
            lines = [f'- :class:`biocrnpyler.components.dna.{module}.{cls}`' for cls in classes]
            sections.append((title, lines))
    write_rst_file(out_file, 'Components', sections)


if __name__ == '__main__':
    # 1) The directory *this script* lives in → that's your docs/ folder
    docs_dir = Path(__file__).resolve().parent

    # 2) The project root is the parent of docs/
    project_root = docs_dir.parent

    # 3) Source code root is that project root
    src_root = project_root

    # 4) Ensure docs/ exists (it does, since the script is already there, but no harm)
    docs_dir.mkdir(exist_ok=True)

    # 5) Generate into docs/
    generate_components_rst(src_root, docs_dir / '_autogen_components.rst')
    generate_mechanisms_rst(src_root, docs_dir / '_autogen_mechanisms.rst')
    generate_mixtures_rst(src_root, docs_dir / '_autogen_mixtures.rst')

    print('RST auto-generation complete! Files written to:', docs_dir)