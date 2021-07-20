import pytest
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from os import listdir, getcwd, pardir
from os.path import join, abspath


def run_notebook(filename, path):
    with open(filename) as f:
        ep = ExecutePreprocessor()
        nb = nbformat.read(filename, nbformat.NO_CONVERT)
        try:
            ep.preprocess(nb, {'metadata': {'path': path}})
        except CellExecutionError:
            msg = f"\nError executing the notebook {join(path, filename)}\n" 
            print(msg)
            raise

cwd = getcwd()
paths = [join(cwd, "examples"), join(cwd, "examples", "Specialized Component Tutorials")]
nb_names = []
for p in paths:
    nb_names += [join(p, f) for f in listdir(p) if f[-6:] == ".ipynb"]

@pytest.mark.parametrize("nb", nb_names)
def test_jupyter_notebooks(nb):
    path = abspath(join(nb, pardir))
    run_notebook(nb, path)
