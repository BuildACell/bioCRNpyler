# make file to automate certain tasks in BioCRNpyler installation and tests
#
# !!!make sure you use TABs instead of 4 spaces once you edit this file!!!
#

# test the core functionality of the toolbox
test :
	python setup.py test

# run the test suite with all dependencies
test_all :
	python setup.py .install[all]
	python -m pytest

# test for default mutable arguments in the code
flake8-mutable:
	flake8  --select M biocrnpyler
	if [[ `flake8  --select M biocrnpyler` ]]; then >&2 echo "default mutable argument detected"; exit 1 ; fi

install :
	python setup.py install

# running the test in Travis requires these packages
get_test_deps :
	pip install codecov
	pip install flake8-mutable
	pip install flake8