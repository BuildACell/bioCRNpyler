# make file to automate certain tasks in BioCRNpyler installation and tests
#
# !!!make sure you use TABs instead of 4 spaces once you edit this file!!!
#
.PHONY: docs
# test the core functionality of the toolbox
test :
	#pip install .[test]
	python setup.py test
	#pytest Tests/Unit
# run the test suite with all dependencies
test_all :
	pip install .[all]
	python setup.py test
	#python setup.py install
	#pytest Tests

# test for default mutable arguments in the code
flake8-mutable :
	flake8  --select M biocrnpyler
	if [[ `flake8  --select M biocrnpyler` ]]; then >&2 echo "default mutable argument detected"; exit 1 ; fi

install :
	python setup.py install

# running the test in Travis requires these packages
get_test_deps :
	pip install codecov
	pip install flake8-mutable
	pip install flake8

docs :
	# this has to be one line to run everything in the docs folder
	cd docs; echo "Running Sphinx docs generator"; sphinx-apidoc -o source/ ../biocrnpyler; python generate_nblinks.py; make clean && make html;

TAGS: biocrnpyler/*.py biocrnpyler/*/*.py biocrnpyler/*/*/*.py docs/*.rst \
    docs/examples/*.ipynb
	ftags $^
