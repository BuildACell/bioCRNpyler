from setuptools import setup

# Get the long description from the README file
with open('README.md') as fp:
    long_description = fp.read()

setup(
    name='biocrnpyler',
    version='0.1',
    author='William Poole',
    author_email='wpoole@caltech.edu',
    url='https://github.com/WilliamIX/biocrnpyler/',
    description='A chemical reaction network compiler for generating large biological circuit models',
    long_description=long_description,
    packages=['biocrnpyler'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
        'Operating System :: POSIX',
        'Operating System :: Unix'
        'Operating System :: MacOS'
    ],
    install_requires=[
          'python-libsbml',
          'numpy',
          'nose',
          'matplotlib',
          'networkx',
          'bokeh',
          'ForceAtlas2',
      ],
    test_suite='nose.collector',
    tests_require=['nose'],
)
