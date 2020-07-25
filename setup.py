from setuptools import setup

# Get the long description from the README file
with open('README.md') as fp:
    long_description = fp.read()

setup(
    name='biocrnpyler',
    version='0.2.1',
    author='BuildACell',
    url='https://github.com/BuildACell/biocrnplyler/',
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
          'matplotlib<=3.2.2',
          'networkx',
          'bokeh>=1.4.0',
          'fa2',
      ],
    setup_requires=["pytest-runner"],
    tests_require=["pytest", "pytest-cov", "nbval"],
)
