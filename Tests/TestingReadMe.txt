In order to run tests:
1. cd to the biocrnpyler/tests folder
2. type: coverage run -m pytest

Building Tests:
Please test all new additions to BioCRNpyler. Tests can be added with existing tests (which are organized by object) or as new test sets for more complex features. 
All tests should be commented so that it is clear what the test is meant to check for.

Tests are divided into two folders:
Unit: these are standard unit test of different parts of the code
Combinatorial: these are more involved tests which compile Mixtures, Components, and Mechanisms together in combinatorial sets. These are much slower to run.