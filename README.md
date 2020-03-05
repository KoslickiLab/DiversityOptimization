# DiversityOptimization
Minimizing biological diversity to improve microbial taxonomic reconstruction

## Intent
The main point of the repository is to implement in python the theoretical results contained in [this manuscript](https://www.biorxiv.org/content/10.1101/2020.01.23.916924v1).

Work is ongoing, so the code may not work as intended until an initial release is made

# Repository stucture
A basic layout of the folder structure is as follows:
```
README.md  # this file
MatlabCode/  # basic matlab code containing a simplified version of https://github.com/dkoslicki/MinimizeBiologicalDiversity
          README.md  # a file describing what is in this folder
          Data/  # data required to run the matlab code
PythonCode/
          README.md  # a file describing what is in this folder
          requirements.txt  # a file containing the necessary requirements to run the python code
          src/  # a folder containing the modules that contain the port of the matlab code
          scripts/  # a folder containing wrappers that call the src/ modules in a command line interface friendly way
          tests/  # a collection of test for the various python components
          data/  # any required data for the python code
```
