# Description
This will contain the core python algorithms and functions. For an example of what this should look like. See [this repo](https://github.com/dkoslicki/PressPurt/tree/master/Python_version/src). In particular, [this file](https://github.com/dkoslicki/PressPurt/blob/master/Python_version/src/NumSwitch.py) shows what such a module might look like structurally. 

# MinDivLP
A basic, regularized version of the MinDivLP algorithm. Reconstructs a vector x given information in the form of y=Ax.

# Form16SSensingMatrix.py
This will form the sensing matrix `A` when given a database of 16S FASTA formatted bacterial genomes.
Note: this will depend on the [dna-utils](https://github.com/dkoslicki/dna-utils) tool. Installation instructions to follow,
but basically
```bash
make 
sudo make install  # may not need sudo if you want only a local install
```