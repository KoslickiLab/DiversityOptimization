import scipy.io as sio
import numpy as np
import numpy.random as rand
import time
import sys
import os
sys.path.append(os.path.abspath("../.."))  # make sure python knows where to find the code
from PythonCode.src.MinDivLP import MinDivLP
from numpy.linalg import norm
from matplotlib import pyplot as plt

## import the data
# set variables
small_k = 4  # smaller k-mer size
large_k = 6  # larger k-mer size

A_k = sio.loadmat('../data/97_otus_subset.fasta_A_%d.mat' % large_k)
A_k_large = A_k['A_k']

A_k = sio.loadmat('../data/97_otus_subset.fasta_A_%d.mat' % small_k)
A_k_small = A_k['A_k']

del A_k

## sub-select the data so things run quickly
cols_vs_rows = 3  # fix 3-times more columns than rows
num_species = cols_vs_rows * 4 ** small_k  # reduce the number of columns of
# the sensing matrix so pictures will be generated in a reasonable amount
# of time.
A_k_small = A_k_small[:, 0:num_species]  # Note: these are already column-normalized to be 1
A_k_large = A_k_large[:, 0:num_species]

## Generate the B matrix
B = (A_k_large > 0)

## Set some variables
q = .1  # fixed, small q value s.t. 0<q<1
support_size = 15  # number of non-zero entries in the simulated ground truth

## Create the simulated ground truth

supp = rand.choice(range(0, num_species), size=(1, support_size), replace=False)  # location of the support
true_x = np.zeros((num_species, 1))  # the true x vector we are trying to reconstruct
true_x[supp] = rand.random((support_size, 1))  # populate with random data

true_x = true_x / sum(true_x)  # normalize to be a probability vector

# Noiseless y-vectors
y_small_true = A_k_small @ true_x
y_large_true = A_k_large @ true_x

# Noisy y-vectors
noise_eps = .00001  # size of noise to add
y_small_noise = (A_k_small @ true_x) + noise_eps * rand.random((A_k_small.shape[0], 1))
y_small_noise = y_small_noise / sum(y_small_noise)

## Noisy computations

const = 10000
start = time.time()
x_star = MinDivLP(A_k_small.toarray(), A_k_large.toarray(), y_small_noise, y_large_true, const, q).reshape(num_species, 1)
end = time.time()

## Measure the reconstruction accuracy
error_l1 = norm(x_star - true_x, 1)
error_l2 = norm(x_star - true_x, 2)
print('L1 error is: %f.' % error_l1)
print('L2 error is: %f.' % error_l2)
print('Time elapsed for MinDivLP: %f seconds.' % (end - start))

## Sanity check plot
plt.plot(range(1, num_species + 1), x_star, 'b.', label='Reconstructed')
plt.plot(range(1, num_species + 1), true_x, 'r.', label='True')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()
plt.show()
