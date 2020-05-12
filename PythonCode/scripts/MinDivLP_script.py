import scipy.io as sio
import numpy as np
import random as rand
from numpy.linalg import norm
from MinDivLP import MinDivLP

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
num_species = cols_vs_rows*4**small_k  # reduce the number of columns of
# the sensing matrix so pictures will be generated in a reasonable amount
# of time.
A_k_small = A_k_small[:, 0:num_species]  # Note: these are already column-normalized to be 1
A_k_large = A_k_large[:, 0:num_species]

## Generate the B matrix
B = (A_k_large > 0)

## Set some variables
q = .1  # fixed, small q value s.t. 0<q<1
support_size = 15  # number of non-zero entries in the simulated ground truth

## create the simulated ground truth

supp = rand.sample(range(0, num_species), support_size)  # location of the support
true_x = np.zeros((num_species, 1))  # the true x vector we are trying to reconstruct

for i in supp:
    true_x[i] = rand.uniform(0, 1)  # populate with random data

true_x = true_x / sum(true_x)  # normalize to be a probability vector

# noisless y-vectors
y_small_true = A_k_small @ true_x
y_large_true = A_k_large @ true_x

# noisy y-vectors
slop_factor = 2  # this is for the non-regularized MinDivLP on noisy data
noise_eps = .00001  # size of noise to add
y_small_noise = A_k_small @ true_x

for i in range(0, A_k_small.shape[0]):
    y_small_noise[i] += noise_eps*rand.uniform(0, 1)  # add only noise to the small y-vector
# why was abs() included? # need to check multiplication of csc_matrix with ndarray
y_small_noise = y_small_noise / sum(y_small_noise)

## Noisy computations

const = 10000

# temporary
from scipy.optimize import nnls

B = A_k_large > 0
f = 1 / np.power(B.T @ y_large_true, 1 - q)

x_star = nnls(np.vstack((f.T, const * A_k_small)), np.append(0, const * y_small_noise))
x_star = x_star / sum(x_star)

x_star = MinDivLP(A_k_small, A_k_large, y_small_noise, y_large_true, const, q)

## Measure the reconstruction accuracy
error_l1 = norm(x_star - true_x, 1)
error_l2 = norm(x_star - true_x, 2)
print('L1 error is: %f\n' % error_l1)
print('L2 error is: %f\n' % error_l2)