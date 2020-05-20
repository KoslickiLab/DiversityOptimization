import scipy.io as sio
import numpy as np
import numpy.random as rand
import time
import sys
import os
import csv
from numpy.linalg import norm
from statistics import mean, median, stdev

sys.path.append(os.path.abspath("../.."))  # make sure python knows where to find MinDivLP
from PythonCode.src.MinDivLP import MinDivLP

## Variables for data acquisition
data_dir = ""  # input directory containing full sensing matrices
no_c = False  # append "_no_c" to exclude reverse complement k-mers

## Parameters of the variables being tested
k_sizes = [4, 6, 8, 10, 12]  # k_mer size for A_k_large
col_counts = [10000, 1000, 100]  # numbers of included columns of A_k (in decreasing order)
support_sizes = [25, 50, 75]  # support sizes of true_x
N = 20  # iterations at each combination of k_size, num_cols, support_size

## Other parameters
q = .1  # fixed, small q value s.t. 0<q<1
const = 10000  # lambda

results = list()

for large_k in k_sizes:  # num_rows is 4^k_size

    ## import the data, set small_k
    small_k = 4  # smaller k-mer size

    A_k = sio.loadmat('{0}/97_otus_subset.fasta_A_{1}{2}.mat'.format(data_dir, large_k, no_c * "_no_c"))
    A_k_large = A_k['A_k']

    A_k = sio.loadmat('{0}/97_otus_subset.fasta_A_{1}{2}.mat'.format(data_dir, small_k, no_c * "_no_c"))
    A_k_small = A_k['A_k']

    del A_k

    for num_cols in col_counts:

        A_k_small = A_k_small[:, 0:num_cols]  # Note: these are already column-normalized to be 1
        A_k_large = A_k_large[:, 0:num_cols]

        for support_size in support_sizes:  # number of non-zero entries in the simulated ground truth
            times = list()
            l1errors = list()
            l2errors = list()

            for i in range(0, N):
                ## Create the simulated ground truth
                supp = rand.choice(range(0, support_size), size=(1, support_size), replace=False)  # location of the
                # support
                true_x = np.zeros((num_cols, 1))  # the true x vector we are trying to reconstruct
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
                start_time = time.time()
                x_star = MinDivLP(A_k_small.toarray(), A_k_large, y_small_noise, y_large_true, const, q).reshape(
                    num_cols, 1)
                times.append(time.time() - start_time)

                ## Measure the reconstruction accuracy
                l1errors.append(norm(x_star - true_x, 1))
                l2errors.append(norm(x_star - true_x, 2))

            info = np.array([4 ** large_k, num_cols, support_size, min(times), max(times), median(times), mean(times),
                             stdev(times), mean(l1errors), mean(l2errors)])

            results.append(info)

    print('finished large_k = {}'.format(large_k))  # Print to ensure running

## Save summary statistics to a csv file for analysis
with open("results.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(
        ['rows', 'columns', 'support size', 'min', 'max', 'median', 'mean', 'std dev', 'mean l1 error', 'mean l2 '
                                                                                                        'error'])
    writer.writerows(results)
