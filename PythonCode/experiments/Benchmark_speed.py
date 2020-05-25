import scipy.io as sio
import numpy as np
import numpy.random as rand
import time
import sys
import os
import csv
from multiprocessing import Process
from statistics import mean, median, stdev

sys.path.append(os.path.abspath("../.."))  # make sure python knows where to find MinDivLP
from PythonCode.src.MinDivLP import MinDivLP

## Data writing information
resultCSV = "results.csv"
overwrite = True

## Variables for data acquisition
data_dir = "../data"  # input directory containing full sensing matrices
no_c = False  # append "_no_c" to exclude reverse complement k-mers

## Parameters of the variables being tested
k_sizes = [4, 6, 8, 10, 12, 13]  # k_mer size for A_k_large to be tested
col_counts = [99322, 50000, 10000, 1000, 100]  # numbers of included columns of A_k (in decreasing order)
support_sizes = [25, 50, 100, 500, 1000, 5000]  # support sizes of true_x
N = 20  # iterations at each combination of k_size, num_cols, support_size

## Other parameters
q = .1  # fixed, small q value s.t. 0<q<1
lamb = 10000  # lambda


## Define testing function to parallelize
def test(large_k, col_counts, support_sizes, N, q, lamb, csvname):
    ## import the data, set small_k
    A_k = sio.loadmat('{0}/97_otus.fasta_A_{1}{2}.mat'.format(data_dir, large_k, no_c * "_no_c"))
    A_k_large = A_k['A_k']

    small_k = 4  # smaller k-mer size
    A_k = sio.loadmat('{0}/97_otus.fasta_A_{1}{2}.mat'.format(data_dir, small_k, no_c * "_no_c"))
    A_k_small = A_k['A_k']

    del A_k

    for num_cols in col_counts:

        A_k_small = A_k_small[:, 0:num_cols]  # Note: these are already column-normalized to be 1
        A_k_large = A_k_large[:, 0:num_cols]

        for support_size in support_sizes:  # number of non-zero entries in the simulated ground truth

            support_size = min(support_size, A_k_small.shape[1])
            times = list()

            for _ in range(N):
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
                y_small_noise = y_small_true + noise_eps * rand.random((A_k_small.shape[0], 1))
                y_small_noise = y_small_noise / sum(y_small_noise)

                ## Noisy computations
                start_time = time.time()
                x_star = MinDivLP(A_k_small.toarray(), A_k_large, y_small_noise, y_large_true, lamb, q).reshape(
                    num_cols, 1)
                times.append(time.time() - start_time)

            ## Append csv file with results
            with open(csvname, "a", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(
                    [4 ** large_k, num_cols, support_size, N, min(times), max(times), median(times), mean(times),
                     stdev(times)])

            print("({0}, {1}, {2})".format(large_k, num_cols, support_size))


if __name__ == '__main__':

    if overwrite:
        ## Save first row of csv
        with open(resultCSV, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                ['rows', 'columns', 'support size', 'N', 'min', 'max', 'median', 'mean', 'std dev'])

    ## Parallelize for each value of large_k
    processes = list()
    for large_k in k_sizes:
        p = Process(target=test, args=(large_k, col_counts, support_sizes, N, q, lamb, resultCSV))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()
