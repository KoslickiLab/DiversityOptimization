#! /usr/bin/env python3
import argparse
import os
import sys
import subprocess
import tempfile
import scipy.io as sio
import numpy as np
from time import time
from scipy.optimize import nnls
from scipy.sparse import vstack
sys.path.append("..")
from src.sparse_nnls import sparse_nnls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reconstructs population proportions with sparse_nnls and nnls and compares",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_file', type=str, help="File name of input database", required=True)
    parser.add_argument('-s', '--small_k', type=int,
                        help="small k-mer size to use (Note: values >14 will probably take too long)", default=6)
    parser.add_argument('-l', '--large_k', type=int,
                        help="large k-mer size to use (Note: values >14 will probably take too long)", default=12)
    parser.add_argument('-c', '--const', type=int, help="lambda (AKA const) value", default=10000)
    parser.add_argument('-q', '--q_value', type=float, help="q value", default=0.1)
    parser.add_argument('--count_complements', action="store_true",
                        help="count compliment of sequences as well", default=False)
    parser.add_argument('-r', '--reference', type=str, help="File name of reference database", required=True)

    args = parser.parse_args()
    input_file = args.input_file
    small_k = args.small_k
    large_k = args.large_k
    const = args.const
    q = args.q_value
    count_rev = args.count_complements
    reference = args.reference

    ## Check if input files exist:

    # input
    if not os.path.exists(input_file):
        raise Exception(f"The input file {input_file} does not appear to exist")
    # reference
    if not os.path.exists(reference):
        raise Exception(f"The reference reference file {reference} does not appear to exist")

    ## Create sensing matrices if they have not been created, and load them

    # large_k
    if not os.path.exists(f"{reference}_A_{large_k}.mat"):
        to_run = f"python ../src/./Form16SSensingMatrix.py -k {large_k} -i {reference} -o {reference}_A_{large_k}.mat {count_rev*'-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large sensing matrix")

    A_k_large = sio.loadmat(f"{reference}_A_{large_k}.mat")['A_k']

    # small_k
    if not os.path.exists(f"{reference}_A_{small_k}.mat"):
        to_run = f"python ../src/./Form16SSensingMatrix.py -k {small_k} -i {reference} -o {reference}_A_{small_k}.mat {count_rev*'-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small sensing matrix")

    A_k_small = sio.loadmat(f"{reference}_A_{small_k}.mat")['A_k']

    ## Create y vectors:

    with tempfile.NamedTemporaryFile() as temp_file:
        t_file = temp_file.name

        # large_k
        to_run = f"python ../src/./Form16SyVector.py -k {large_k} -i {input_file} -o {t_file} {count_rev*'-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large y vector")

        y_large = sio.loadmat(t_file)['y'].T

        # small_k
        to_run = f"python ../src/./Form16SyVector.py -k {small_k} -i {input_file} -o {t_file} {count_rev * '-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small y vector")

        y_small = sio.loadmat(t_file)['y'].T

    B = A_k_large > 0

    epsilon = 0.0001
    denom = np.power(B.T @ y_large, 1 - q) + epsilon
    f = 1/denom

    # nnls vs sparse_nnls

    start = time()
    x_nnls = nnls(vstack((f.T, const * A_k_small)).toarray(), np.append(0, const * y_small))[0]
    end = time()
    print(f"Time for nnls is {end - start}.")
    start = time()
    x_snnls = sparse_nnls(vstack((f.T, const * A_k_small)).tocsr(), np.append(0, const * y_small))
    end = time()
    print(f"Time for sparse_nnls is {end-start}.")
    print(f"l1 difference is {np.linalg.norm(x_nnls - x_snnls, 1)}.")
    print(f"l2 difference is {np.linalg.norm(x_nnls - x_snnls, 2)}.")
