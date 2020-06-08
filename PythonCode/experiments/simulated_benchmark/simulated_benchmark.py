import sys
import os
import tempfile
import subprocess
import csv
import numpy as np
import scipy.io as sio
import argparse
import traceback

sys.path.append(os.path.abspath("../../.."))
from PythonCode.experiments.simulated_benchmark.make_true_x import make_true_x
from PythonCode.src.MinDivLP import MinDivLP
from scipy.linalg import norm

## Data writing information
resultCSV = "big_results.csv"

## Parameters for parallelization
large_k_vals = [8, 10, 12]  # larger k-mer size
noisy = False

# Location of bbmap
#bbmap_dir = "bbmap"  # Chase's system
#bbmap_dir = "/home/dkoslicki/Documents/bbmap"  # David's system


def calculateDiv(x, A, q):
    y = A @ x
    B = A > 0
    f = 1 / np.power(B.T @ y, 1 - q)

    return np.power(x @ f, 1 / (1 - q))


def simulate(large_k, mock_maker, csv_name, bbmap_dir):
    i = 0

    ## Iterations at each set of parameters
    N = 1

    # Optimization Constants
    q_vals = [0.1, 0.01, 0.001]
    lambda_vals = [1000, 10000, 20000]

    # k values
    small_k_vals = [4, 6, 8]  # smaller k-mer size

    # Simulation parameters

    coverage_vals = [20, 40, 60]
    supportsize_vals = [1, 2, 5, 15, 25]
    full_reference_genome = "../../data/97_otus_subset.fasta"
    reference_genome = f"./mock_reference.fa"

    A_k_large_full = sio.loadmat(f"../../data/97_otus.fasta_A_{large_k}.mat")['A_k'][:, 0:10000]

    for small_k in small_k_vals:
        # faster to load than create
        A_k_small_full = sio.loadmat(f"../../data/97_otus.fasta_A_{small_k}.mat")['A_k'][:, 0:10000]

        for support_size in supportsize_vals:

            for coverage in coverage_vals:

                for _ in range(N):

                    ## Create ground truth mock metagenome and mock reference
                    to_run = f"{mock_maker} {bbmap_dir}/./randomreads.sh {coverage} {support_size}"
                    res = subprocess.run(to_run,
                                         shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    if res.returncode != 0:
                        tb = traceback.format_exc()
                        error_type, error, _ = sys.exc_info()
                        raise Exception(f"Metagenome and reference creation failed due to {error}")

                    for q in q_vals:

                        for const in lambda_vals:

                            # Create y_small and y_large
                            with tempfile.NamedTemporaryFile() as temp_file:
                                output_file = temp_file.name

                                to_run = f"python ../../src/Form16SyVector.py -k {small_k} -i mock_16S_metagenome.fa -o {output_file}"
                                res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL,
                                                     stderr=subprocess.DEVNULL)
                                if res.returncode != 0:
                                    raise Exception("Failed to form y_small")

                                y_small = sio.loadmat(output_file)['y'].T

                                to_run = f"python ../../src/Form16SyVector.py -k {large_k} -i mock_16S_metagenome.fa -o {output_file}"
                                res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL,
                                                     stderr=subprocess.DEVNULL)
                                if res.returncode != 0:
                                    raise Exception("Failed to form y_large")

                                y_large = sio.loadmat(output_file)['y'].T

                            ## Full reference database calculations

                            if full_reference_genome:  # Do so only if location is entered

                                # Create true_x
                                try:
                                    true_x = make_true_x("mock_16S_metagenome.fa", full_reference_genome, bbmap_dir)
                                except:
                                    tb = traceback.format_exc()
                                    error_type, error, _ = sys.exc_info()
                                    raise Exception(f"error_type: {error_type}\n error: {error}")


                                # Reconstruct true_x with MinDivLP
                                x_star = MinDivLP(A_k_small_full.toarray(), A_k_large_full, y_small, y_large, const, q)

                                # Calculate error and support data
                                l1error = norm(true_x - x_star, 1)
                                supp_x_star = np.where(x_star > 0)[0]
                                supp_true_x = np.where(true_x > 0)[0]
                                support_in_common = np.intersect1d(supp_true_x, supp_x_star)

                                ## Append csv file with results
                                with open(csv_name, "a", newline="") as f:
                                    writer = csv.writer(f)
                                    writer.writerow(
                                        [small_k, large_k, support_size, coverage, q, const, 1, l1error,
                                         len(supp_x_star), len(support_in_common),
                                         calculateDiv(true_x, A_k_large_full, q),
                                         calculateDiv(x_star, A_k_large_full, q)])

                            ## Mock reference database calculations (database containing only sequences in mock metagenome)

                            if reference_genome:  # Do so only if location is entered

                                # Create true_x
                                true_x = make_true_x(f"mock_16S_metagenome.fa", reference_genome, bbmap_dir)

                                # Create A_k_small and A_k_large from mock reference database
                                with tempfile.NamedTemporaryFile() as temp_file:
                                    output_file = temp_file.name

                                    to_run = f"python ../../src/Form16SSensingMatrix.py -k {small_k} -i {reference_genome} -o {output_file}"
                                    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL,
                                                         stderr=subprocess.DEVNULL)
                                    if res.returncode != 0:
                                        raise Exception("Failed to form A_k_small")

                                    A_k_small = sio.loadmat(output_file)['A_k']

                                    to_run = f"python ../../src/Form16SSensingMatrix.py -k {large_k} -i {reference_genome} -o {output_file}"
                                    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL,
                                                         stderr=subprocess.DEVNULL)
                                    if res.returncode != 0:
                                        raise Exception("Failed to form A_k_large")

                                    A_k_large = sio.loadmat(output_file)['A_k']

                                # Reconstruct true_x with MinDivLP
                                x_star = MinDivLP(A_k_small.toarray(), A_k_large, y_small, y_large, const, q)

                                # Calculate error and support data
                                l1error = norm(true_x - x_star, 1)
                                supp_x_star = np.where(x_star > 0)[0]
                                supp_true_x = np.where(true_x > 0)[0]
                                support_in_common = np.intersect1d(supp_true_x, supp_x_star)

                                ## Append csv file with results
                                with open(csv_name, "a", newline="") as f:
                                    writer = csv.writer(f)
                                    writer.writerow(
                                        [small_k, large_k, support_size, coverage, q, const, 0, l1error,
                                         len(supp_x_star), len(support_in_common),
                                         calculateDiv(true_x, A_k_large, q), calculateDiv(x_star, A_k_large, q)])

                            # keep counts
                            i += 1
                            print(f"large_k = {large_k}: {i}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Uses bbmap to create simulated metagenomes and then benchmarks the MinDivLP for a variety of parameter values",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bbmap_loc', type=str, help="location of bbmap directory on your system", required=True)
    args = parser.parse_args()
    bbmap_dir = args.bbmap_loc
    if not os.path.exists(bbmap_dir):
        raise Exception(f"The directory {bbmap_dir} does not appear to exist. Please check and try again.")

    ## Save first row of csv if it hasn't been written
    if not os.path.exists(resultCSV):
        with open(resultCSV, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(
                ['small k', 'large k', 'support size', 'coverage', 'q', 'const', 'database', 'l1 error',
                 'constructed support size', 'common support size', 'true diversity', 'constructed diversity'])

    for large_k in large_k_vals:
        simulate(large_k, "./make_mock_metagenomes.sh", resultCSV, bbmap_dir)
