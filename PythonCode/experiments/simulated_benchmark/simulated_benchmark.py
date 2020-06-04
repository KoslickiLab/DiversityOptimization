import sys
import os
import scipy.io as sio
import tempfile
import subprocess

sys.path.append(os.path.abspath("../../.."))
from PythonCode.experiments.simulated_benchmark.make_true_x import make_true_x
from PythonCode.src.MinDivLP import MinDivLP
from scipy.linalg import norm

# Constants
q = .1
const = 10000

# k values
small_k = 4  # smaller k-mer size
large_k = 6  # larger k-mer size

## Use noiseless then noisy mock metagenomes
for mock_maker in ["make_mock_metagenomes.sh"]:  # TODO: create make_noisy_mock_metagenomes.sh or option in make_mock_metagenomes.sh

    ## Create ground truth mock metagenome and mock reference
    res = subprocess.run(f"bash {mock_maker} ./bbmap/randomreads.sh", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Metagenome and reference creation failed")

    # Create y_small and y_large
    with tempfile.NamedTemporaryFile() as temp_file:
        output_file = temp_file.name

        to_run = f"py ../../src/Form16SyVector.py -k {small_k} -i mock_16S_metagenome.fa -o {output_file}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form y_small")

        y_small = sio.loadmat(output_file)['y'].T

        to_run = f"py ../../src/Form16SyVector.py -k {large_k} -i mock_16S_metagenome.fa -o {output_file}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form y_large")

        y_large = sio.loadmat(output_file)['y'].T

    ## Use mock_reference database then the full reference database
    for reference_genome in ["./mock_reference.fa", "../../data/97_otus_subset.fasta"]:

        # Create true_x
        true_x = make_true_x("mock_16S_metagenome.fa", reference_genome)

        # Create A_k_small and A_k_large from reference database
        with tempfile.NamedTemporaryFile() as temp_file:
            output_file = temp_file.name

            to_run = f"py ../../src/Form16SSensingMatrix.py -k {small_k} -i {reference_genome} -o {output_file}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception("Failed to form A_k_small")

            A_k_small = sio.loadmat(output_file)['A_k']

            to_run = f"py ../../src/Form16SSensingMatrix.py -k {large_k} -i {reference_genome} -o {output_file}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception("Failed to form A_k_large")

            A_k_large = sio.loadmat(output_file)['A_k']

        # Reconstruct true_x with MinDivLP
        x_star = MinDivLP(A_k_small.toarray(), A_k_large, y_small, y_large, const, q)

        # Print error
        print(f"l1 error = {norm(true_x - x_star,1)}")
        print(f"l2 error = {norm(true_x - x_star,2)}")
