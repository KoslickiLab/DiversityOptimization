import subprocess
import tempfile
import scipy.io as sio
import numpy as np

# No revcomp
with tempfile.NamedTemporaryFile() as temp_file:
    output_file = temp_file.name

    to_run = f"../src/./Form16SSensingMatrix.py -k 4 -i ../data/97_otus_subset.fasta -o {output_file}"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Test failed to run")

    known_result = sio.loadmat("../data/97_otus_subset.fasta_A_4.mat")
    test_result = sio.loadmat(output_file)
    thresh = 0.00001
    assert np.sum(np.abs(known_result['A_k'] - test_result['A_k']))



print("Tests passed successfully!")

