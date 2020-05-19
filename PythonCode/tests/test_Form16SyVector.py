import subprocess
import tempfile
import scipy.io as sio
import numpy as np
import os

# Revcomp
with tempfile.NamedTemporaryFile() as temp_file:
    output_file = temp_file.name

    to_run = f"../src/./Form16SyVector.py -k 4 -i ../data/mock_16S_metagenome.fa -o {output_file}"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Test failed to run")

    # run KMC
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.NamedTemporaryFile() as KMC_out_file:
            with tempfile.NamedTemporaryFile() as KMC_dump_file:
                to_run = f"kmc -b -fm -k4 -ci0 -cs1000000 ../data/mock_16S_metagenome.fa {os.path.join(temp_dir, KMC_out_file.name)} ."
                res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
                to_run = f"kmc_dump -ci0 -cs1000000 {os.path.join(temp_dir, KMC_out_file.name)} {os.path.join(temp_dir, KMC_dump_file.name)}"
                res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)

                # do the comparison
                test_result = sio.loadmat(output_file)['y']
                kmc_results = []
                with open(KMC_dump_file.name) as fid:
                    for line in fid.readlines():
                        kmc_results.append(int(line.strip().split()[1]))
                kmc_results = np.array(kmc_results)
                kmc_results = kmc_results/np.sum(kmc_results)
                thresh = 0.00001
                assert np.sum(np.abs(kmc_results.flatten() - test_result.flatten())) < thresh



print("Tests passed successfully!")

