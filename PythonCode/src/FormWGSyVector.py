#! /usr/bin/env python
import argparse
import os
import subprocess
import tempfile
import itertools
import numpy as np
from CMash import MinHash as MH
from scipy.sparse import csc_matrix, save_npz

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a y-vector (i.e. sample vector) when presented with a fasta or fastq input WGS metagenome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-k', '--k_size', type=int,
                        help="k-mer size to use")
    parser.add_argument('-m', '--max_ram', type=int,
                        help="max amount of RAM in GB", default=12)
    parser.add_argument('--ci', type=int,
                        help="minimum count of appearances for a k-mer to be included", default=0)
    parser.add_argument('--cs', type=int,
                        help="maximum count of appearances recorded for a k-mer", default=256)
    parser.add_argument('-c', '--count_complements', action="store_true",
                        help="count compliment of sequences as well", default=False)
    parser.add_argument('-i', '--input_file', type=str, help="File name of input data")
    parser.add_argument('-t', '--training_prefix', type=str,
                        help="File path to training files (only prefix)", required=True)
    parser.add_argument('-o', '--output_file', type=str,
                        help="Output file of the y-vector in .mat format.",
                        required=True)

    ## Read in the arguments
    args = parser.parse_args()
    k_size = args.k_size
    max_ram = args.max_ram
    ci = args.ci
    cs = args.cs
    count_rev = args.count_complements
    input_file_name = args.input_file
    output_file_name = args.output_file
    training_prefix = args.training_prefix

    ## Existence checks (input, kmc, kmc_tools, kmc_dump)

    # check if the input exists
    if not os.path.exists(input_file_name):
        raise Exception(f"The input file {input_file_name} does not appear to exist")

    # check if kmc is installed
    res = subprocess.run("kmc", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception(
            "It appears that kmc is not installed. Please consult the README, install kmc, and try again.")

    # check if kmc_tools is installed
    res = subprocess.run("kmc_tools", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception(
            "It appears that kmc_tools is not installed. Please consult the README, install kmc_tools, and try again.")

    # check if kmc_dump is installed
    res = subprocess.run("kmc_dump", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception(
            "It appears that kmc_dump is not installed. Please consult the README, install kmc_dump, and try again.")

    ## Run KMC and intersect
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.NamedTemporaryFile() as kmc_output:

            # Count k-mers (note: the output is named f"{kmc_output.name}.pre" and f"{kmc_output.name}.suf")
            to_run = f"kmc -k{k_size} {~count_rev * '-b '}-ci{ci} -cs{cs} -fm -m{max_ram} {input_file_name} {kmc_output.name} {temp_dir}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception(
                    "An unexpected error was encountered while running kmc, please check the input FASTA file is in the "
                    "correct format. If errors persist, contact the developers.")

            with tempfile.NamedTemporaryFile() as intersect_file:

                # Intersect with training kmers, keeping counts of sample kmers
                to_run = f"kmc_tools simple {training_prefix} {kmc_output.name} intersect {intersect_file.name} -ocright"
                res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
                if res.returncode != 0:
                    raise Exception("An unexpected error was encountered while running kmc_tools simple.")

                # Dump the intersected counted k-mers for reading into matrix
                to_run = f"kmc_dump -ci0 -cs100000 {intersect_file.name} {intersect_file.name}"
                res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE)
                if res.returncode != 0:
                    raise Exception("An unexpected error was encountered while running kmc_tools.")

                ## Load in list of k-mers
                database = MH.import_multiple_from_single_hdf5(training_prefix + ".h5")
                kmers = sorted(set(itertools.chain.from_iterable(genome._kmers for genome in database)))

                ## Iterate through KMC's dump file to extract k-mers and their counts while determining
                #  their corresponding index
                indices = []
                data = []
                for line in intersect_file:
                    info = line.split()
                    try:
                        indices.append(kmers.index(info[0].decode("utf-8")))
                        data.append(int(info[1]))
                    except ValueError:
                        print("k-mer mismatch") ## TODO: Identify why there are so many mismatches


    ## Sort the indices and data
    sorter = sorted(range(len(indices)), key=indices.__getitem__)
    indices = np.array([indices[i] for i in sorter])
    data = np.array([data[i] for i in sorter])
    data = data / np.sum(data)
    indptr = np.array([0, len(indices)])

    ## Create and save a csc_matrix as .npz file
    save_npz(output_file_name, csc_matrix((data, indices, indptr), shape=(len(indices), 1)), compressed=True)
