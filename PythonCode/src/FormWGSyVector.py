#! /usr/bin/env python
import argparse
import os
import subprocess
import tempfile
import numpy as np
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
                        help="minimum count of appearances for a k-mer to be included", default=2)
    parser.add_argument('--cs', type=int,
                        help="maximum count of appearances recorded for a k-mer", default=256)
    parser.add_argument('-c', '--count_complements', action="store_true",
                        help="count compliment of sequences as well", default=False)
    parser.add_argument('-i', '--input_file', type=str, help="File name of input database")
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

    # check if fastq or fasta
    with open(input_file_name, 'r') as fid:
        line = fid.readline()
        first_char = line[0]
        if first_char == '>':
            is_fasta = True
        else:
            is_fasta = False

    ## Run KMC
    with tempfile.TemporaryDirectory() as temp_dir:
        with tempfile.NamedTemporaryFile() as kmc_output:

            # Count k-mers (note: the output is named f"{kmc_output.name}.pre" and f"{kmc_output.name}.suf")
            if is_fasta:
                to_run = f"kmc -k{k_size} {~count_rev * '-b '}-ci{ci} -cs{cs} -fa -m{max_ram} {input_file_name} {kmc_output.name} {temp_dir}"
            else:
                to_run = f"kmc -k{k_size} {~count_rev * '-b '}-ci{ci} -cs{cs} -m{max_ram} {input_file_name} {kmc_output.name} {temp_dir}"

            res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception(
                    "An unexpected error was encountered while running kmc, please check the input FASTA file is in the "
                    "correct format. If errors persist, contact the developers.")

            # Sort the counted k-mers
            to_run = f"kmc_tools transform {kmc_output.name} sort {kmc_output.name}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception("An unexpected error was encountered while running kmc_tools.")

            # Dump the counted k-mers for reading into matrix
            to_run = f"kmc_dump -ci0 -cs100000 {kmc_output.name} {kmc_output.name}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception("An unexpected error was encountered while running kmc_tools.")


            ## Iterate through KMC's dump file to extract k-mers and their counts while converting k-mers to
            #  their corresponding index
            bytetoB4 = {65: 0, 67: 1, 71: 2, 84: 3}  # dict to convert byte values of ACGT to their base 4 representation
            indices = []
            data = []
            for line in kmc_output:
                info = line.split()
                index = 0
                for i in range(k_size):
                    index += bytetoB4[info[0][i]] * 4 ** (k_size - (i + 1))  # Convert base4 (ACTG <-> 0123) to base 10
                indices.append(index)
                data.append(int(info[1]))


    ## Convert extracted data to csc_matrix (numpy arrays are infeasible)
    indices = np.array(indices)
    data = np.array(data)
    data = data / np.sum(data)
    indptr = np.array([0, len(indices)])

    y = csc_matrix((data, indices, indptr), shape=(4 ** k_size, 1))

    ## Save the csc_matrix as .npz file (.mat does not appear to work well with very large k_size)
    save_npz(output_file_name, y, compressed=True)
