#! /usr/bin/env python
import argparse
import os
import subprocess
import numpy as np
from scipy.sparse import coo_matrix
import scipy.io as sio
from sklearn.preprocessing import normalize
import tempfile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Creates a y-vector (i.e. sample vector) when presented with a fasta or fastq input 16S metagenome.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-k', '--k_size', type=int,
                        help="k-mer size to use (Note: values >14 will probably take too long)")
    parser.add_argument('-c', '--count_complements', action="store_true",
                        help="count compliment of sequences as well", default=False)
    parser.add_argument('-i', '--input_file', type=str, help="File name of input database")
    parser.add_argument('-o', '--output_file', type=str,
                        help="Output file of the y-vector in .mat format.",
                        required=True)

    # read in the arguments
    args = parser.parse_args()
    k_size = args.k_size
    count_rev = args.count_complements
    input_file_name = args.input_file
    output_file_name = args.output_file

    # check if the input exists
    if not os.path.exists(input_file_name):
        raise Exception(f"The input file {input_file_name} does not appear to exist")

    # check if dna-utils is installed
    res = subprocess.run("kmer_counts_per_sequence -h", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception(
            "It appears that dna-utils is not installed. Please consult the README, install dna-utils, and try again.")

    # check if fastq or fasta
    with open(input_file_name, 'r') as fid:
        line = fid.readline()
        first_char = line[0]
        if first_char == '>':
            is_fasta = True
        else:
            is_fasta = False

    if count_rev:
        if is_fasta:
            to_run = f"kmer_total_count -i {input_file_name} -k {k_size} -c"
        else:
            to_run = f"sed -n '1~4s/^@/>/p;2~4p' {input_file_name} | kmer_total_count -k {k_size} -c"
    else:
        if is_fasta:
            to_run = f"kmer_total_count -i {input_file_name} -k {k_size}"
        else:
            to_run = f"sed -n '1~4s/^@/>/p;2~4p' {input_file_name} | kmer_total_count -k {k_size}"

    res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("An unexpected error was encountered, please check the input FASTA file is in the correct format. If errors persist, contact the developers.")

    y = np.fromstring(res.stdout.decode('utf-8'), sep='\n', dtype=int)
    y_norm = y / np.sum(y)
    sio.savemat(output_file_name, {"y": y_norm}, do_compression=True)
