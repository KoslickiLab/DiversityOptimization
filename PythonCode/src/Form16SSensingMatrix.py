#! /usr/bin/env python
import argparse
import os
import subprocess
import numpy as np
from scipy.sparse import coo_matrix
import scipy.io as sio
from sklearn.preprocessing import normalize

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description="Will form the sensing matrix `A` when given a database of microbial 16S genomes in FASTA format.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-k', '--k_size', type=int, help="k-mer size to use (Note: values >14 will probably take too long)")
	parser.add_argument('-c', '--count_complements', action="store_true",
						help="count compliment of sequences as well", default=False)
	parser.add_argument('-i', '--input_file', type=str, help="File name of input database")
	parser.add_argument('-o', '--output_file', type=str,
						help="Output file of sparse representation of sensing matrix `A` in .mat format.", required=True)

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
		raise Exception("It appears that dna-utils is not installed. Please consult the README, install dna-utils, and try again.")

	if count_rev:
		res = subprocess.run(f"kmer_counts_per_sequence -i {input_file_name} -k {k_size} -c -s", shell=True, stdout=subprocess.PIPE)
	else:
		res = subprocess.run(f"kmer_counts_per_sequence -i {input_file_name} -k {k_size} -s", shell=True, stdout=subprocess.PIPE)

	if res.returncode == 0:
		print("Sparse matrix created, now converting to .mat format")
	else:
		print("An unexpected error was encountered, please check the input FASTA file is in the correct format. If errors persist, contact the developers.")

	sparse_text_form = res.stdout.decode('utf-8')  # convert from bytes to string
	J, I, V = np.array(list(map(int,sparse_text_form.split()))).reshape((-2, 3)).transpose()  # pull out the indicies
	sparse_matrix_form = coo_matrix((V, (I, J)), shape=(4 ** k_size, J[-1] + 1))  # convert to sparse matrix format
	sparse_matrix_form_norm = normalize(sparse_matrix_form, norm='l1', axis=0)  # normalize by rows
	sio.savemat(output_file_name, {"A_k": sparse_matrix_form_norm}, do_compression=True)
