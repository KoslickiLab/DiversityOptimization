#! /usr/bin/env python
import argparse
import os
import subprocess
import tempfile
import scipy.io as sio
from MinDivLP import MinDivLP  # Location may need adjustment
from convertToTaxonomy import convertToTaxonomy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Reconstructs population proportions of a sample",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_file', type=str, help="File name of input database", required=True)
    parser.add_argument('-o', '--output_dir', type=str,
                        help="Output file directory.",
                        required=True)
    parser.add_argument('-s', '--small_k', type=int,
                        help="small k-mer size to use (Note: values >14 will probably take too long)", default=6)
    parser.add_argument('-l', '--large_k', type=int,
                        help="large k-mer size to use (Note: values >14 will probably take too long)", default=12)
    parser.add_argument('-c', '--const', type=int, help="lambda (AKA const) value", default=10000)
    parser.add_argument('-q', '--q_value', type=float, help="q value", default=0.1)
    parser.add_argument('--count_complements', action="store_true",
                        help="count compliment of sequences as well", default=False)
    parser.add_argument('-r', '--reference', type=str, help="File name of reference database", required=True)
    parser.add_argument('-t', '--taxonomy', type=str, help="File name of reference taxonomy", required=True)

    args = parser.parse_args()
    input_file = args.input_file
    output_dir = args.output_dir
    small_k = args.small_k
    large_k = args.large_k
    const = args.const
    q = args.q_value
    count_rev = args.count_complements
    reference = args.reference
    taxonomy = args.taxonomy

    ## Check if input files exist:

    # input
    if not os.path.exists(input_file):
        raise Exception(f"The input file {input_file} does not appear to exist")
    # reference
    if not os.path.exists(reference):
        raise Exception(f"The reference reference file {reference} does not appear to exist")
    # taxonomy
    if not os.path.exists(taxonomy):
        raise Exception(f"The reference taxonomy file {taxonomy} does not appear to exist")

    ## Create sensing matrices if they have not been created, and load them

    # large_k
    if not os.path.exists(f"{reference}_A_{large_k}.mat"):
        to_run = f"python /media/sf_Shared/MinDivLP/src/./Form16SSensingMatrix.py -k {large_k} -i {reference} -o {reference}_A_{large_k}.mat {count_rev*'-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large sensing matrix")

    A_k_large = sio.loadmat(f"{reference}_A_{large_k}.mat")['A_k']

    # small_k
    if not os.path.exists(f"{reference}_A_{small_k}.mat"):
        to_run = f"python /media/sf_Shared/MinDivLP/src/./Form16SSensingMatrix.py -k {small_k} -i {reference} -o {reference}_A_{small_k}.mat {count_rev*'-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small sensing matrix")

    A_k_small = sio.loadmat(f"{reference}_A_{small_k}.mat")['A_k']

    ## Create y vectors:

    with tempfile.NamedTemporaryFile() as temp_file:
        t_file = temp_file.name

        # large_k
        to_run = f"python /media/sf_Shared/MinDivLP/src/./Form16SyVector.py -k {large_k} -i {input_file} -o {t_file} {count_rev * '-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large y vector")

        y_large = sio.loadmat(t_file)['y'].T

        # small_k
        to_run = f"python /media/sf_Shared/MinDivLP/src/./Form16SyVector.py -k {small_k} -i {input_file} -o {t_file} {count_rev * '-c'}"
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small y vector")

        y_small = sio.loadmat(t_file)['y'].T
		

    ## Run MinDivLP

    x = MinDivLP(A_k_small.toarray(), A_k_large, y_small, y_large, const, q)
	
    ## Convert to TSV, then to BIOM
    tsv_file = f"{output_dir}/rep_seqs_tax_assignments.txt"

    convertToTaxonomy(x, reference, taxonomy, tsv_file)

    to_run = f"biom convert -i {tsv_file} -o {output_dir}/table.biom --to-json --process-obs-metadata taxonomy"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Failed to convert TSV to BIOM format")