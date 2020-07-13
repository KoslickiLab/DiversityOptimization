#! /usr/bin/env python
import argparse
import os
import sys
import subprocess
import tempfile
import scipy.io as sio
from pandas import read_csv
from csv import writer
from time import time
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))  # make sure python knows where to find the code
from src.MinDivLP import MinDivLP
from src.ConvertXToTaxonomicProfile import convertToTaxonomy

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
    parser.add_argument('-p', '--prevent_output', action="store_true",
                        help="stop output", default=False)

    args = parser.parse_args()
    input_file = os.path.abspath(args.input_file)
    output_dir = os.path.abspath(args.output_dir)
    small_k = args.small_k
    large_k = args.large_k
    const = args.const
    q = args.q_value
    count_rev = args.count_complements
    reference = args.reference
    taxonomy = args.taxonomy
    prevent_output = args.prevent_output

    start = time()

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
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/Form16SSensingMatrix.py'))
        to_run = f"""python "{path}" -k {large_k} -i "{reference}" -o "{reference}_A_{large_k}.mat" {count_rev * '-c'}"""
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large sensing matrix")

    A_k_large = sio.loadmat(f"{reference}_A_{large_k}.mat")['A_k']

    # small_k
    if not os.path.exists(f"{reference}_A_{small_k}.mat"):
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/Form16SSensingMatrix.py'))
        to_run = f"""python "{path}" -k {small_k} -i "{reference}" -o "{reference}_A_{small_k}.mat" {count_rev * '-c'}"""
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small sensing matrix")

    A_k_small = sio.loadmat(f"{reference}_A_{small_k}.mat")['A_k']

    ## Create y vectors:

    with tempfile.NamedTemporaryFile() as temp_file:
        t_file = temp_file.name
        path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../src/Form16SyVector.py'))

        # large_k
        to_run = f"""python "{path}" -k "{large_k}" -i "{input_file}" -o "{t_file}" {count_rev * '-c'}"""
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form large y vector")

        y_large = sio.loadmat(t_file)['y'].T

        # small_k
        to_run = f"""python "{path}" -k {small_k} -i "{input_file}" -o "{t_file}" {count_rev * '-c'}"""
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to form small y vector")

        y_small = sio.loadmat(t_file)['y'].T

    ## Run MinDivLP

    print("Running MinDivLP and saving.")

    x = MinDivLP(A_k_small, A_k_large, y_small, y_large, const, q)

    if not prevent_output:
        ## Convert to TSV, then to BIOM
        tsv_file = os.path.join(output_dir, "taxonomy.tsv")
        biom_file = os.path.join(output_dir, "table.biom")

        # Check for metadata in sample-metadata.tsv and extract the sample ID
        if os.path.exists(os.path.join(os.path.dirname(input_file), "sample-metadata.tsv")):
            metadata_df = read_csv(os.path.join(os.path.dirname(input_file), "sample-metadata.tsv"), sep='\t', header=0)
            sample_id = metadata_df.iloc[0]["SampleID"]
        else:
            sample_id = "sample"

        convertToTaxonomy(x, reference, taxonomy, sample_id, tsv_file)

        to_run = f"""biom convert -i "{tsv_file}" -o "{biom_file}" --to-json --process-obs-metadata taxonomy"""
        res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
        if res.returncode != 0:
            raise Exception("Failed to convert TSV to BIOM format")

        print("Completed.")
        timer = time() - start

        # Print info useful for checking progress of code
        print([output_dir.split('/')[-4], small_k, large_k, const, q, timer])

        # Record times
        with open('times.csv', 'a') as fd:
            writer(fd).writerow([output_dir.split('/')[-4], small_k, large_k, const, q, timer])
