#!/usr/bin/python
import argparse
import os
import subprocess
import tempfile
from itertools import chain
from scipy.sparse import csc_matrix, save_npz
from CMash import MinHash as MH

## TODO: adjust for very large files by only loading portions of hdf5 database file at a time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Forms the sensing matrix `A` when given a database of WGS genomes in FASTA format. Also"
                    "runs KMC on a temporary FASTA file consisting of each kmer counted by CMash for intersection"
                    "with sample FASTA files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-k', '--k_size', type=int,
                        help="k-mer size to use (Default: 21)",
                        default=21)
    parser.add_argument('-c', '--cmash_loc', type=str,
                        help="Location of CMash for accessing MakeStreamingDNADatabase.py.",
                        required=True)
    parser.add_argument('-i', '--input_file', type=str,
                        help="Path to file containing absolute file names of training genomes",
                        required=True)
    parser.add_argument('-o', '--output_file', type=str,
                        help="Output file of sparse representation of sensing matrix `A` in .mat format.",
                        default="TrainingDatabase",
                        required=True)


    ## Read in the arguments
    args = parser.parse_args()
    k_size = args.k_size
    input_file_names = os.path.abspath(args.input_file)
    if not os.path.exists(input_file_names):
        raise Exception("Input file %s does not exist." % input_file_names)
    output_file_name = os.path.abspath(args.output_file)
    cmash_loc = args.cmash_loc

    ## Run CMash and load database
    print("Running CMash.")

    database_file_name = output_file_name + ".h5"
    to_run = f"python {os.path.join(cmash_loc, 'scripts/MakeStreamingDNADatabase.py')} -k {k_size} {input_file_names} {database_file_name}"
    res = subprocess.run(to_run, shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception("Error running CMash (MakeStreamingDNADatabase.py)")

    print("Loading in database.")

    database = MH.import_multiple_from_single_hdf5(database_file_name)

    ## Extract list of kmers
    print("Gathering kmers.")

    # (convert to set for removing duplicates, sort for consistent ordering)
    kmers = sorted(set(chain.from_iterable(genome._kmers for genome in database)))

    ## Gather data for, create, and save csc_matrix, each column representing the counts of a particular genome's kmers
    print("Creating matrix.")

    data = []
    indices = []
    indptr = [0]

    for genome in database:
        column_indices = [kmers.index(kmer) for kmer in genome._kmers]  # Find indices of a particular genome's kmers
        sorter = sorted(range(len(column_indices)), key=column_indices.__getitem__)  # argsort the indices
        # Concatenate sorted indices and corresponding data (counts) for canonical csc_matrix format
        indices += [column_indices[i] for i in sorter]
        data += [genome._counts[i] for i in sorter]
        indptr.append(len(indices))

    save_npz(output_file_name, csc_matrix((data, indices, indptr), shape=(len(kmers), len(database))), compressed=True)

    ## Run KMC on artificial FASTA file for future intersection with y vectors

    print("Creating FASTA and running KMC.")

    # Check if kmc is installed
    res = subprocess.run("kmc", shell=True, stdout=subprocess.DEVNULL)
    if res.returncode != 0:
        raise Exception(
            "It appears that kmc is not installed. Please consult the README, install kmc, and try again.")

    with tempfile.TemporaryDirectory() as temp_dir:
        with open("test.fasta", "w+") as training_out_file:

            to_write = ""
            iterator = 0
            for kmer in kmers:
                to_write += ">seq%d\n" % iterator
                to_write += "%s\n" % kmer
                iterator += 1
            training_out_file.write(to_write)

            # Run KMC on training fasta
            to_run = f"kmc -k{k_size} -sm -fm -ci0 -cs3 {training_out_file.name} {output_file_name} {temp_dir}"
            res = subprocess.run(to_run, shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            if res.returncode != 0:
                raise Exception("An unexpected error was encountered while running kmc.")

    print("Completed.")
