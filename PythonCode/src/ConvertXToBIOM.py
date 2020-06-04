import numpy as np
from biom.table import Table

## For now only one sample in one environment with an estimate of population from read counts
def convertToBiom(x, trainingfile, samplefile, samplename):
    """ convertToBiom
        Creates a BIOM table from population proportions, a reference database fasta file, and a sample fasta file
        Call via:
        table = convertToBiom(x, trainingfile, samplefile, samplename)
        Parameters are:
        x is the vector of population proportions
        training file is the string location of the reference database fasta file (assumed to list OTU id after '>' for
            the following sequence)
        samplefile is the string location of the sample fasta file of the sample
        samplename is the name of the given sample

        Returns:
        table: a BIOM table with estimates

    """
    observe_ids = list()

    sample_len = 0
    persequence_len = list()

    # Create OTU ids with those provided in FASTA
    with open(trainingfile, 'r') as file:
        for line in file:
            if line[0] == '>':
                observe_ids.append(line[1:-1])
            else:  # Count the number of bases in each full sequence, i.e. the length of each sequence
                persequence_len.append(len(line)-1)  # Accounting for '>'

    with open(samplefile, 'r') as file:
        for line in file:
            if line[0] != '>':  # Count the number of bases in the population
                sample_len += (len(line)-1)  # Accounting for '>'

    unit_len = np.array(persequence_len) @ x  # unit_len is the count of bases per normalized population
    const = sample_len / unit_len  # Calculate the number of unit populations to reach the number of bases in the
    # sample, sample_len. Likely an underestimate; number of bases in the population should be greater than or equal
    # to the number of bases in the sample

    counts = np.round(const * x).reshape((-1, 1))

    return Table(counts, observe_ids, [samplename])
