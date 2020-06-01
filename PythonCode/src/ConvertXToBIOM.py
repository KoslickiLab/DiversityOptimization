import numpy as np
from biom.table import Table

## For now only one sample in one environment with an estimate of population from read counts
def convertToBiom(x, trainingfile, yfile, samplename):

    observe_ids = list()

    sample_reads = 0
    sequence_reads = list()

    # Create OTU ids with those provided in FASTA
    with open(trainingfile, 'r') as file:
        for line in file:
            if line[0] == '>':
                observe_ids.append(line[1:-1])
            else:  # Count the number of mers in each full sequence
                sequence_reads.append(len(line)-1)  # Accounting for '>'

    with open(yfile, 'r') as file:
        for line in file:
            if line[0] != '>':  # Count the number of mers in the population
                sample_reads += (len(line)-1)  # Accounting for '>'

    unit_reads = np.array(sequence_reads) @ x  # unit_reads is the count of reads per unit population
    const = sample_reads / unit_reads  # Calculate the number of unit populations to reach the number of reads in the
    # sample. Likely an underestimate; number of reads in the population should be greater than or equal to the
    # number of reads in the sample

    counts = np.round(const * x).reshape((-1, 1))

    return Table(counts, observe_ids, [samplename])
