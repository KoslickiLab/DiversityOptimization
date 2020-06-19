import numpy as np
import pandas as pd
import csv


def convertToTaxonomy(x, trainingfasta, trainingtaxonomy, filename):
    """ convertToTaxonomy
        Creates a tsv table with a taxonomic profile from population proportions, a reference database fasta
        file, and the database's taxonomic identifiers
        Call via:
        table = convertToBiom(x, trainingfasta, trainingtaxonomy, filename)
        Parameters are:
        x is the vector of population proportions
        trainingfasta is the string location of the reference database fasta file (assumed to list OTU id after '>' for
            the following sequence, )
        trainingtaxonomy is the string location of the txt file containing taxonomic IDs of each OTU/sequence in
            trainingfasta (assumed to list one OTU/sequence id tab separated with the corresponding taxonomic id on
            each line)
        filename is the string location of the tsv file to be saved

        Outputs:
        A tsv table containing the taxonomic profile of x
    """

    ## Determine OTUs that are in the new sample
    support = np.where(x > 0)[0]
    sample_OTUs = list()

    with open(trainingfasta, "r") as fasta:
        i = 0
        for line in fasta:
            if line[0] == ">":
                if i in support:
                    sample_OTUs.append(line[1:-1])
                i += 1

    ## Find taxa of nonzero OTUs and place them in the same order
    sample_taxIDs = sample_OTUs.copy()

    with open(trainingtaxonomy, "r") as tax:
        tax_reader = csv.reader(tax, delimiter='\t')

        for line in tax_reader:
            if line[0] in sample_OTUs:
                index = sample_OTUs.index(line[0])
                sample_taxIDs[index] = line[1]  # at same index as OTU

    ## Create data frame
    data = [[sample_OTUs[i], x[support][i], sample_taxIDs[i]] for i in range(len(support))]

    df = pd.DataFrame(data, columns=["#OTU ID", "MockHiSeq.even", "taxonomy"])

    df.to_csv(filename, sep='\t', index=False)
