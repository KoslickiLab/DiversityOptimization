import numpy as np
import pandas as pd
import os
import csv


def convertToTaxonomy(x, trainingfasta, trainingtaxonomy, filename):
    """ convertToTaxonomy
        Creates or appends a tsv table with a taxonomic profile from population proportions, a reference database fasta
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
        A tsv table containing the taxonomic profiles OTU of samples already included in filename and the taxonomic
            profile of x
    """

    x = np.round(1000 * x)  # temporary constant to create counts rather than relative abundances

    ## Import tsv data frame or initialize empty dataframe
    if os.path.exists(filename):
        df = pd.read_csv(filename, delimiter='\t')
    else:
        df = pd.DataFrame(columns=["#NAME"])

    ## Determine OTUs that are in the new sample
    support = np.where(x > 0)[0]
    includedIDs = list()

    with open(trainingfasta, "r") as fasta:
        i = 0
        for line in fasta:
            if line[0] == ">":
                if i in support:
                    includedIDs.append(line[1:-1])
                i += 1

    ## Find taxa of nonzero OTUs and place them in the same order
    taxID_sample = includedIDs

    with open(trainingtaxonomy, "r") as tax:
        tax_reader = csv.reader(tax, delimiter='\t')

        for line in tax_reader:
            if line[0] in includedIDs:
                index = includedIDs.index(line[0])
                taxID_sample[index] = line[1]  # at same index as otu

                # indiv_tax = [a[3:] for a in line[1].split() if len(a) > 4]  # Don't include "k__", "p__", etc.
                # taxID_sample[index] = "".join(indiv_tax)

    ## Create new data frame
    new_data = list()
    indices_to_add = list(range(len(support)))  # Indices within the support of OTUs not included in original data frame

    # Add new column to existing rows
    for i in range(df.shape[0]):
        if df["#NAME"][i] in taxID_sample:  # Add abundance if in support
            index = taxID_sample.index(df["#NAME"][i])  # Index of already included otu in support
            del indices_to_add[index]
            new_data.append(list(df.loc[i]) + [x[support][index]])
        else:  # Add zero if not in  support
            new_data.append(list(df.loc[i]) + [0])

    # Add new rows
    for i in indices_to_add:
        new_data.append([taxID_sample[i]] + [0 for _ in range(df.shape[1] - 1)] + [x[support][i]])

    # Create new data frame
    new_df = pd.DataFrame(new_data, columns=list(df.columns) + [f"Sample{df.shape[1]}"])

    new_df.to_csv(filename, sep='\t', index=False)
