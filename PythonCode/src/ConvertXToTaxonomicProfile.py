import numpy as np
import pandas as pd
import os
import csv


def convertToTaxonomy(x, trainingfasta, trainingtaxonomy, sample_id, filename, append=True):
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

    ## Import tsv data frame or initialize empty dataframe
    if os.path.exists(filename) and append:
        df = pd.read_csv(filename, delimiter='\t')
    else:
        df = pd.DataFrame(columns=["#OTU ID", "taxonomy"])

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

    ## Create new data frame
    new_data = list()
    indices_to_add = list(range(len(support)))  # Indices within the support of OTUs not included in original data frame

    # Add new column to existing rows
    for i in range(df.shape[0]):
        current_row = list(df.loc[i])
        if str(df["#OTU ID"][i]) in sample_OTUs:  # Add abundance if in support
            index = sample_OTUs.index(str(df["#OTU ID"][i]))  # Index of already included OTU in support
            indices_to_add.remove(index)
            new_data.append(current_row[:-1] + [x[support][index]] + [current_row[-1]])
        else:  # Add zero if not in  support
            new_data.append(current_row[:-1] + [0] + [current_row[-1]])

    # Add new rows
    for i in indices_to_add:
        new_data.append([sample_OTUs[i]] + [0 for _ in range(df.shape[1] - 2)] + [x[support][i]] + [sample_taxIDs[i]])

    # Create new data frame
    column_names = list(df.columns)
    new_df = pd.DataFrame(new_data, columns=column_names[:-1] + [sample_id] + [column_names[-1]])

    new_df.to_csv(filename, sep='\t', index=False)
