import sys
import os
import csv
import numpy as np

sys.path.append(os.path.abspath("../.."))
from PythonCode.src.ConvertXToTaxonomicProfile import convertToTaxonomy

## This test runs convertToTaxonomy and verifies that the last column is added correctly

# Locations of files used
reference_genome = "../data/97_otus.fasta"
reference_taxonomy = "../data/97_otu_taxonomy.txt"
tsv_name = "test.tsv"

# Create a randomized x value
x = np.zeros(99322)
x[np.random.choice(range(10), size=3, replace=False)] = [1, 1, 1]
support, = np.where(x > 0)

# Gather possible OTU IDs for comparison
possible_otu_ids = list()
with open(reference_genome, "r") as fasta:
    i = 0
    for line in fasta:
        if line[0] == ">":
            possible_otu_ids.append(line[1:-1])
            i += 1
        if i > max(support):  # End when beyond support
            break

# Gather taxIDs of each OTU
possible_taxa = np.zeros(max(support)+1, dtype = object)
with open(reference_taxonomy, "r") as tax:
    tax_reader = csv.reader(tax, delimiter='\t')

    for line in tax_reader:
        if line[0] in possible_otu_ids:
            index = possible_otu_ids.index(line[0])
            possible_taxa[index] = line[1]

# Gather taxIDs present in tsv
previous_taxa = list()

if os.path.exists(tsv_name):
    with open(tsv_name, "r") as tax:
        tax_reader = csv.reader(tax, delimiter='\t')
        next(tax_reader)  # Skip header

        for line in tax_reader:
            previous_taxa.append(line[-1])
        

# Print each taxID in sample in order
print("Intended taxonomies to be added or appended:")
print(possible_taxa[support])
print("Intended OTU IDs to be added or appended:")
print(np.array(possible_otu_ids)[support])
print("Support of x: ", support)


## Run function
convertToTaxonomy(x, reference_genome, reference_taxonomy, tsv_name)

# Print each taxID included in the tsv file and check for correctness
print("\nTSV table:")

with open(tsv_name, "r") as tax:
    tax_reader = csv.reader(tax, delimiter='\t')
    next(tax_reader)  # Skip header

    for line in tax_reader:
        print(line)

        if line[-1] in possible_taxa[support]:  # Check to make sure correct counts are added
            index_in_support, = np.where(possible_taxa[support] == line[-1])[0]
            assert line[-2] == str(100 * x[support][index_in_support]), "Nonzero count not added correctly"

            if not line[-1] in previous_taxa:  # Check to make sure new rows only add nonzero count to final column
                for val in line[1:-2]:
                    assert val == '0.0', "Nonzero count added to previous samples"

        else:  # Check to make sure counts that aren't supposed to be added are not added
            assert line[-2] == '0.0', "Nonzero count added that should be zero"

print("\nThe TSV table was appended correctly.")
