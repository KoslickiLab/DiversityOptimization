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
x[np.random.choice(range(10), size=3)] = [1, 1, 1]
support, = np.where(x > 0)
convertToTaxonomy(x, reference_genome, reference_taxonomy, tsv_name)

# Gather possible OTU IDs for comparison
possible_otu_ids = list()
with open(reference_genome, "r") as fasta:
    i = 0
    for line in fasta:
        if line[0] == ">":
            possible_otu_ids.append(line[1:-1])
            i += 1
        if i >= 10:
            break

# Gather the taxID of each OTU
intended_taxa = np.zeros(10, dtype = object)
with open(reference_taxonomy, "r") as tax:
    tax_reader = csv.reader(tax, delimiter='\t')

    for line in tax_reader:
        if line[0] in possible_otu_ids:
            index = possible_otu_ids.index(line[0])
            intended_taxa[index] = line[1]

# Print each taxID in order
print("Intended taxonomies to be added or appended:")
print(intended_taxa[np.where(x>0)[0]])


# Print each taxID included in the tsv file and check for correctness
print("\nTSV table:")
correct = True
with open(tsv_name, "r") as tax:
    tax_reader = csv.reader(tax, delimiter='\t')
    next(tax_reader)  # Skip header

    for line in tax_reader:
        print(line)

        if line[0] in intended_taxa[support]:  # Check to make sure correct counts are added
            index_in_support, = np.where(intended_taxa[support] == line[0])[0]
            if line[-1] != str(100 * x[support][index_in_support]):  # Constant may need to change
                correct = False

        else:  # Check to make sure counts that aren't supposed to be added are not added
            if line[-1] != '0.0':
                correct = False

if correct:
    print("\nThe TSV table was appended correctly.")
else:
    print("\nError.")
