# Description
This directory contains files used for benchmarking MinDivLP against QIIME implemented classification methods with the classifier evaluation tool [tax-credit](https://github.com/caporaso-lab/tax-credit-code).

# classify_mindivlp.py
A command line script to create or load sensing matrices, create y vectors, call MinDivLP, and convert the output to TSV and BIOM formats.

# convertToTaxonomy
A function to convert the reconstructed population into TSV format with columns being OTU IDs, proportions, and taxonomies of a single sample

# taxonomy-assignment-mindivlp.ipynb
An iPython Notebook to run a parameter sweep of MinDivLP tests and store them with other classification method results. (stored in tax-credit-data/ipynb/mock-community)

# evaluate-classification-accuracy.ipynb
An iPython Notebook to evaluate qualities of classifications by various measures and plot comparisons by classification method and taxonomic level. (stored in tax-credit-data/ipynb/mock-community)