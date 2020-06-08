#!/bin/bash
set -e

# This script uses bbmap to create simulated metagenomes.
# Usage: make_mock_metagenomes.sh <location of randomreads.sh>
# Output: 1. mock_16S_metagenome.fa: FASTA file of the simulated mock metagenome
#         2. mock_16S_metagenome.fastq: FASTQ file of the simulated mock metagenome
#         3. mock_reference.fa: the reference database used to generate the simulated mock metagenomes (curently set to the first numOrganisms in bbmapRef in the source code)

# requires bbmap: https://sourceforge.net/projects/bbmap/

#bbmapRandomReads="~/Documents/bbmap/randomreads.sh"
bbmapRandomReads=$1  # read the first input to this script, which should be the location of randomreads.sh

if [ -f "$bbmapRandomReads" ]; then
    echo "$bbmapRandomReads exist"
else
    echo "$bbmapRandomReads does not exist.
Usage: make_mock_metagenomes.sh <location of randomreads.sh>
Example usage: ./make_mock_metagenomes.sh ~/Documents/bbmap/randomreads.sh <coverage> <support_size>
Output: 1. mock_16S_metagenome.fa: FASTA file of the simulated mock metagenome
        2. mock_16S_metagenome.fastq: FASTQ file of the simulated mock metagenome
        3. mock_reference.fa: the reference database used to generate the simulated mock metagenomes (curently set to the first numOrganisms in bbmapRef in the source code)"
    exit 1
fi


bbmapRef='mock_reference.fa'
referenceDatabase='../../data/97_otus_subset.fasta'
mockCommunity="mock_16S_metagenome.fa"
coverage=$2
numOrganisms=$3
expDist="f"  # if "f", than an even-ish distribution will be used for x. If "t", then an exponential-ish distribution will be used for x


# Next, make a silly mock community
echo "mock community organisms:"
# TODO: first, randomly shuffle the reference so the order is different each time, or make this a command line option. See bbmap/shuffle.sh
head -n $((2 * ${numOrganisms})) ${referenceDatabase} >> "${bbmapRef}_temp"
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "${bbmapRef}_temp" | tail -n +2 > "${bbmapRef}"
rm "${bbmapRef}_temp"
# FASTA
eval ${bbmapRandomReads} ref=${bbmapRef} out="${mockCommunity}_temp.fa" metagenome=f adderrors=false simplenames=t minlength=250 maxlength=250 maxsnps=0 coverage=${coverage} maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0 maxinslen=0 maxdellen=0 maxsublen=0 maxnlen=0 mininslen=0 mindellen=0 minsublen=0 minnlen=0
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < "${mockCommunity}_temp.fa" | tail -n +2 > "${mockCommunity}"
# FASTQ
eval ${bbmapRandomReads} ref=${bbmapRef} out="${mockCommunity%.fa}.fastq" metagenome=f adderrors=false simplenames=t minlength=250 maxlength=250 maxsnps=0 coverage=${coverage} maxsnps=0 maxinss=0 maxdels=0 maxsubs=0 maxns=0 maxinslen=0 maxdellen=0 maxsublen=0 maxnlen=0 mininslen=0 mindellen=0 minsublen=0 minnlen=0
rm "${mockCommunity}_temp.fa"
rm -rf ref 2> /dev/null
#rm -rf organism_files 2> /dev/null
