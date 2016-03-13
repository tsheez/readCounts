Read Counts
by Tim Shaw

Usage
readcounts.py [Input File Location] [Output File Location] [Length of 3' Random-mer]

This program was created for the Lykke-Andersen lab at UCSD for the purpose of processing sequencing results so they can
used by other downstream scripts.

readCounts takes an input FASTQ file, trims the random-mer from the 3' end, and counts how many copies of each read there
are. In order to remove PCR artifacts, it also compares the 3' random-mer. Matching 3' random-mers are considered to be
PCR artifacts and don't contribute to the number of unique reads.

Output
read without random-mer, total number of reads, number of unique reads