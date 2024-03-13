#!/usr/bin/env python
# count_and_print_kmers2.py
seq = 'GCCGGCCCTCAGACAGGAGTGGTCCTGGATGTGGATG'
kmer_length = 6

# Initialize a k-mer dictionary
kmer_dictionary = {}

stop = len(seq) - kmer_length + 1

# Iterate over the positions
for start in range(0, stop):
    # Get the substring at a specific start and end position
    kmer = seq[start:start + kmer_length]

    # Increase the count of this kmer by 1
    # if the kmer is new, start at 0 before incrementing by 1
    kmer_dictionary[kmer] = kmer_dictionary.get(kmer, 0) + 1 

# Iterate over the k-mers and counts in the dictionary
for kmer, count in kmer_dictionary.items():
    print("{0}\t{1}".format(kmer, count))
    # there's also an f-string option
    # print(f'{kmer}\t{count}')
