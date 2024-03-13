
"""
In this program we are finding the first Open Reading Frame in the given fasta file.
In this Program we are receiving dmel-all_chromosome-r6.17.fasta file and finding first
open reading frame
"""

import argparse
import re
from Bio.Seq import Seq
from Bio import SeqIO


def get_args():
    """Return parsed command-line arguments."""
    # Create the parser object
    parser = argparse.ArgumentParser(
        description="This program translates the first ORF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get the FASTA file of sequences
    # Add arguments to the parser
    parser.add_argument('filename',  # variable to access this data later: args.filename
                        metavar='FASTA',  # shorthand to represent the input value
                        help='Provide name and path to FASTA file to process.',
                        # message to the user, it goes into the help menu
                        type=str)
    parser.add_argument('-p', '--pattern',  # access with args.pattern
                        help='Provide a regex pattern for filtering FASTA entries',
                        default='^\d{1}\D*$')  # default works for Drosophila chromosomes

    return parser.parse_args()


def find_first_orf(rna):
    """Return first open-reading frame of RNA sequence as a Bio.Seq object.

    Must start with AUG
    Must end with UAA, UAG, or UGA
    Must have even multiple of 3 RNA bases between
    """
    try:
        # Update regex to find the ORF
        orf = re.search('AUG([AUGC]{3})+?(UAA|UAG|UGA)', str(rna)).group()
    except AttributeError:  # if no match found, orf should be empty
        orf = ""
    return Seq(orf)


def translate_first_orf(dna):
    """This function transcribe the given dna sequence to rna and finds the
     first orf and translate it into protein
    Assumes input sequences is a Bio.Seq object.
    """

    # Transcribe the DNA, find the first ORF, translate said ORF
    rna = dna.transcribe()

    orf = find_first_orf(rna)

    protein = orf.translate()
    return protein


if __name__ == "__main__":
    # Get command-line arguments
    args = get_args()
    # print(args.filename)

    for record in SeqIO.parse(args.filename, "fasta"):
        # re.search  matches object with some special character ^ starts with, \d to find a numeric digit
        # {1} one time and $ ends with
        if re.match("^\d{1}\D*$", record.id):
            protein = translate_first_orf(record.seq)
        # if the FASTA record's ID matches the regex pattern,
        # then print out its record ID then a tab space then the translated first ORF
            print("%s:\t %s" % (record.id, protein))

