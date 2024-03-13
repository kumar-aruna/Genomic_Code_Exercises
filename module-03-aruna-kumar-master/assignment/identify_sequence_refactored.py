#!/usr/bin/env python
"""
Identifying the given file consists aminoacid/nucleic acid/lower case aminoacid/rna

"""

import argparse


def get_args():
    """Return parsed command-line arguments."""

    parser = argparse.ArgumentParser(
        description="Parsing the file names",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # get a list of text files to process
    parser.add_argument('file_list',  # variable to access this data later: args.file_list
                        metavar='FILE', # shorthand to represent the input value
                        help='Provide file name to process. For multiple files, separate their names with spaces.', # message to the user, it goes into the help menu
                        type=str,
                        # required=False,
                        nargs="+" # will combine multiple textfile inputs into a list
                        )

    return parser.parse_args()

def identify_sequence(sequence):
    """Identifying the sequence is aminoacid/rna/lower case nucleic acid/nucleic acid """

    # Create a flag for nucleic acid/aminoacid present/lowe-nucleic acid/rna sequence
    nucleic_acid_present = True
    amino_acid_present = True
    lower_nucleic_acid_present = True
    rna_sequence_present = True
    """Identify a DNA sequence in lowercase"""

    #  Standard nucleotides
    standard_nucleic_acid = ["A", "T", "G", "C"]
    standard_lower_case_nucleic_acid = ["a", "t", "g", "c"]
    standard_rna_sequence = ["A", "U", "G", "C"]

    # Standard/Essential Amino acid abbreviation
    # source: https://www.technologynetworks.com/applied-sciences/articles/essential-amino-acids-chart-abbreviations-and
    # -structure-324357
    standard_amino_acid = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W",
                           "Y",
                           "V", "O", "U", "B", "Z", "X"]

    # Identify the Nucleic acid sequence/ Aminoacid sequence/rna  in the file
    for value in sequence:
        if value not in standard_nucleic_acid:
            nucleic_acid_present = False
        if value not in standard_amino_acid:
            amino_acid_present = False
        if value.lower() not in standard_lower_case_nucleic_acid:
            lower_nucleic_acid_present = False
        if value not in standard_rna_sequence:
            rna_sequence_present = False

    sequence_type = ""

    if nucleic_acid_present:
        sequence_type = "nucleic acid"
    elif lower_nucleic_acid_present:
        sequence_type = "lower nucleic acid"
    elif rna_sequence_present:
        sequence_type = "rna sequence"
    elif amino_acid_present:
        sequence_type = "amino acid"

    return sequence_type


if __name__ == "__main__":
    # get the arguments
    args = get_args()
    print(args.file_list)

    for file_path in args.file_list:
        # Initialize a line counter
        line_count = 0

        # Initialize sequence value
        seq = ''

        data_path = "../data/"

        # Open the genome file
        with open(data_path + file_path) as sequence:
            # Iterate over the lines in the file
            for line in sequence:
                # Increment line count by one
                line_count += 1
                # Appending sequence line to string
                seq += line

        #  Identifying the length of the Raw sequence
        print(f"The length of the sequence is : {len(seq)} ")

        # Removing the new line character in the Raw sequence
        cleaned_data = seq.replace("\r", "").replace("\n", "")

        sequence_type = identify_sequence(cleaned_data)

        print(f"Filename is {file_path} and sequence type is {sequence_type}")