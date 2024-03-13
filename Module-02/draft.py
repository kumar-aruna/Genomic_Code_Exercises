def scoreMotif(seq, motif, base_idx):
    """Scores a 13bp sequence based on the motif matrix."""
    score = 0
    for i, base in enumerate(seq):
        score += motif[base_idx[base]][i]
    return score

def find_orf(sequence, start_index):
    """Find ORF starting from a start codon (ATG or GTG) to one of the stop codons (TAA, TAG, TGA) or end of sequence."""
    stop_codons = ['TAA', 'TAG', 'TGA']
    for i in range(start_index, len(sequence) - 2, 3):  # Iterate over sequence in steps of 3 (codon size)
        codon = sequence[i:i+3]
        if codon in stop_codons:
            orf_length = i - start_index
            if orf_length >= 60:
                return sequence[start_index:i], i  # Return ORF and position of stop codon
    # If no stop codon found but length criteria is met, return the sequence as a potential gene
    if len(sequence) - start_index >= 60:
        return sequence[start_index:], len(sequence)
    return None, None


def process_contig_sequence(sequence, motif, base_idx, header, output_file):
    """Process a single contig sequence to find ORFs in high-scoring 13bp segments, considering ATG and GTG as start codons."""
    for i in range(len(sequence) - 12):  # Go through each 13bp segment
        segment = sequence[i:i+13]
        score = scoreMotif(segment, motif, base_idx)
        if score > 7.25:
            # Search for start codon in the high-scoring segment
            for j in range(i, len(sequence) - 2):
                if sequence[j:j+3] == 'ATG' or sequence[j:j+3] == 'GTG':
                    # Find ORF starting from this start codon
                    orf, stop_pos = find_orf(sequence, j)
                    if orf:
                        orf_length = stop_pos - j
                        output_file.write(f">{header}|Length {orf_length}|at position {j}\n{orf}\n")
                        break  # Stop searching after finding a valid ORF for this start codon

def process_contigs_file(input_file, motif, base_idx, output_file_path):
    with open(input_file, 'r') as file, open(output_file_path, 'w') as output_file:
        header = ''
        sequence = ''
        for line in file:
            if line.startswith('>'):
                if sequence:  # If there's an accumulated sequence, process it
                    process_contig_sequence(sequence, motif, base_idx, header, output_file)
                    sequence = ''  # Reset sequence for next contig
                header = line.strip().split('|')[0][1:]  # Update header to contain only contig number
            else:
                sequence += line.strip()
        # Don't forget to process the last contig in the file
        if sequence:
            process_contig_sequence(sequence, motif, base_idx, header, output_file)

# Define the base index and motif matrix
base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
motif = [
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # A
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],      # T
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],    # C
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]    # G
]

# Example usage (
import random

l = 150  # read length
k = 100
bases = ["A", "T", "C", "G"]


def generateReads(l, k):  # prints k reads of length l

    for i in range(k):
        seq = ""
        for j in range(l):
            seq += bases[random.randint(0, 3)]
        print(">sequence " + str(i + 1))
        print(seq)


generateReads(l, k)

s = "ATG ATC ATA CGT CCG ATC GTT CAA TGC GTG TAG CCC CGG GAC CGT GAG TGT GGT ACT AGT AGC CGT GTC AT"

print(len(s))