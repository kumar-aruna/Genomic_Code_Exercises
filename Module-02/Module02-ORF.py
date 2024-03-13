# Define the scoring matrix and base index
base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
motif = [
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # For A
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],      # For T
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],    # For C
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]    # For G
]

# Function to score a 13bp sequence based on the motif matrix
def score_motif(seq, motif, base_idx):
    score = 0
    for i, base in enumerate(seq):
        if base in base_idx:
            score += motif[base_idx[base]][i]
    return score

# Function to find valid ORFs in a sequence
def find_valid_orfs(sequence, base_idx, motif, min_length=60):
    stop_codons = ['TAA', 'TAG', 'TGA']
    valid_orfs = []

    for i in range(len(sequence) - 12):  # Adjust for motif length
        segment = sequence[i:i+13]
        if score_motif(segment, motif, base_idx) > 7.25:
            print(f"Motif: {segment} at position {i + 1} scored")  # Prin
            if segment[9:12] == "ATG":
                for j in range(i + 9, len(sequence), 3):
                    if sequence[j:j+3] in stop_codons:
                        if j - (i + 9) >= min_length:
                            valid_orfs.append((i + 10, sequence[i+9:j+3]))
                        break
                else:
                    if len(sequence) - (i + 9) >= min_length:
                        valid_orfs.append((i + 10, sequence[i+9:]))
                break  # Assuming only the first valid ORF is needed per contig

    return valid_orfs

# Function to process each contig sequence, find valid ORFs, and write to the output file
def process_sequence(sequence, header, outfile, base_idx, motif):
    valid_orfs = find_valid_orfs(sequence, base_idx, motif)
    for start_pos, orf in valid_orfs:
        contig_number = header.split()[0].lstrip('>')
        outfile.write(f"> {contig_number}|Length {len(orf)}|at position {start_pos}\n{orf}\n\n")

# Function to process all contigs from the input file and write results to the output file
def process_contigs(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence, header = '', ''
        for line in infile:
            if line.startswith('>'):
                if sequence:  # Process the previous sequence before reading the next header
                    process_sequence(sequence, header, outfile, base_idx, motif)
                header = line.strip()
                sequence = ''
            else:
                sequence += line.strip()
        if sequence:  # Don't forget to process the last sequence in the file
            process_sequence(sequence, header, outfile, base_idx, motif)

# Specify the path to your input file and the output file
input_file_path = 'spaceSeq.fa'  # Update this path to your actual input file
output_file_path = 'output_orfs.txt'  # The results will be written here

# Execute the processing of contigs to find and write ORFs
process_contigs(input_file_path, output_file_path)

print("Processing complete. Valid ORFs have been written to the output file.")
