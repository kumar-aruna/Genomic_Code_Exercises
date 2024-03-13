# # Define the motif scoring matrix and base index
# base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
# motif = [
#     [.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # A
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],      # T
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],    # C
#     [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]    # G
# ]
#
# def scoreMotif(seq, motif, base_idx):
#     """Scores a 13bp sequence based on the motif matrix."""
#     score = 0
#     for i, base in enumerate(seq):
#         if base in base_idx:  # Ensure base is recognized
#             score += motif[base_idx[base]][i]  # Add score for this base at position i
#     return score
#
# def process_contig_sequence(sequence, motif, base_idx, header, output_file):
#     """Process a single contig sequence to find and write 13bp segments with score > 7.25 to an output file."""
#     with open(output_file, 'a') as file:  # Open file in append mode
#         for i in range(len(sequence) - 12):  # Iterate through each 13bp segment
#             segment = sequence[i:i+13]
#             score = scoreMotif(segment, motif, base_idx)
#             if score > 7.25:
#                 print(segment)
#                 file.write(f"{header}: Position {i+1}, Segment {segment}, Score {score}\n")
#
# def process_contigs_file(input_file, output_file):
#     with open(input_file, 'r') as file:
#         header = ''
#         sequence = ''
#         for line in file:
#             if line.startswith('>'):  # Header line
#                 if sequence:  # If there's an accumulated sequence, process it
#                     process_contig_sequence(sequence, motif, base_idx, header, output_file)
#                 header = line.strip()  # Update header with current contig name
#                 sequence = ''  # Reset sequence for the next contig
#             else:
#                 sequence += line.strip()  # Accumulate sequence lines
#         # Process the last contig in the file
#         if sequence:
#             process_contig_sequence(sequence, motif, base_idx, header, output_file)
#
# # Specify the input and output file paths
# input_file_path = 'spaceSeq.fa'  # Update this path as needed
# output_file_path = 'scored_segments.txt'  # Output file for scored segments
#
# # Clear the output file before starting, to ensure it's fresh for each run
# with open(output_file_path, 'w') as file:
#     file.write("")  # Clear existing content
#
# # Process the contigs and write the results to the output file
# process_contigs_file(input_file_path, output_file_path)
#
# print("Processing complete. Scored segments are written to", output_file_path)
#


# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index
# Define the motif scoring matrix and base index


# Define the motif scoring matrix and base index
# Step 1: Define the scoring matrix and base index
base_idx = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
motif = [
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, 2, -99, -99, .5],  # A
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, 2, -99, 0],      # T
    [0, 0, 0, 0, 0, 0, 0, 0, 0, -99, -99, -99, 0],    # C
    [.5, .5, .5, .5, 0, 0, 0, 0, 0, .5, -99, 2, 0]    # G
]

# Step 2: Scoring function for a 13bp sequence
# Assume base_idx and motif are defined as previously stated

# Provided contig sequence for demonstration
contig_sequence = "AGCATTATTATGATCATACGTCCGATCGTTCAATGCGTGTAGCCCCGGGACCGTGAGTGTGGTACTAGTAGCCGTGTCAT"

# Function to score a sequence segment based on the motif matrix
def score_motif(seq, motif, base_idx):
    score = 0
    for i, base in enumerate(seq):
        if base in base_idx:
            score += motif[base_idx[base]][i]
    return score

# Function to find valid ORFs with no stop codon within the first 60 bases
def find_valid_orfs(sequence, base_idx, motif, min_length=60):
    stop_codons = ['TAA', 'TAG', 'TGA']
    valid_orfs = []

    # Iterate through the sequence to find "ATG" and ensure it's part of a scored motif
    for i in range(len(sequence) - 12):  # Adjust for motif length
        segment = sequence[i:i+13]
        score = score_motif(segment, motif, base_idx)
        if score > 7.25 and segment[9:12] == "ATG":
            print(segment)# Ensuring "ATG" at correct positions
            # Start scanning for ORFs from "ATG", ensuring the first 60 bases don't contain stop codons
            orf = ""
            stop_codon_found = False
            for j in range(i + 9, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in stop_codons:
                    if j - (i + 9) >= min_length:  # Check if ORF length before stop codon is at least min_length
                        valid_orfs.append(sequence[i+9:j+3])  # Include the sequence up to the stop codon
                    stop_codon_found = True
                    break
            if not stop_codon_found and len(sequence) - (i + 9) >= min_length:
                valid_orfs.append(sequence[i+9:])  # Include ORF if it extends to the end without stop codon

            break  # Process only the first valid ORF that meets criteria

    return valid_orfs

# Applying the function to the given sequence
valid_orfs = find_valid_orfs(contig_sequence, base_idx, motif)

# Printing the found valid ORFs
for orf in valid_orfs:
    print(f"Valid ORF found: {orf}")
    print(f"Length: {len(orf)}\n")


