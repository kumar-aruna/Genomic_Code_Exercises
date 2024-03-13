"""
This function takes Seqread file as an input file, find Maximum overlaps
 and generates output file with Contigs.
"""


def readReads(reads_file):
    """This function Reads DNA sequences from the input seqread_file
    and returns them as a list."""
    with open(reads_file, 'r') as file:
        # Read each line, strip to remove  whitespace,
        # and return a list of non-empty lines
        return [line.strip() for line in file.readlines() if line.strip()]


def findMaxOverlap(contig, reads, k):
    """
    Find the read with the maximum overlap with the contig.
    :param contig: The current contig sequence being built using reads.
    :param reads: A list of remaining reads to check for maximum overlap.
    :param k: Minimum length of overlap to consider and generate contigs.
    :return: The read with the maximum overlap and the length of that overlap happened.
    """
    max_overlap = 0  # Initialize the maximum overlap length as Zero
    max_read = None  # Initialize the read with the maximum overlap

    # Iterate through each read in the list
    for read in reads:
        # Find the maximum possible overlap
        overlap = max(len(contig), len(read))

        # Check for overlap from k to maximum possible overlap to make a contig
        for i in range(k, overlap):
            # Check if the current contig ends or starts with a substring of the read
            if contig.endswith(read[:i]) or contig.startswith(read[-i:]):
                # Update max_overlap and max_read if a larger overlap is found
                if i > max_overlap:
                    max_overlap, max_read = i, read

    return max_read, max_overlap


def assembleReads(reads, k):
    """
    Assembles reads into contigs based on overlapping sequences.
    :param reads: A list of DNA sequence reads.
    :param k: Minimum length of overlap required to merge reads.
    :return: A list of assembled contigs.
    """
    contigs = []  # Initialize the empty list of contigs

    # Iterate until all reads are processed
    while reads:
        # Start a new contig with the first read in the list
        contig = reads.pop(0)

        # Continuously try to find and merge reads with maximum overlap
        while True:
            max_read, max_overlap = findMaxOverlap(contig, reads, k)
            if max_read:
                # Merge the contig with the read with maximum overlap
                if contig.endswith(max_read[:max_overlap]):
                    contig += max_read[max_overlap:]
                else:
                    contig = max_read[:-max_overlap] + contig

                # Remove the read that was merged into the contig
                reads.remove(max_read)
            else:
                # No more overlaps found, break the loop
                break

        # Add the final contig to the list of contig
        contigs.append(contig)

    return contigs


def writeContigs(contigs, output_filename):
    """Writes the list of contigs to an output file in FASTA format."""
    # Open the output file for writing
    with open(output_filename, 'w') as file:
        # Write each contig to the file in FASTA format
        for i, contig in enumerate(contigs, 1):
            file.write(f">Contig_{i}\n{contig}\n")



# Input file path where reads are located
input_seq_file_path = 'seqReadFile2023.txt'
# Output file path where contigs are located
output_file_path = 'contigs.fasta'
k = 10  # Minimum overlap length for assembling reads

# Read the DNA sequences from the file
reads = readReads(input_seq_file_path)
# Assemble the reads into contigs
contigs = assembleReads(reads, k)
# Write the assembled contigs to the output file
writeContigs(contigs, output_file_path)
