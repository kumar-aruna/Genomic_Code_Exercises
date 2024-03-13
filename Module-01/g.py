# The code provided earlier calculates coverage but not density.
# Density for a contig can be calculated as the number of reads covering each base divided by the total number of reads.

# First, let's modify the calculate_coverage function to calculate density instead.

def calculate_density(contigs, reads):
    # Calculate coverage first
    coverage = {contig: [0] * len(contig) for contig in contigs}
    for read in reads:
        for contig in contigs:
            index = contig.find(read)
            while index != -1:
                for i in range(index, index + len(read)):
                    coverage[contig][i] += 1
                index = contig.find(read, index + 1)

    # Now calculate density based on the coverage
    density = {contig: [0] * len(contig) for contig in contigs}
    for contig in contigs:
        total_reads_covering = sum(coverage[contig])
        for i in range(len(contig)):
            if total_reads_covering > 0:  # To avoid division by zero
                density[contig][i] = coverage[contig][i] / total_reads_covering
    return density


# Next, we need to modify the plot_coverage function to plot density instead of coverage.

def plot_density(density):
    plt.figure(figsize=(15, 5))
    # Plot each contig's density on the same graph for comparison
    for contig, dens in density.items():
        plt.plot(dens, label=contig[:10] + '...')  # Use the first 10 characters of the contig as the label
    plt.legend()
    plt.xlabel('Position in contig')
    plt.ylabel('Density')
    plt.title('Read Density across all contigs')
    plt.show()


# Now let's use the modified functions with your file paths
contigs = read_fasta('contigs.fasta')  # Replace with your actual contigs file path
reads = read_reads('')  # Replace with your actual reads file path

# Calculate density
density = calculate_density(contigs, reads)

# Plot the density
plot_density(density)

