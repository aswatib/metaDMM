import os
import csv
import random
import numpy as np
from Bio import SeqIO
import pandas as pd
from scipy.special import gammaln, digamma

# Function to generate fastq reads from a genome sequence using fixed read length
def generate_random_reads(genome_sequence, read_length, num_reads):
    genome_length = len(genome_sequence)
    reads = []
    for _ in range(num_reads):
        start_position = random.randint(0, genome_length - read_length)
        read = genome_sequence[start_position:start_position + read_length]
        reads.append(read)
    return reads

# Function to generate counts using sparse Dirichlet distribution
def generate_sparse_dirichlet_counts(num_species, total_reads, sparsity_factor):
    proportions = np.random.dirichlet([sparsity_factor] * num_species)
    scaled_counts = np.round(np.random.multinomial(total_reads, proportions)).astype(int)
    return scaled_counts

# Function to calculate Dirichlet-Multinomial parameters using Method of Moments
def calculate_alpha(data):
    n = np.sum(data, axis=1)
    mean_proportions = np.mean(data / n[:, None], axis=0)
    variance_proportions = np.var(data / n[:, None], axis=0, ddof=1)
    alpha = mean_proportions * (((mean_proportions * (1 - mean_proportions)) / variance_proportions) - 1)
    alpha_0 = np.sum(alpha)
    return alpha, alpha_0

# Specify parameters
genome_folder = "entrez.genomes"
fixed_read_length = 200  # Fixed read length
fixed_num_reads_per_genome = 1000  # Fixed number of reads per genome
num_samples = 20
sample_folder = "Dirichlet_alpha1"
stats_file = "stats_Dirichlet_alpha1.csv"
sparsity_factor = 1  # Change this value to test different alphas

# Initialize a dictionary to store genome sequences
genomes = {}
for genome_file in os.listdir(genome_folder):
    species_name = os.path.splitext(genome_file)[0]
    with open(os.path.join(genome_folder, genome_file), "r", encoding="latin1") as f:
        genome_sequence = str(next(SeqIO.parse(f, "fasta")).seq)
        genomes[species_name] = genome_sequence

# Initialize a list to store statistics for all samples
all_samples_stats = []

# Generate reads for each sample
for sample_num in range(1, num_samples + 1):
    sample_reads = {}
    sample_stats = {}

    # Generate the number of reads for each species using Sparse Dirichlet distribution
    num_reads_species = generate_sparse_dirichlet_counts(len(genomes), fixed_num_reads_per_genome, sparsity_factor)

    # Iterate over species and generate reads
    for i, (species_name, genome_sequence) in enumerate(genomes.items()):
        num_reads = num_reads_species[i]
        reads = generate_random_reads(genome_sequence, fixed_read_length, num_reads)
        sample_reads[species_name] = reads
        sample_stats[species_name] = len(reads)

    # Save the statistics for the current sample
    all_samples_stats.append(sample_stats)

    # Create sample fastq file
    os.makedirs(sample_folder, exist_ok=True)
    sample_filename = os.path.join(sample_folder, f"Sample{sample_num}.fastq")
    with open(sample_filename, "w") as sample_file:
        for species_name, reads in sample_reads.items():
            for i, read in enumerate(reads, start=1):
                read_id = f"Sample{sample_num}_{species_name}_Read{i}"
                fastq_entry = f"@{read_id}\n{read}\n+\n{'I'*len(read)}\n"
                sample_file.write(fastq_entry)

# Convert sample statistics to a DataFrame for alpha calculation
sample_stats_df = pd.DataFrame(all_samples_stats)
sample_stats_df.to_csv(os.path.join(sample_folder, stats_file), index=False)

# Calculate alpha values for the generated reads
data = sample_stats_df.values
alpha_values, alpha_0 = calculate_alpha(data)

# Output the calculated alpha values
print("Alpha values for the generated reads:")
print(alpha_values)
print("Sum of alpha values (alpha_0):", alpha_0)
print("Mean of alpha values:", np.mean(alpha_values))

# Write all sample statistics to a CSV file
all_stats_filename = os.path.join(sample_folder, stats_file)
fieldnames = ['Sample'] + list(all_samples_stats[0].keys())
with open(all_stats_filename, "w", newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for idx, sample_stats in enumerate(all_samples_stats, start=1):
        sample_stats['Sample'] = f'Sample {idx}'
        writer.writerow(sample_stats)

print(f"All sample statistics saved to {all_stats_filename}")
