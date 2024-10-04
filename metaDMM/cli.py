import os
import click
import pandas as pd
import numpy as np
from Bio import SeqIO
from metaDMM.read_generator import generate_random_reads, generate_sparse_dirichlet_counts
from metaDMM.utils import calculate_alpha

@click.command()
@click.option('--genome-folder', '-g', required=True, type=click.Path(exists=True), help='Path to the folder containing genome FASTA files')
@click.option('--read-length', '-r', default=200, type=int, help='Fixed read length')
@click.option('--num-reads', '-n', default=1000, type=int, help='Number of reads per genome')
@click.option('--num-samples', '-s', default=20, type=int, help='Number of samples to generate')
@click.option('--sample-folder', '-o', required=True, type=click.Path(), help='Output folder for the generated samples')
@click.option('--stats-file', '-f', default='stats.csv', help='Filename for the sample statistics')
@click.option('--sparsity-factor', '-a', default=1.0, type=float, help='Sparsity factor for the Dirichlet distribution')
def cli(genome_folder, read_length, num_reads, num_samples, sample_folder, stats_file, sparsity_factor):
    try:
        # Initialize a dictionary to store genome sequences
        genomes = {}
        for genome_file in os.listdir(genome_folder):
            if not genome_file.endswith(('.fasta', '.fa', '.fna')):
                click.echo(f"Skipping non-FASTA file: {genome_file}")
                continue
            
            species_name = os.path.splitext(genome_file)[0]
            file_path = os.path.join(genome_folder, genome_file)
            
            try:
                with open(file_path, "r", encoding="latin1") as f:
                    records = list(SeqIO.parse(f, "fasta"))
                    if not records:
                        click.echo(f"Warning: No sequences found in {genome_file}")
                        continue
                    genome_sequence = str(records[0].seq)
                    genomes[species_name] = genome_sequence
                    click.echo(f"Successfully loaded genome: {species_name} (length: {len(genome_sequence)})")
            except Exception as e:
                click.echo(f"Error reading genome file {genome_file}: {str(e)}")

        if not genomes:
            click.echo("Error: No valid genome sequences were loaded. Please check your input files.")
            return

        click.echo(f"Number of genomes loaded: {len(genomes)}")

        # Initialize a list to store statistics for all samples
        all_samples_stats = []

        # Generate reads for each sample
        for sample_num in range(1, num_samples + 1):
            click.echo(f"Generating sample {sample_num}")
            sample_reads = {}
            sample_stats = {}

            # Generate the number of reads for each species using Sparse Dirichlet distribution
            num_reads_species = generate_sparse_dirichlet_counts(len(genomes), num_reads, sparsity_factor)
            click.echo(f"Number of reads per species: {num_reads_species}")

            # Iterate over species and generate reads
            for i, (species_name, genome_sequence) in enumerate(genomes.items()):
                num_reads = num_reads_species[i]
                reads = generate_random_reads(genome_sequence, read_length, num_reads)
                sample_reads[species_name] = reads
                sample_stats[species_name] = len(reads)
                click.echo(f"Generated {len(reads)} reads for {species_name}")

            # Save the statistics for the current sample
            all_samples_stats.append(sample_stats)

            # Create sample FASTQ file
            os.makedirs(sample_folder, exist_ok=True)
            sample_filename = os.path.join(sample_folder, f"Sample{sample_num}.fastq")
            with open(sample_filename, "w") as sample_file:
                for species_name, reads in sample_reads.items():
                    for i, read in enumerate(reads, start=1):
                        read_id = f"Sample{sample_num}_{species_name}_Read{i}"
                        fastq_entry = f"@{read_id}\n{read}\n+\n{'I'*len(read)}\n"
                        sample_file.write(fastq_entry)
            click.echo(f"Saved sample to {sample_filename}")

        # Convert sample statistics to a DataFrame for alpha calculation
        sample_stats_df = pd.DataFrame(all_samples_stats)
        stats_file_path = os.path.join(sample_folder, stats_file)
        sample_stats_df.to_csv(stats_file_path, index=False)
        click.echo(f"Saved statistics to {stats_file_path}")

        # Calculate alpha values for the generated reads
        data = sample_stats_df.values
        alpha_values, alpha_0 = calculate_alpha(data)

        # Output the calculated alpha values
        click.echo("Alpha values for the generated reads:")
        click.echo(alpha_values)
        click.echo(f"Sum of alpha values (alpha_0): {alpha_0}")
        click.echo(f"Mean of alpha values: {np.mean(alpha_values)}")

    except Exception as e:
        click.echo(f"An error occurred: {str(e)}")

if __name__ == '__main__':
    cli()