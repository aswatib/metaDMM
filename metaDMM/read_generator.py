import os
import random
import numpy as np
from Bio import SeqIO
import pandas as pd
from scipy.special import gammaln, digamma

def generate_random_reads(genome_sequence, read_length, num_reads):
    """
    Generate random reads from a genome sequence using a fixed read length.

    Args:
        genome_sequence (str): The input genome sequence.
        read_length (int): The fixed length of the generated reads.
        num_reads (int): The number of reads to generate.

    Returns:
        list: A list of the generated reads.
    """
    genome_length = len(genome_sequence)
    reads = []
    for _ in range(num_reads):
        start_position = random.randint(0, genome_length - read_length)
        read = genome_sequence[start_position:start_position + read_length]
        reads.append(read)
    return reads

def generate_sparse_dirichlet_counts(num_species, total_reads, sparsity_factor):
    """
    Generate counts using a sparse Dirichlet distribution.

    Args:
        num_species (int): The number of species or taxa.
        total_reads (int): The total number of reads to generate.
        sparsity_factor (float): The sparsity factor for the Dirichlet distribution.

    Returns:
        numpy.ndarray: The generated counts for each species.
    """
    proportions = np.random.dirichlet([sparsity_factor] * num_species)
    scaled_counts = np.round(np.random.multinomial(total_reads, proportions)).astype(int)
    return scaled_counts
