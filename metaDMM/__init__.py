"""
metaDMM: A tool for generating synthetic DNA reads using a Dirichlet Multinomial distribution.

This package provides a command-line interface and functions to generate synthetic DNA reads
from complete genome sequences, using a Dirichlet Multinomial distribution to simulate
the base composition and read length distribution.
"""

__all__ = [
    'generate_random_reads',
    'generate_sparse_dirichlet_counts',
    'calculate_alpha',
    'MetaDMMCLI',
]

from .read_generator import generate_random_reads, generate_sparse_dirichlet_counts
from .utils import calculate_alpha
from .cli import cli

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
