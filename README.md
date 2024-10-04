# metaDMM

metaDMM is a tool for generating synthetic DNA reads using a Dirichlet Multinomial distribution (DMM for Dirichlet Multinomial Modeling). It simulates metagenomic data based on complete bacterial genomes with different concentration parameters (alpha).

## Installation

To install metaDMM, follow these steps:

1. Clone the repository:
   ```
   git clone https://github.com/aswatib/metaDMM.git
   ```

2. Navigate to the project directory:
   ```
   cd metaDMM
   ```

3. Install the package:
   ```
   pip install -e .
   ```

## Usage

To use metaDMM, run the following command:

```
metadmm --genome-folder /path/to/genomes --sample-folder /path/to/output --num-samples 5 --num-reads 1000 --read-length 150 --sparsity-factor 0.5
```

Arguments:
- `--genome-folder`: Path to the folder containing genome FASTA files
- `--sample-folder`: Output folder for the generated samples
- `--num-samples`: Number of samples to generate (default: 20)
- `--num-reads`: Number of reads per genome (default: 1000)
- `--read-length`: Fixed read length (default: 200)
- `--sparsity-factor`: Sparsity factor for the Dirichlet distribution (default: 1.0)

## Output

The tool generates:
1. FASTQ files for each sample in the specified output folder.
2. A `stats.csv` file containing sample statistics.

## License

This project is licensed under the APACHE 2.0 License (LICENSE).

## Acknowledgements

This tool was developed as the next step following the master's thesis on "Exploring Microbial Diversity in Metagenomics".

Development of this tool was assisted by Claude, an AI assistant created by Anthropic.

## Contact
Project Link: https://github.com/aswatib/metaDMM
