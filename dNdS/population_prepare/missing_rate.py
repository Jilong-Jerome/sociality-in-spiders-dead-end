import sys
from Bio import SeqIO
import matplotlib.pyplot as plt

def calculate_missing_fraction(sequence):
    total_length = len(sequence)
    missing_count = sequence.upper().count('N')
    return missing_count / total_length

def main(fasta_file):
    gene_missing_fractions = []
    high_missing_rate_genes = []

    with open(fasta_file, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            missing_fraction = calculate_missing_fraction(str(record.seq))
            gene_missing_fractions.append(missing_fraction)
            if missing_fraction > 0.1:
                high_missing_rate_genes.append((record.id, missing_fraction))

    # Draw distribution plot
    plt.hist(gene_missing_fractions, bins=50)
    plt.xlabel('Missing Fraction')
    plt.ylabel('Frequency')
    plt.title('Distribution of Missing Nucleotide Fractions')
    plt.savefig('missing_fraction_distribution.png')

    # Print gene names with missing rate > 10%
    print("\nGene names with a missing rate higher than 10%:")
    for gene_name, missing_rate in high_missing_rate_genes:
        print(f"{gene_name}: {missing_rate * 100:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python missing_rate.py <fasta_file>")
    else:
        fasta_file = sys.argv[1]
        main(fasta_file)

