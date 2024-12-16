import pysam
import csv
import sys
from collections import defaultdict

# Define the input and output files
vcf_file = sys.argv[1]  # "input.vcf"
output_file = sys.argv[2]  # "output.tsv"
fraction_file = sys.argv[3]  # "fraction.tsv"

# Function to classify individuals into populations
def classify_individuals(samples):
    popA = [i for i, sample in enumerate(samples) if sample.startswith("a")]
    popB = [i for i, sample in enumerate(samples) if sample.startswith("b")]
    return popA, popB

# Function to count alleles in a population
def count_alleles_in_population(genotypes, pop_indices):
    ref_count = 0
    alt_count = 0
    na_count = 0

    for idx in pop_indices:
        genotype = genotypes[idx].split(":")[0]  # Extract genotype
        if genotype == "0/0":
            ref_count += 2
        elif genotype == "0/1" or genotype == "1/0":
            ref_count += 1
            alt_count += 1
        elif genotype == "1/1":
            alt_count += 2
        else:
            na_count += 2  # Missing or ambiguous genotypes

    return ref_count, alt_count, na_count

# Open the VCF file
vcf = pysam.VariantFile(vcf_file)
samples = list(vcf.header.samples)
popA_indices, popB_indices = classify_individuals(samples)

# Initialize counters for population-specific SNPs per chromosome
chromosome_data = defaultdict(lambda: {
    "popA_specific_snps": 0,
    "popB_specific_snps": 0,
    "shared_snps": 0,
    "total_snps": 0,
})

total_popA_specific_snps = 0
total_popB_specific_snps = 0
total_shared_snps = 0
total_snps = 0

# Open the output file for writing
with open(output_file, "w", newline="") as outfile:
    writer = csv.writer(outfile, delimiter="\t")

    # Write header
    writer.writerow([
        "chrom", "pos",
        "popA_ref_count", "popA_alt_count", "popA_na_count",
        "popB_ref_count", "popB_alt_count", "popB_na_count",
        "snp_specificity"
    ])

    # Process each variant in the VCF
    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        genotypes = [sample["GT"] for sample in record.samples.values()]

        # Handle missing genotype data
        genotypes = [
            "./." if g is None else "/".join(map(str, g)) 
            for g in genotypes
        ]

        # Count alleles for each population
        popA_counts = count_alleles_in_population(genotypes, popA_indices)
        popB_counts = count_alleles_in_population(genotypes, popB_indices)

        # Determine SNP specificity
        if popA_counts[1] > 0 and popB_counts[1] == 0:  # popA-specific
            snp_specificity = "popA_specific"
            chromosome_data[chrom]["popA_specific_snps"] += 1
            total_popA_specific_snps += 1
        elif popB_counts[1] > 0 and popA_counts[1] == 0:  # popB-specific
            snp_specificity = "popB_specific"
            chromosome_data[chrom]["popB_specific_snps"] += 1
            total_popB_specific_snps += 1
        elif popA_counts[1] > 0 and popB_counts[1] > 0:  # Shared SNP
            snp_specificity = "shared"
            chromosome_data[chrom]["shared_snps"] += 1
            total_shared_snps += 1
        else:  # No alternative alleles in either population
            continue

        # Increment total SNPs for the chromosome and genome
        chromosome_data[chrom]["total_snps"] += 1
        total_snps += 1

        # Write results to the output file
        writer.writerow([
            chrom, pos,
            *popA_counts,
            *popB_counts,
            snp_specificity
        ])

# Write the fraction file
with open(fraction_file, "w", newline="") as fracfile:
    writer = csv.writer(fracfile, delimiter="\t")

    # Write header
    writer.writerow([
        "chrom",
        "popA_specific_snps", "popB_specific_snps", "shared_snps", "total_snps",
        "fraction_popA_specific", "fraction_popB_specific", "fraction_shared", "fraction_total_population_specific"
    ])

    # Write data for each chromosome
    for chrom, data in chromosome_data.items():
        chrom_total_snps = data["total_snps"]
        popA_specific = data["popA_specific_snps"]
        popB_specific = data["popB_specific_snps"]
        shared_snps = data["shared_snps"]
        total_population_specific_snps = popA_specific + popB_specific
        fraction_popA = popA_specific / chrom_total_snps if chrom_total_snps > 0 else 0
        fraction_popB = popB_specific / chrom_total_snps if chrom_total_snps > 0 else 0
        fraction_shared = shared_snps / chrom_total_snps if chrom_total_snps > 0 else 0
        fraction_total_population_specific = total_population_specific_snps / chrom_total_snps if chrom_total_snps > 0 else 0
        writer.writerow([
            chrom, popA_specific, popB_specific, shared_snps, chrom_total_snps,
            f"{fraction_popA:.4f}", f"{fraction_popB:.4f}", f"{fraction_shared:.4f}", f"{fraction_total_population_specific:.4f}"
        ])

    # Write whole-genome level data
    genome_total_population_specific_snps = total_popA_specific_snps + total_popB_specific_snps
    genome_fraction_popA = total_popA_specific_snps / total_snps if total_snps > 0 else 0
    genome_fraction_popB = total_popB_specific_snps / total_snps if total_snps > 0 else 0
    genome_fraction_shared = total_shared_snps / total_snps if total_snps > 0 else 0
    genome_fraction_total_population_specific = genome_total_population_specific_snps / total_snps if total_snps > 0 else 0
    writer.writerow([
        "whole_genome",
        total_popA_specific_snps, total_popB_specific_snps, total_shared_snps, total_snps,
        f"{genome_fraction_popA:.4f}", f"{genome_fraction_popB:.4f}", f"{genome_fraction_shared:.4f}", f"{genome_fraction_total_population_specific:.4f}"
    ])

print(f"Output written to {output_file}")
print(f"Fraction report written to {fraction_file}")

