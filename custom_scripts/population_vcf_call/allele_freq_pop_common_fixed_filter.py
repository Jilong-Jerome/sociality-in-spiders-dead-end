import pysam
import sys
from collections import defaultdict


def get_sample_indices(header_line, population_samples):
    """
    Get the indices of the samples for a specific population from the VCF header line.
    """
    header_parts = header_line.strip().split("\t")
    sample_indices = {sample: i for i, sample in enumerate(header_parts[9:], start=9)}
    return [sample_indices[sample] for sample in population_samples]


def write_population_vcf(vcf_header, output_file, population_samples, records):
    """
    Write a VCF file for a specific population.
    Includes only the given samples in the header and their genotypes in the records.
    """
    with open(output_file, 'w') as out:
        header_lines = str(vcf_header).splitlines()
        sample_indices = []

        # Write the header, filtering the samples
        for line in header_lines:
            if line.startswith("#CHROM"):
                # Modify the sample list in the header
                header_parts = line.split("\t")
                sample_indices = get_sample_indices(line, population_samples)
                filtered_header = "\t".join(header_parts[:9] + [header_parts[i] for i in sample_indices])
                out.write(filtered_header + "\n")
            else:
                out.write(line + "\n")

        # Write the records
        for record in records:
            fields = record.strip().split("\t")  # Use strip() to remove trailing newline
            genotypes = fields[9:]  # All genotype fields
            filtered_genotypes = [genotypes[i - 9] for i in sample_indices]
            out.write("\t".join(fields[:9] + filtered_genotypes) + "\n")


def write_common_summary(common_summary_file, chromosome_summary, thresholds):
    """
    Write the summary file for common variants.
    """
    with open(common_summary_file, 'w') as out:
        out.write("Chromosome\tThreshold\tTotal_SNPs\tCommon_Variants\tFraction_Common\n")
        total_snps_genome = 0
        total_common_by_mac = defaultdict(int)

        for chrom, counts in chromosome_summary.items():
            total_snps_genome += counts["total_snps"]
            for threshold in thresholds:
                common_variants = counts["common_variants_by_mac"][threshold]
                total_common_by_mac[threshold] += common_variants
                fraction_common = common_variants / counts["total_snps"] if counts["total_snps"] > 0 else 0
                out.write(f"{chrom}\t{threshold}\t{counts['total_snps']}\t{common_variants}\t{fraction_common:.4f}\n")

        for threshold in thresholds:
            genome_common = total_common_by_mac[threshold]
            genome_fraction_common = genome_common / total_snps_genome if total_snps_genome > 0 else 0
            out.write(f"Genome\t{threshold}\t{total_snps_genome}\t{genome_common}\t{genome_fraction_common:.4f}\n")


def write_fixed_summary(fixed_summary_file, chromosome_summary):
    """
    Write the summary file for fixed SNPs.
    """
    with open(fixed_summary_file, 'w') as out:
        out.write("Chromosome\tFixed_Differently_SNPs\tTotal_SNPs\tFraction_Fixed_Differently\n")
        total_fixed_genome = 0
        total_snps_genome = 0

        for chrom, counts in chromosome_summary.items():
            chrom_fixed = counts["fixed_differently"]
            total_fixed_genome += chrom_fixed
            total_snps_genome += counts["total_snps"]
            fraction_fixed = chrom_fixed / counts["total_snps"] if counts["total_snps"] > 0 else 0
            out.write(f"{chrom}\t{chrom_fixed}\t{counts['total_snps']}\t{fraction_fixed:.4f}\n")

        genome_fraction_fixed = total_fixed_genome / total_snps_genome if total_snps_genome > 0 else 0
        out.write(f"Genome\t{total_fixed_genome}\t{total_snps_genome}\t{genome_fraction_fixed:.4f}\n")


def process_vcf(
    vcf_file,
    individual_a,
    individual_b,
    output_file,
    common_summary_file,
    fixed_summary_file,
    filtered_vcf_file,
    pop_a_vcf_file,
    pop_b_vcf_file,
    mac_upper_limit
):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)

    # Separate populations by ID prefixes
    population_a_ids = [i for i, s in enumerate(samples) if s.startswith("a")]
    population_b_ids = [i for i, s in enumerate(samples) if s.startswith("b")]
    population_a_samples = [samples[i] for i in population_a_ids]
    population_b_samples = [samples[i] for i in population_b_ids]

    if individual_a not in samples or individual_b not in samples:
        print("Error: Specified individuals not found in VCF.")
        sys.exit(1)

    individual_a_idx = samples.index(individual_a)
    individual_b_idx = samples.index(individual_b)

    # Store population-specific records
    pop_a_records = []
    pop_b_records = []

    # Create filtered VCF file
    with open(filtered_vcf_file, 'w') as filtered_vcf:
        filtered_vcf.write(str(vcf.header))

        # Initialize counters for common variants and fixation
        chromosome_summary = defaultdict(lambda: {
            "total_snps": 0,
            "fixed_differently": 0,
            "common_variants_by_mac": defaultdict(int)
        })

        thresholds = list(range(2, mac_upper_limit + 1))

        with open(output_file, 'w') as out:
            out.write(
                "Chrom\tPos\tRef\tAlt\tGenotype_A\tGenotype_B\t"
                "A_Ref_Count\tA_Alt_Count\tA_Missing_Count\tA_Alt_Freq\t"
                "B_Ref_Count\tB_Alt_Count\tB_Missing_Count\tB_Alt_Freq\t"
                "Total_Ref_Count\tTotal_Alt_Count\tMinor_Allele_Count\tFixed_Differently\n"
            )

            for record in vcf.fetch():
                ref_allele = record.ref
                alt_allele = record.alts[0]

                # Get genotypes for the two individuals
                genotype_a = record.samples[individual_a].alleles
                genotype_b = record.samples[individual_b].alleles

                # Skip if either individual has missing alleles
                if genotype_a is None or genotype_b is None or None in genotype_a or None in genotype_b:
                    continue

                # Filter for polymorphic sites between the two individuals
                if set(genotype_a) == set(genotype_b):
                    continue

                # Handle missing genotypes by replacing None with "."
                genotype_a_str = "/".join([allele if allele is not None else "." for allele in genotype_a])
                genotype_b_str = "/".join([allele if allele is not None else "." for allele in genotype_b])

                # Calculate allele counts for populations A and B
                genotypes_all = [record.samples[sample].alleles for sample in samples]
                a_ref = sum(1 for idx in population_a_ids for allele in genotypes_all[idx] if allele == ref_allele)
                a_alt = sum(1 for idx in population_a_ids for allele in genotypes_all[idx] if allele == alt_allele)
                a_missing = len(population_a_ids) * 2 - (a_ref + a_alt)

                b_ref = sum(1 for idx in population_b_ids for allele in genotypes_all[idx] if allele == ref_allele)
                b_alt = sum(1 for idx in population_b_ids for allele in genotypes_all[idx] if allele == alt_allele)
                b_missing = len(population_b_ids) * 2 - (b_ref + b_alt)

                # Determine minor allele count (MAC)
                mac = min(a_ref + b_ref, a_alt + b_alt)

                # Determine if ref and alt are fixed differently in populations A and B
                fixed_diff = (a_ref > 0 and a_alt == 0 and b_ref == 0 and b_alt > 0) or \
                             (a_ref == 0 and a_alt > 0 and b_ref > 0 and b_alt == 0)

                # Update chromosome-level summary
                chrom = record.chrom
                chromosome_summary[chrom]["total_snps"] += 1
                if fixed_diff:
                    chromosome_summary[chrom]["fixed_differently"] += 1
                    # Write to filtered VCF
                    filtered_vcf.write(str(record))

                # Write to population-specific record lists
                record_str = str(record)
                if a_ref == 0 and a_alt > 0:  # Alt fixed in population A
                    pop_a_records.append(record_str)
                if b_ref == 0 and b_alt > 0:  # Alt fixed in population B
                    pop_b_records.append(record_str)

                for threshold in thresholds:
                    if mac >= threshold:
                        chromosome_summary[chrom]["common_variants_by_mac"][threshold] += 1

                # Write detailed results to output.tsv
                out.write(f"{record.chrom}\t{record.pos}\t{ref_allele}\t{alt_allele}\t")
                out.write(f"{genotype_a_str}\t{genotype_b_str}\t")
                out.write(f"{a_ref}\t{a_alt}\t{a_missing}\t{a_alt / (a_ref + a_alt):.4f}\t")
                out.write(f"{b_ref}\t{b_alt}\t{b_missing}\t{b_alt / (b_ref + b_alt):.4f}\t")
                out.write(f"{a_ref + b_ref}\t{a_alt + b_alt}\t{mac}\t{fixed_diff}\n")

    # Write population-specific VCFs
    write_population_vcf(vcf.header, pop_a_vcf_file, population_a_samples, pop_a_records)
    write_population_vcf(vcf.header, pop_b_vcf_file, population_b_samples, pop_b_records)

    # Write summary files
    write_common_summary(common_summary_file, chromosome_summary, thresholds)
    write_fixed_summary(fixed_summary_file, chromosome_summary)


if __name__ == "__main__":
    if len(sys.argv) != 11:
        print("Usage: python filter_vcf.py <vcf_file> <individual_a> <individual_b> <output_file> <common_summary_file> <fixed_summary_file> <filtered_vcf_file> <pop_a_vcf_file> <pop_b_vcf_file> <mac_upper_limit>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    individual_a = sys.argv[2]
    individual_b = sys.argv[3]
    output_file = sys.argv[4]
    common_summary_file = sys.argv[5]
    fixed_summary_file = sys.argv[6]
    filtered_vcf_file = sys.argv[7]
    pop_a_vcf_file = sys.argv[8]
    pop_b_vcf_file = sys.argv[9]
    mac_upper_limit = int(sys.argv[10])

    process_vcf(
        vcf_file,
        individual_a,
        individual_b,
        output_file,
        common_summary_file,
        fixed_summary_file,
        filtered_vcf_file,
        pop_a_vcf_file,
        pop_b_vcf_file,
        mac_upper_limit
    )

