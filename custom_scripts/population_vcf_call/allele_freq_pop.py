import pysam
import sys

def calculate_allele_counts(population_ids, genotypes_all, ref_allele, alt_allele):
    ref_count = 0
    alt_count = 0
    missing_count = 0

    for idx in population_ids:
        genotype = genotypes_all[idx]
        if genotype is None:  # Missing genotype
            missing_count += 1
        else:
            # Map nucleotides to reference/alternative alleles
            for allele in genotype:
                if allele == ref_allele:
                    ref_count += 1
                elif allele == alt_allele:
                    alt_count += 1
                else:  # Handle unknown alleles
                    missing_count += 1

    total_alleles = ref_count + alt_count
    alt_freq = alt_count / total_alleles if total_alleles > 0 else 0
    return ref_count, alt_count, missing_count, alt_freq

def process_vcf(vcf_file, individual_a, individual_b, output_file):
    vcf = pysam.VariantFile(vcf_file)
    samples = list(vcf.header.samples)

    # Separate populations by ID prefixes
    population_a_ids = [i for i, s in enumerate(samples) if s.startswith("a")]
    population_b_ids = [i for i, s in enumerate(samples) if s.startswith("b")]
    all_population_ids = population_a_ids + population_b_ids

    if individual_a not in samples or individual_b not in samples:
        print("Error: Specified individuals not found in VCF.")
        sys.exit(1)

    individual_a_idx = samples.index(individual_a)
    individual_b_idx = samples.index(individual_b)

    with open(output_file, 'w') as out:
        out.write(
            "Chrom\tPos\tRef\tAlt\tGenotype_A\tGenotype_B\t"
            "A_Ref_Count\tA_Alt_Count\tA_Missing_Count\tA_Alt_Freq\t"
            "B_Ref_Count\tB_Alt_Count\tB_Missing_Count\tB_Alt_Freq\t"
            "Total_Ref_Count\tTotal_Alt_Count\tTotal_Missing_Count\tTotal_Alt_Freq\n"
        )

        for record in vcf.fetch():
            ref_allele = record.ref
            alt_allele = record.alts[0]

            # Get genotypes for all individuals
            genotypes_all = [record.samples[sample].alleles for sample in samples]
            genotype_a = genotypes_all[individual_a_idx]
            genotype_b = genotypes_all[individual_b_idx]

            # Skip if either selected individual has missing alleles
            if genotype_a is None or genotype_b is None or None in genotype_a or None in genotype_b:
                continue

            # Handle missing genotypes by replacing None with "."
            genotype_a_str = "/".join([allele if allele is not None else "." for allele in genotype_a])
            genotype_b_str = "/".join([allele if allele is not None else "." for allele in genotype_b])

            # Check if the site is polymorphic between the two individuals
            alleles = set(genotype_a + genotype_b)
            if ref_allele in alleles and alt_allele in alleles:
                # Calculate allele counts and frequencies for populations A, B, and both
                a_ref, a_alt, a_missing, a_alt_freq = calculate_allele_counts(
                    population_a_ids, genotypes_all, ref_allele, alt_allele
                )
                b_ref, b_alt, b_missing, b_alt_freq = calculate_allele_counts(
                    population_b_ids, genotypes_all, ref_allele, alt_allele
                )
                total_ref, total_alt, total_missing, total_alt_freq = calculate_allele_counts(
                    all_population_ids, genotypes_all, ref_allele, alt_allele
                )

                # Write results to file
                out.write(f"{record.chrom}\t{record.pos}\t{ref_allele}\t{alt_allele}\t")
                out.write(f"{genotype_a_str}\t{genotype_b_str}\t")
                out.write(f"{a_ref}\t{a_alt}\t{a_missing}\t{a_alt_freq:.4f}\t")
                out.write(f"{b_ref}\t{b_alt}\t{b_missing}\t{b_alt_freq:.4f}\t")
                out.write(f"{total_ref}\t{total_alt}\t{total_missing}\t{total_alt_freq:.4f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python filter_vcf.py <vcf_file> <individual_a> <individual_b> <output_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    individual_a = sys.argv[2]
    individual_b = sys.argv[3]
    output_file = sys.argv[4]

    process_vcf(vcf_file, individual_a, individual_b, output_file)

