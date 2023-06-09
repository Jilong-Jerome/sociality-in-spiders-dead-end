# PSMC

## Creating artifitial diplod individual

For social species, we concatnated the fastq file of DNAseq from two individual into two fastq file R1.fastq and R2.fastq.

For subsocial species, R1 and R2 fastq from a single individual is used.

Then the concatened paired reads R1 and R2 are aligned to species reference genome in the same way as shown in dN/dS DNAseq alignment [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/dna_align/workflow.py)

## Prepare PSMC input

The complete workflow for running PSMC separately for X chromosomes and autosomes can be found at [workflow](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/PSMC/run_psmc/workflow.py)

1. Generate consensus fastq for each chromsome separately

```
#inputs
{ref} = species reference fasta
{chrom} = chromosome_id
{bam_in} = bam file for the individual
{name} = individual id

#Run samtools to obtain consensus fastq
/home/jilong/software/samtools-1.12/samtools mpileup -Q 30 -q 30 -u -v -f {ref} -r {chrom} {bam_in} |  bcftools call -c -V indels |  vcfutils.pl vcf2fq -d 5 -D 80 -Q 30 > {name}_{chrom}.fq
```

2. Combine consensus fastq from chromosmes into autosomes set and X chromosomes set

```
{file_list} = list of fastq file of chromsomes belonging to autosomes/X chromosomes

file_string = ""
    for filename in file_list:
        file_string = file_string + " {filename}".format(filename=filename)
cat {file_string} > {combine_name}
```

3. Convert consensus fastq to psmcfa format for psmc

```
fq2psmcfa {fq_file} > {combine_name}
```

## Run PSMC

```
psmc -p "4+25*2+4+6" -o {psmc_out} {psmcfa}
```

## Filling values and calculating Ne_X/Ne_A
The codes for parsing and combining results from autosomes and X chromsomes can be found at [fill_value.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/PSMC/fill_value.py)
```
#For autosomes
$PSMC_HOME/utils/psmc_plot.pl -R -u 5e-09 -g 1 -p {species_prefix} {species_prefix}.psmc
#For X chromosomes
$PSMC_HOME/utils/psmc_plot.pl -R -u 3.8e-09 -g 1 -p {species_prefix} {species_prefix}.psmc

#Parse outputs and calculating ratios
python fill_value.py {Autosome_results} {X_results} {ratios_output}
```
