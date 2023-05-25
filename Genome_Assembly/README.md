# Genome assembly

## Hifiasm
We used Pacbio HiFi long reads and HiC sequence for each species as the input for Hifiasm.

The output of Hifiasm is two assembly graphs in GFA format, each represent a haplotype.

We use awk to retrive the fasta sequence from the GFA file of the longer haplotype
## 3D-DNA scaffolding
