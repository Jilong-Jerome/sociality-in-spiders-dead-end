# sociality-in-spiders-dead-end
Pipelines and scripts for data analysis of the project "Sociality in spiders is an evolutionary dead end"

## [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md)
Please see details in [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md)
### Hifiasm
We use [hifiasm](https://www.nature.com/articles/s41592-020-01056-5) version 0.16.1-r375 to assemble inital halpotype-resolved contigs with Pacbio HiFi long reads and the HiC reads
### 3D-DNA scaffolding
We use [Juicer](https://www.sciencedirect.com/science/article/pii/S2405471216302198?via%3Dihub), [Juicebox](https://www.sciencedirect.com/science/article/pii/S240547121500054X?via%3Dihub), and [3D-DNA](https://github.com/aidenlab/3d-dna) pipeline to align the HiC reads to the assembled contigs and futher scafflod contigs into chromosome-level scaffolds according to the HiC contact pattern.

## Genome Annotation
### RepeatModeling and Repeatmasking
We use [RepeatModeler2](https://www.pnas.org/doi/10.1073/pnas.1921046117) to generate repeat library for each species. Species-specific repeat library is combined with the Arthropoda repeat library from Repbase for RepeatMasker to mask the reference genome of each specices.
### RNA-seq alignment
### BRAKER
### Gene Synteny

## dN/dS estimation
### Orthogroups filtering
### PAML branch-wise dN/dS
### PAML pairwise dN/dS

## PSMC

## Codes for plotting
