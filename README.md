# sociality-in-spiders-dead-end
Pipelines and scripts for data analysis of the project "Sociality in spiders is an evolutionary dead end"

## [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md)
Please see details in [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md)
### Hifiasm
We use [hifiasm](https://www.nature.com/articles/s41592-020-01056-5) version 0.16.1-r375 to assemble inital halpotype-resolved contigs with Pacbio HiFi long reads and the HiC reads
### 3D-DNA scaffolding
We use [Juicer](https://www.sciencedirect.com/science/article/pii/S2405471216302198?via%3Dihub), [Juicebox](https://www.sciencedirect.com/science/article/pii/S240547121500054X?via%3Dihub), and [3D-DNA](https://github.com/aidenlab/3d-dna) pipeline to align the HiC reads to the assembled contigs and futher scafflod contigs into chromosome-level scaffolds according to the HiC contact pattern.

## [Genome Annotation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Annotation/README.md)
Please see details in [Genome Annotation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Annotation/README.md)
### RepeatModeling and Repeatmasking
We use [RepeatModeler2](https://www.pnas.org/doi/10.1073/pnas.1921046117) to generate repeat library for each species. Species-specific repeat library is combined with the Arthropoda repeat library from Repbase for RepeatMasker to mask the reference genome of each specices.
### RNA-seq alignment
We use [STAR](https://github.com/alexdobin/STAR) to align RNA-seq to each species reference genome as hints for gene prediction from transcriptome.
### BRAKER
We ran [BRAKER2](https://github.com/Gaius-Augustus/BRAKER) pipeline to do gene prediction using hints from RNA transcriptome and protein homology simultaneously for each chromosome of each species separtely. Then the results from chromosmes are combined into species annotations.

Special cases were handled for HiC_scaffold_11 of *S.dumicola* and half of the HiC_scaffold_16 of *S.mimosarum*, where we failed to ran BRAKER2 ETP mode thourgh them.

### Gene Synteny

## dN/dS estimation
### Orthogroups filtering
### PAML branch-wise dN/dS
### PAML pairwise dN/dS

## PSMC

## Codes for plotting
