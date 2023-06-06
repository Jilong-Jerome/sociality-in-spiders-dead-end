# sociality-in-spiders-dead-end
Pipelines and scripts for data analysis of the project "Sociality in spiders is an evolutionary dead end"

## [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md)
Please see details in [Genome Assembly](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/README.md) page. 
### Hifiasm
We use [hifiasm](https://www.nature.com/articles/s41592-020-01056-5) version 0.16.1-r375 to assemble inital halpotype-resolved contigs with Pacbio HiFi long reads and the HiC reads
### 3D-DNA scaffolding
We use [Juicer](https://www.sciencedirect.com/science/article/pii/S2405471216302198?via%3Dihub), [Juicebox](https://www.sciencedirect.com/science/article/pii/S240547121500054X?via%3Dihub), and [3D-DNA](https://github.com/aidenlab/3d-dna) pipeline to align the HiC reads to the assembled contigs and futher scafflod contigs into chromosome-level scaffolds according to the HiC contact pattern.

## [Genome Annotation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Annotation/README.md)
Please see details in [Genome Annotation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Annotation/README.md) page. 
### RepeatModeling and Repeatmasking
We use [RepeatModeler2](https://www.pnas.org/doi/10.1073/pnas.1921046117) to generate repeat library for each species. Species-specific repeat library is combined with the Arthropoda repeat library from [Repbase](https://www.girinst.org/repbase/) for [RepeatMasker](https://www.repeatmasker.org) to mask the reference genome of each specices.
### RNA-seq alignment
We use [STAR](https://github.com/alexdobin/STAR) to align RNA-seq to each species reference genome as hints for gene prediction from transcriptome.
### BRAKER
We ran [BRAKER2](https://github.com/Gaius-Augustus/BRAKER) pipeline to do gene prediction using hints from RNA transcriptome and protein homology simultaneously for each chromosome of each species separtely. Then the results from chromosmes are combined into species annotations.

Special cases were handled for HiC_scaffold_11 of *S.dumicola* and half of the HiC_scaffold_16 of *S.mimosarum*, where we failed to ran BRAKER2 ETP mode thourgh them.

The final species annotation completeness is checked with [BUSCO](https://busco.ezlab.org) together with a set of 1013 Arthropoda genes from OrthoDB v10.
### Gene Synteny
We use [GENESPACE](https://github.com/jtlovell/GENESPACE) to plot the gene synteny across our *de novo* assemblies. During the GENESPACE analysis, [OrthoFinder](https://github.com/davidemms/OrthoFinder) was performed for finding ortholougous groups among the six species: *S.lineatus*, *S.mimosarum*, *S.dumicola*, *S.tentoriicola*, *S.sarasinorum* and *S.bicolor*. 

## [dN/dS estimation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/README.md)
Please see details in [dN/dS estimation](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/README.md) page. 

### Orthogroups filtering
After GENESPACE/OrthoFinder, we obtained the single-copy orthologs across the 6 species where we have chromosome level assemblies. To integrate the RNA-seq from *S.africanus* and short DNA-seq from *S.pacificus*, we apply extra steps to filter and prepare high-quality codon-alignment across all the 8 species.

#### *S.africanus*
For *S.africanus*, where we have only RNA-seq, we use Trinity to assemble transcripts from the RNA-seq. The assembled transcripts are aligned to all the single-copy orthologs found in the 6 species with chromosome-level assemblies. Every single-copy orthologs that obtained a hit in the assembled transcripts with an anverage simialrity over 75% is kept.

#### *S.pacificus*
For *S.pacificus*, we align the short DNA-seq to the reference genome of *S.sarasinorum*, which is the closest species to *S.pacificus*. We create a consensus reference sequence for *S.pacificus* based on the alignment. The genes belonging to a certain single-copy ortholog groups are retreived based on genome annotation of *S.sarasinorum*

### Alignment
We use [MACSEv2](https://academic.oup.com/mbe/article/35/10/2582/5079334) to align the retrive coding region of single-copy orthologs across all 8 species. The alignment is further filtered for continious size and fraction of polymorphisim of local alignment blocks to avoid artefact from local mis-alignments. The final alignments of each single-copy orthologs is converted into phylip format with BioPython for further analysis.

### PAML branch-wise dN/dS
We random select 500, 100 genes out of the 2302 autosomal genes, 347 X chromosome genes for 500 times. In each bootstrapping, the alignments of the selected genes were concatenated as the input for CodeML in [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) to do branch-wise dN/dS estimation.

### PAML pairwise dN/dS
We start with mapping reads from each population to the species reference genome and calling SNPs with [Playtypus](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data). The reference genome of each population is created by substituing the reference genome with alternative nucleotdie.
 
We random select 500, 100 genes out of the 2302 autosomal genes, 347 X chromosome genes for 500 times. In each bootstrapping, sequence of the coding region from the selected genes were retrived and concatenated from the two population reference genome as the input for CodeML in PAML to do pair-wise dN/dS estimation.

### Solving for social transition time

### Species divergence time

## PSMC


## Codes for plotting
