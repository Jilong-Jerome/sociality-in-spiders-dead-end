# Genome annotation

## Repeatmasking
1. We used RepeatModeler2 to model a species-specfic repeat liberay from the scaffolded chromosome-level assembly.

The workflow in gwf for analysising all species is in the attached [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/hifiasm/workflow.py)

Showcase using *S.sarasinorum* as an example
```
#Specifying path to species genome and results output 
sp = "SARA"
out = sp
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)

# Run RepeatModeler
BuildDatabase -name {out} -engine ncbi {genome}
RepeatModeler -pa 24 -engine ncbi -database {out}
```
2. The arthropoda repeat library from Repbase and the modeled species-specific library are combined as the repeat library input for RepeatMasker to do soft-masking.

Showcase using *S.sarasinorum* as an example
```
#Prepare inputs
sp = "SARA"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
lib1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/repbase/repbase_arthropoda.fa"
lib2 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/model/{sp}-families.fa".format(sp=sp)
lib = "{sp}_repbase.fa".format(sp=sp)

#Run RepeatMasker (softmasking)
cat {lib1} {lib2} > {lib}
RepeatMasker -e ncbi -pa 24 -xsmall -dir {path} -lib {lib} {fasta}
```
## Gathering evidence for gene prediction

### Evidence from RNA-seq

1. We use STAR to align RNA-seq from single individual to the softmasked species reference genome.

Showcase of the workflow to align individuals of a single species.
```

# indexing genome using STAR
STAR --runThreadN 18 --runMode genomeGenerate --genomeDir {path} --genomeFastaFiles {ref}

# Setting inputs
for ind in sp_dict: 
    fq1 = sp_dict[ind][0] 
    fq2 = sp_dict[ind][1] 
    indname = ind.replace("-","_") 

# Run STAR mapping
STAR --genomeDir {index_path} --runThreadN 16 --readFilesCommand zcat --readFilesIn {fq1} {fq2} --outSAMattrRGline ID:{indname} --outFileNamePrefix {align_path}/{indname}/{indname} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
```
2. We use samtools to merge and sort the individal bam files into one single bam file for each species. 

Showcase of the workflow
```
# Merge bams
/home/jilong/software/samtools-1.12/samtools merge -@ 16 -O BAM {path}/{outname}.bam {files}
# Sort reads
/home/jilong/software/samtools-1.12/samtools sort -@ 16 -O BAM -o {bam_out}.bam {bam_in}.bam
# Index bam
/home/jilong/software/samtools-1.12/samtools index -@ 16 -b {bam_out}.bam {bam_out}.bai
```
### Evidence from protein homology
1. For protein homology evidence, we used protein sequence from a previous [S.dumicola annotation from NCBI](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Stegodyphus_dumicola/100/) and [arthropoda protein from OrthoDB v10](https://academic.oup.com/nar/article/47/D1/D807/5160989)

## BRAKER

1. We split the chromosome-level assemblies into single chromosome to be ran through the BRAKER pipeline.

Showcase of the chromsome splitting
```
seqkit grep -p HiC_scaffold_{chrom_id} {genome} -o {out}
```

2. We ran BRAKER2 pipe line on each single chromsomes, with the RNA alignment bam file and the comined protein sequence fasta as hints for gene prediction.

Showcase of running BRAKER2 in ETP mode (taking transcriptome and protein homology evidence at the same time).
```
# Run BRAKER2 for each chromosomes
braker.pl --species={aug_sp} --genome={genome} --prot_seq={protein} --bam {bam} --etpmode --softmasking --cores=12 --gff3

# Combine BRAKER results for species annotations
cat {string_gff} > {sp}_braker.gff3
```
### Special Cases
#### HiC_scaffold_11 of *S.dumicola*
We run BRAKER2 with only uning hints from transcriptome on the whole genome of *S.dumicola*, then we retreive only the gene predicted with full support from transcriptome data.
Example codes
```
# Run BRAKER RNA mode
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine/DUM_hifi_hic_scaffolded_trim.fa.masked"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam"
braker.pl --species={species} --genome={genome} --bam={bam} --softmasking on --cores=24 --gff3

#Check the supports of gene hints from transcriptome
gtf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/augustus.gtf"
hint = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/hintsfile.gff"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/support"
/home/jilong/software/BRAKER/scripts/predictionAnalysis/selectSupportedSubsets.py {gtf} {hint} --fullSupport full --anySupport any --noSupport no
```
