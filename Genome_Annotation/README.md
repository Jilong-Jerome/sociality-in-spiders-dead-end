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
