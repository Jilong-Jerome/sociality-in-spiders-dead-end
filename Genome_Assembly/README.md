# Genome assembly

## Hifiasm
We used Pacbio HiFi long reads and HiC sequence for each species as the input for Hifiasm.

The output of Hifiasm is two assembly graphs in GFA format, each represent a haplotype.

The workflow for all species is attached [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/workflow.py)

Showcase using S.dumicola as an example
```
#Specifying path to HiFi reads, HiC reads and results output 
hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/DUM/DUM_HiFi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/DUM/DUM_hifi"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_2.fq"

# Run hifiasm
hifiasm -o {hifi_out} -t32 --h1 {hic_1} --h2 {hic_2} {hifi_reads}

```
We use awk to retrive the fasta sequence from the GFA file of the longer haplotype
## 3D-DNA scaffolding
