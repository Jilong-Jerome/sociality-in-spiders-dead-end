# PSMC

## Orthogroups filtering
After OrthoFinder2 and GENESPACE analysis, we already obtained the results of the coding region sequence of each orthologs in all six species with the chromosome-level assembly

1. Build transcriptome for *S.africanus*

See complete codes in [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/AFR/workflow.py) 

```
# Prepare inputs
read1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/intact_africanus_R1.fq"
read2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/intact_africanus_R2.fq"
outname = "africanums_trinity"

# Run Trinity
Trinity --seqType fq --trimmomatic --max_memory 250G --left {read1} --right {read2} --CPU 32 --output "{outname}"
```
