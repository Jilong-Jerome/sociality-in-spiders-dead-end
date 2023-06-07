# dN/dS estimation

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

2. We used Diamond to build sequence database of all the protein sequence retreived and translated from the single-copy orthologus of all the six species. Then we use blastx option to query every transcipt sequence from the *S.africanus* trinity transcriptome.
See complete codes in [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/diamond/workflow.py) 


```
# Prepare inputs
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/single_ortho_sequences/single_copy_pep.fa"
out = "stegodyphus_single_copy"

# Make diamond database
diamond makedb --in {fasta} -d {out}

# Prepare inputs
query = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/trinity/africanus/africanums_trinity.Trinity.fasta"
db = "stegodyphus_single_copy"
out = "africanus_hit"
n = 1

# Run diamond blastx
diamond blastx -p 24 -k {n} -q {query} -d {db} -o {out}.tsv --ultra-sensitive
```
3. After the Diamond blast, we matched the transcripts from *S.africanus* to the single copy orthologs found across the other six species. Length of the matched *S.africanus* transcripts and average transcript length of the orthologs is checked. For each *S.africanus* transcript, the longest alignment with any orthologs is kept and the alignment length should be over 200 amino acid. The sequence identity between *S.africanus* and other ortholog sequences are evaluated and a minmum threshold of 0.75 is applied for any matched *S.africanus* transcript to be valid.
See complete codes in [og_filter.Rmd](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/AFR/og_filter.Rmd)


