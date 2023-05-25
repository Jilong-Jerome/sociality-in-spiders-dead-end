# Genome assembly

## Hifiasm
1. We used Pacbio HiFi long reads and HiC sequence for each species as the input for Hifiasm.

The output of Hifiasm is two assembly graphs in GFA format, each represent a haplotype.

The workflow in gwf for analysising all species is in the attached [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/hifiasm/workflow.py)

Showcase using *S.dumicola* as an example
```
#Specifying path to HiFi reads, HiC reads and results output 
hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/DUM/DUM_HiFi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/DUM/DUM_hifi"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_2.fq"

# Run hifiasm
hifiasm -o {hifi_out} -t32 --h1 {hic_1} --h2 {hic_2} {hifi_reads}

```
2. We use awk to retrive the fasta sequence from the GFA file of the longer haplotype resloved.

Showcase using *S.dumicola* as an example
```
awk '/^S/{print ">"$2"\n"$3}' /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/DUM/DUM_hifi.tmp.hic.hap2.p_ctg.gfa > /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/DUM/DUM_hifi.tmp.hic.hap2.p_ctg.fa
```
## 3D-DNA scaffolding

The workflow in gwf for analysising all species is in the attached [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Assembly/hic_scaffold/workflow.py)

1. First we index the contig fasta files with bwa and use the [Juicer](https://github.com/aidenlab/juicer) pipeline to aligned paired sequenced HiC reads to the indexed contigs. 

Showcase of the setting parameters for Juicer aligning process
```
#Specifiying species name and hifiasm assembled contigs
species = "DUM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg.fa"
#Run Juicer pipeline
/home/jilong/software/juicer/scripts/juicer.sh -d /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/{species} -D /home/jilong/software/juicer -p /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/{species}/chrom.sizes -s none -z {fasta} -q short -Q 12:00:00 -l normal -L 24:00:00 -t 36 > {species}_juicer.log
```
