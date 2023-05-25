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

2. We use the Juicer alignment and the contigs fasta as the input for 3D-DNA scaffolding pipeline. 

Showcase of the setting parameters for Juicer aligning process
```
# Specifiying inputs
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/3d_dna"
species = "DUM"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/aligned/merged_nodups.txt"
r = 0

#Run 3D_DNA assemblign process
mkdir -p {folder}
cd {folder}
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/tmp
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline.sh -r {r} --early-exit --editor-repeat-coverage 30 --editor-coarse-stringency 20 --splitter-input-size 500000 --splitter-coarse-resolution 500000 --splitter-coarse-stringency 20 {fasta} {merged} > {folder}/{species}_3ddna.log.tmp
mv {folder}/{species}_3ddna.log.tmp {folder}/{species}_3ddna.log
```

3. After 3D-DNA scaffolding process, we manually review and curate the contigs orders and orientation in the megascaffold from the output of 3D-DNA. After manual review in Juicebox, the megascaffold is split visually into chromsome-level scaffolds according to the HiC contact pattern.

Showcase of file name of scaffold assembly before and after the manual curation
```
#The megascaffold from 3D-DNA scaffolding, without rounds of misjoint correction.
DUM_hifi.tmp.hic.hap2.p_ctg.0.assembly
#After manual curation
DUM_hifi.tmp.hic.hap2.p_ctg.0.review.assembly
```

4. We use the manual reviewd assembly to finalize scafflod fasta sequence with 3D-DNA post-review pipeline.
Showcase of exporting the fasta sequence of each chromsome-level scaffold for *S.dumicola*

```
#Specifying inputs
review = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/3d_dna/DUM_hifi.tmp.hic.hap2.p_ctg.0.review.assembly"
draft = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/reference/DUM_hifi.tmp.hic.hap2.p_ctg_wrapped.fa"
merged = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/DUM/aligned/merged_nodups.txt"
folder = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM"
species = "DUM"
# Run 3D-DNA post-review for finalizing the reference output
cd {folder}
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/tmp
/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/3d_dna/3d-dna/run-asm-pipeline-post-review.sh --build-gapped-map --sort-output -s finalize -r {review} {draft} {merged} > {folder}/{species}_export.log.tmp
mv {folder}/{species}_export.log.tmp {folder}/{species}_export.log
```
