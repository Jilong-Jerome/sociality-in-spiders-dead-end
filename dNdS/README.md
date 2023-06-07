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

#### Length comparison
![Length of *S.africanus* transcripts and thier corresponding average length of orthologs in six species with chromosome-level assemblies](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/AFR/afr_og_length_comparison.jpeg)
#### Distribution of sequence identity among blast hits longer than 200 amino acids
![The distribution of sequence identity among blast hits longer than 200 amino acids](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/AFR/Distribution_of_sequence_identity_of_AFR_hits.jpeg)

4. We created a consensus sequence reference from a single *S.pacificus* individual. We aligned the paird short reads from *S.pacificus* to the reference of *S.sarasinorum* and then called a consensus sequence ignoring the indels to keep the consensus reference at the same length with *S.sarasinorum* so that we can use the same annotation to retrive gene sequence.

See complete codes in [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/PAC/workflow.py)

```
# Prepare inputs
id = "PAC_128"
bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/pacificus_consensus/pac128.bam"
fasta_out = id+"_consensus.fa"

# Run samtools consensus
samtools consensus -a --show-del yes --show-ins no -f fasta -o {fasta_out} {bam_in} -@ 8

```

## Alignments

### Multiple alignment of each single-copy orthologs across the 8 species in analysis

 According to the above filtering, we ending up with a selection of ortholog groups that has a corresponding coding sequence in each species of all the eight species included in the analysis. We retrieve the relevant fasta coding sequence and create a multiple alignment for each selected ortholog using MACSE.

See complete codes in [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/8_sp/workflow.py)

```
# Prepare inputs
group = #Orthologs ID#
# Run samtools consensus
macse -prog alignSequences -seq {group}_unalign.fasta -out_AA {group}_AA.fasta -out_NT {group}_NT.fasta -fs 10 -stop 10
```

### Concatenate alignments of a random set of genes with bootstrapping

See complete codes for alignment process and branch-wise dN/dS in [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/workflow.py)

We do bootstrapping estimation of dN/dS for X chromosomes and autosomes separately. We randomly sample 500 or 100 ortholog groups out of the 2302 autosomal genes or 347 X chromosome genes respectively. In the following example codes, we show how a single round of boostrapping being performed.

1. Select a random set of genes

#### Showcase of random selecting a set of 500 autosomal genes from the 2302 autosomal single-copy orthologs

```
# Prepare inputs
n = 500
out = "auto_{n}_{j}".format(n=n,j=i+1)
target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/auto_og_pass.txt"

# random select 500 genes from autosomal orthologs
shuf -n {n} {target} |cut -f1 >  {out}_id.txt 
```

2. Concatenate the alignments from the selected set of genes and then filter out regions that are prone to local mis-alignment. Continuous alignments without gaps are considerd as a single alignment block, only the alignment block longer than 300 nucleotide and polymorphysim fraction lower than 15% is kept.

Codes for alignment filtering can be found in [filter_region.R](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/filter_region.R) [find_region.R](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/find_region.R) [align_filter.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/align_filter.py), and [align_filter_region_fa.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/align_filter_region_fa.py)

#### Showcase of concatenating the alignments from concate

```
# Prepare inputs
cat_string  = ""
ids = open(path+"/"+id_list+"_id.txt")
for og_id in ids:
    phy_file = ALN_PATH+"/{og_id}_NT_fix.fasta".format(og_id=og_id.strip("\n"))
    cat_string = cat_string + phy_file + " "

# Concatenate alignments from the selected genes
conda activate goalign
goalign concat -i {cat_string} > {out}.cat

conda activate biopython
python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/align_filter.py {out}.cat {out}.region
echo raw_region_done

conda activate gwf
Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/find_region.R {out}.region {out}.region.all
echo region_summary_done

conda activate gwf
Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/filter_region.R {out}.region.all {out}.region.filtered
echo region_retrieve_done

conda activate biopython
python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/align_filter_region_fa.py {out}.cat {out}.region.filtered.tsv {out}.filtered.phy
echo retreive_region_align_done

conda activate goalign
goalign concat -i {out}.filtered.phy -p > {out}.filtered.concat.phy
conda activate clustalo
trimal -in {out}.filtered.concat.phy -phylip3.2 -out {out}.filtered.concat.paml.phy
```

## Branch-wise dN/dS estimation

The concatenated alignments after filtering are ready for estimating dN/dS in branch-wise for all the branches in the phylogeny of 8 *Stegodyphus* species.

The control file used for codeml brach-wise dN/dS can be found at [codeml.ctl](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/codeml_ofree.ctl)
The tree phylogeny used for estimation can be found at [tree.txt](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/boostrap/tree_ofree.txt)
#### Showcase of running codeml in PAML for brach-wise dN/dS
```
#Prepare inputs
out = "auto_{n}_{j}".format(n=n,j=i+1)
codeml_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/codeml"
mode = "ofree"
reform_script = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/reform_results.py"
alignment = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}/{out}.fa".format(out=out)

#Run CodeML and parse outputs
cp {codeml_path}/codeml_{mode}.ctl codeml.ctl
cp {codeml_path}/tree_{mode}.txt tree.txt
ln -f -s {group}.fa.filtered.concat.paml.phy seq.txt
conda activate paml
codeml     
echo {group} > {group}.name
paste {group}.name rst1 > {group}.res
rm {group}.name
tail -n 10 results.txt > results.raw
conda activate ete3
python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/combine_OG/boostrap/paml2tab.py results.raw {group}.tab {group}
python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/boostrap/pair_social.py {group}.tab {group}_social.tab
```


## Pair-wise dN/dS estimation
For each social species, we create two separate population reference using 1 single individual respectively. The pair-wise dN/dS estimation is performed with boostrapping, which resembles the branch-wise dN/dS process

1. Create reference genome for two separate populations.

We start with mapping a single individual from each population to the reference genome. We then call SNPs (ignoring indels) for each individual. The population reference is then created by substitute nucleotide at the SNP sites being the alternative nucleotide.

Codes for alignment of single individual paird short read DNAseq can be found at [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/dna_align/workflow.py)

Workflow for creating all population reference and retriveing coding sequence of orthologs can be found at [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/population_prepare/workflow.py)
#### Showcase of creating population reference
```
# Call variants
platypus callVariants --bamFiles={bam} --refFile={ref} --output={name}.vcf

# Remove indels
vcftools --vcf {vcfin} --remove-indels --recode --out {vcfout}

# Create new population reference
bgzip -c {vcfin} > {name}.noindel.vcf.gz
bcftools index {name}.noindel.vcf.gz
bcftools consensus -f {refin} --sample {name} {name}.noindel.vcf.gz > {refout}

# Retreive coding sequence from selected set of genes
agat_sp_extract_sequences.pl -g {gffin} -f {refin} -t cds -o {faout}.temp 
python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/pep_reformat_sp.py {sp} {faout}.temp {faout}
```

2. Boostraping pair-wise dN/dS estimation for random set of genes (autosomes only)

The full process resembles the same procedure as branch-wise dN/dS estimation, but with only a pair of sequence. Since the two population reference genome shares the same coordinate system, we did not filter for local mis-alignment as in previous multiple genome alignments.

The workflow for the population pair-wise dN/dS estimation can be found at [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/population_dnds/workflow.py)

The control file for codeml in PAML can be found at [codeml.ctl](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/dNdS/population_dnds/codeml.ctl)


