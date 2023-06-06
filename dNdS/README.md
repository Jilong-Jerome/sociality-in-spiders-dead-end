# dN/dS estimation

## Orthogroups filtering
1. We used RepeatModeler2 to model a species-specfic repeat liberay from the scaffolded chromosome-level assembly.

The workflow in gwf for analysising all species is in the attached [workflow.py](https://github.com/Jilong-Jerome/sociality-in-spiders-dead-end/blob/main/Genome_Annotation/repeat_masking/workflow.py)

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
