from workflow_targets import *
from workflow_templates import *

gwf = Workflow()

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs"
#MIM
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/MIM/MIM_hifi_hic_scaffolded_trim.fa"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/MIM/M_Antana_132_final.bam"
name = "M_Antana_132"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/MIM/M_Wee_146_final.bam"
name = "M_Wee_146"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/MIM/MIM_hifi_hic_scaffolded_trim.fa"
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/M_Antana_132_noindel.recode.vcf"
name = "M_Antana_132"
run_altref(gwf,path,vcfin,refin,name)
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/M_Wee_146_noindel.recode.vcf"
name = "M_Wee_146"
run_altref(gwf,path,vcfin,refin,name)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta"
gffin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/annotations/MIM_Ste_All_gene.gff3"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/M_Antana_132.fa"
name = "M_Antana_132"
run_agat(gwf,path,refin,gffin,name)
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/M_Wee_146.fa"
name = "M_Wee_146"
run_agat(gwf,path,refin,gffin,name)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cat_fasta"

fain = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta/M_Wee_146_cds_auto.fa"
faout = "M_Wee_146_cds_auto_cat.fa"
name = "M_Wee_146_auto"
run_concat(gwf,path,fain,faout,name)

fain = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta/M_Wee_146_cds_x.fa"
faout = "M_Wee_146_cds_x_cat.fa"
name = "M_Wee_146_x"
run_concat(gwf,path,fain,faout,name)

fain = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta/M_Antana_132_cds_auto.fa"
faout = "M_Antana_132_cds_auto_cat.fa"
name = "M_Antana_132_auto"
run_concat(gwf,path,fain,faout,name)

fain = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta/M_Antana_132_cds_x.fa"
faout = "M_Antana_132_cds_x_cat.fa"
name = "M_Antana_132_x"
run_concat(gwf,path,fain,faout,name)

#DUM
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/DUM_hifi_hic_scaffolded_trim.fa"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/DUM/D_karasburg_166_final.bam"
name = "D_karasburg_166"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/DUM/D_otawi_122_final.bam"
name = "D_otawi_122"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/DUM_hifi_hic_scaffolded_trim.fa"
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/D_karasburg_166_noindel.recode.vcf"
name = "D_karasburg_166"
run_altref(gwf,path,vcfin,refin,name)
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/D_otawi_122_noindel.recode.vcf"
name = "D_otawi_122"
run_altref(gwf,path,vcfin,refin,name)
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta"
gffin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/annotations/DUM_Ste_All_gene.gff3"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/D_karasburg_166.fa"
name = "D_karasburg_166"
run_agat(gwf,path,refin,gffin,name)
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/D_otawi_122.fa"
name = "D_otawi_122"
run_agat(gwf,path,refin,gffin,name)

#SARA
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/SARA/SARA_Sri_final.bam"
name = "SARA_Sri"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/SARA/SARA_Hym_final.bam"
name = "SARA_Hym"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/SARA_Hym_noindel.recode.vcf"
name = "SARA_Hym"
run_altref(gwf,path,vcfin,refin,name)
vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/SARA_Sri_noindel.recode.vcf"
name = "SARA_Sri"
run_altref(gwf,path,vcfin,refin,name)
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta"
gffin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/annotations/SARA_Ste_All_gene.gff3"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/SARA_Hym.fa"
name = "SARA_Hym"
run_agat(gwf,path,refin,gffin,name)
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/SARA_Sri.fa"
name = "SARA_Sri"
run_agat(gwf,path,refin,gffin,name)

#refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs/PAC.fa"
refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/consensus/fasta/PAC_128_consensus.fa"
name = "PAC_128_concensus"
run_agat(gwf,path,refin,gffin,name)

#PAC
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/PAC/PAC_128_psmc_consensus_trim.fa"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/PAC/PAC_130_final.bam"
name = "PAC_130"
run_bam2vcf(gwf,path,bam,ref,name)
run_noindel(gwf,path,name)
#path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/alt_refs"
#refin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
#vcfin = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/vcfs/SARA_Hym_noindel.recode.vcf"
#name = "SARA_Hym"
#run_altref(gwf,path,vcfin,refin,name)
