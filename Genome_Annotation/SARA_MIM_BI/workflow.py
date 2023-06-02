from gwf import *
from workflow_dicts import *
from workflow_templates import *
from workflow_targets import *
gwf = Workflow()
## STAR index
for sp in ["SARA","MIM","BI","MIM16"]:
    build_STAR_index_sp(sp,gwf)
## STAR align
logs_dict ={}
for sp in ["SARA","MIM","BI","MIM16"]:
    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_align".format(sp=sp)
    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_index".format(sp=sp)
    logs = STAR_align_dict(sp_fq_dict[sp],sp,align_path,index_path,gwf)
    gwf.target_from_template(
        name = "STAR_align_{sp}".format(sp=sp),
        template = logs_sum(logs,LOG_PATH+"/align_STAR_{sp}.DONE".format(sp=sp))
        )
    logs_dict[sp]=logs
for sp in ["SARA","MIM","BI","MIM16"]:
    path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged".format(sp=sp)
    filepath = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_align".format(sp=sp)
    filenames = sp_fq_dict[sp].keys() ##change points if need to swtich species
    log = merge_STAR_bam(path,sp,filepath,filenames,logs_dict[sp],gwf)
    path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged".format(sp=sp)
    bam_in = "STAR_merged_{sp}".format(sp=sp)
    bam_out = "STAR_merged_{sp}_sort".format(sp=sp)
    log = sort_STAR_bam(path,bam_in,bam_out,sp,gwf) 

## repeat_mask genome
for sp in ["SARA","MIM","BI"]:
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/model".format(sp=sp)
    out = sp
    model_sp(genome,sp,path,out,gwf)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/mask".format(sp=sp)
    lib1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/repbase/repbase_arthropoda.fa"
    lib2 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/model/{sp}-families.fa".format(sp=sp)
    lib = "{sp}_repbase.fa".format(sp=sp)
    repeat_sp(genome,lib1,lib2,lib,path,sp,gwf)

## split masked genome by chrom
sp = "SARA"
numbers = 19
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
sp = "MIM"
numbers = 19
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
sp = "BI"
numbers = 18
logs = grep_chrom_fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "split_masked_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/Split_masked_chrom_{sp}.DONE".format(sp=sp))
    )
## braker split masked genome
sp = "SARA"
numbers = 17
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
sp = "MIM"
numbers = 19
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
lists = ["Ste_3","16.part_01","16.75per","16.round"]
logs = braker_sp_chrom_list(sp,lists,gwf)
### Try for MIM scaffold 16 part 2
sp="MIM16"
chrom_id="16_part_1"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/braker/HiC_scaffold_16_part_1"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/repeat_masker/MIM_hifi_hic_scaffolded_trim.fa.masked.16.part_01"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/STAR_merged_MIM16_sort.bam"
protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
task_name = "braker_MIM16_part_1"
aug_sp = "MIM_braker_MIM16_part_1"
gwf.target_from_template(
    name=task_name,
    template = braker_combine_mode(path,sp,aug_sp,chrom_id,genome,protein,bam)
)

sp="MIM16"
chrom_id="16_part_2"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/braker/HiC_scaffold_16_part_2"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/repeat_masker/MIM_hifi_hic_scaffolded_trim.fa.masked.16.part_02"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/STAR_merged_MIM16_sort.bam"
protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.formated.fa"
task_name = "braker_MIM16_part_2"
aug_sp = "MIM_braker_MIM16_part_2"
gwf.target_from_template(
    name=task_name,
    template = braker_protein_mode(path,sp,aug_sp,chrom_id,genome,protein,bam)
)

sp = "BI"
numbers = 17
logs = braker_sp_chrom(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/BRAKER_{sp}.DONE".format(sp=sp))
    )
## gff2fa
sp = "SARA"
numbers = chrom_id_dict[sp]
logs = agat_gff2fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_gff2fa_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/gff2fa_{sp}.DONE".format(sp=sp))
    )
sp = "MIM"
numbers = chrom_id_dict[sp]
logs = agat_gff2fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_gff2fa_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/gff2fa_{sp}.DONE".format(sp=sp))
    )
sp = "BI"
numbers = chrom_id_dict[sp]
logs = agat_gff2fa(sp,numbers,gwf)
gwf.target_from_template(
    name = "braker_gff2fa_{sp}".format(sp=sp),
    template = logs_sum(logs,LOG_PATH+"/gff2fa_{sp}.DONE".format(sp=sp))
    )
## cat all results from braker
sp = "SARA"
numbers = chrom_id_dict[sp]
combine_results_braker_fa(sp,numbers,gwf)
## BUSCO check annotation completeness
sp = "SARA"
busco_sp(sp,gwf)
for sp in ["DUM","TENT","MIM","SARA","BI","LIN"]:
    busco_per_sp(sp,gwf)
