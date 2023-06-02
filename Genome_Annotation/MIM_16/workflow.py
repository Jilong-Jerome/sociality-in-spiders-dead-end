from gwf import *
from workflow_targets import *
gwf = Workflow()

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/blat"
db = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/repeat_masker/MIM_hifi_hic_scaffolded_trim.fa.masked.16.part_02"
query = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/BI_mrna.fa"
outname = "MIM16_part2_blat"
run_blat(gwf,path,db,query,outname)

for num_id in ["001","002","003","004","005","006","007",
               "008","009","010","011","012","013","014"]:
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/blat"
    db = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/repeat_masker/MIM_hifi_hic_scaffolded_trim.fa.masked.16.part_02"
    query = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/BI_mrna.fa.split/BI_mrna.part_{num}.fa".format(num=num_id)
    outname = "MIM16_part2_blat_{num_id}".format(num_id=num_id)
    run_blat(gwf,path,db,query,outname)
    run_filter_blat(gwf,path,outname)

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/MIM16_part2"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/repeat_masker/MIM_hifi_hic_scaffolded_trim.fa.masked.16.part_02"
cfg = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/augustus/extrinsic.bug.cfg"
hints = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM16/star_merged/MIM16_augustus.hints.gff"
outname = "augu_MIM16_part2"
run_augustus(gwf,path,fasta,cfg,hints,outname)
