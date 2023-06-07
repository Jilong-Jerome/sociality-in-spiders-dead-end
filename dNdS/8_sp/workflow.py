from workflow_targets import *
from workflow_templates import *

gwf = Workflow()
groups = "all_og_uniq.txt"
#groups = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/sp_7_OG.demo"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/single_ortho"
group_id = open(groups)
for group in group_id:
    build_group_fasta(gwf,path,group.strip("\n"))
    run_group_align(group.strip("\n"),path,gwf) 
