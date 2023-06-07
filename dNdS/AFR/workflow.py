from gwf import *
from workflow_targets import *
gwf = Workflow()
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/trinity/africanus"
read1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/intact_africanus_R1.fq"
read2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/intact_africanus_R2.fq"
outname = "africanums_trinity"
run_trinity_single(path,read1,read2,outname,gwf)


