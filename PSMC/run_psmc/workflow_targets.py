from workflow_templates import *
from workflow_dicts import *

BAM_PATH="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/"
PSMC_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/PSMC"
def run_psmc_fq_ind(gwf,sp,ind,chrom):
    path =  PSMC_PATH+"/{sp}/{ind}".format(sp=sp,ind=ind)
    bam_in = BAM_PATH+"{sp}/{ind}_final.bam".format(sp=sp,ind=ind)
    fa_ref = sp2ref_dict[sp]
    gwf.target_from_template(
        name = "consensus_fq_{ind}_{chrom}".format(ind=ind,chrom=chrom),
        template = get_consensus_fq_chr(path,bam_in,chrom,fa_ref,ind)
    )
    return LOG_PATH+"/consensus_fq_{name}_{chrom}.DONE".format(name=ind,chrom=chrom)

def run_cat_list(gwf,file_list,logs,combine_name,path):
    gwf.target_from_template(
        name = "cat_{combine_name}".format(combine_name=combine_name),
        template = cat_files(file_list,logs,combine_name,path)
    )
    return LOG_PATH+"/cat_{combine_name}.DONE".format(combine_name=combine_name)

def run_fq2psmcfa(gwf,path,fq_file,combine_name):
    gwf.target_from_template(
        name = "psmcfa_{combine_name}".format(combine_name=combine_name),
        template = fq2psmcfa(path,fq_file,combine_name)
    )
    return LOG_PATH+"/fq2psmcfa_{combine_name}.DONE".format(combine_name=combine_name)

def run_psmc(gwf,path,combine_name,psmcfa,psmc_out):
    gwf.target_from_template(
        name = "psmc_{combine_name}".format(combine_name=combine_name),
        template = psmc(path,psmcfa,psmc_out)
    )
    return LOG_PATH+"/psmc_{combine_name}.DONE".format(combine_name=psmc)
