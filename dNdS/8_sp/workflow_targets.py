from workflow_templates import *
from gwf import *
def get_modified_sara_value(file_path):
    with open(file_path, "r") as file:
        line = file.readline()
        values = line.strip('\n').split('\t')

        for value in values:
            if value.startswith("SARA"):
                modified_value = value + ".t1"
                return modified_value
def build_group_fasta(gwf,path,group):
    group_list = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/7_dnds/single_ortho/{group}.list".format(group=group)
    pac_name = get_modified_sara_value(group_list)
    pac_cds = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/cds_fasta/PAC_128_concensus_cds.fa"
    gwf.target_from_template(
            name = "retrieve_{group}_fasta".format(group=group),
            template = retrieve_ortholog_fasta(path,pac_cds,pac_name,group)
        )
    return LOG_PATH+"/ortho_fasta_{group}.DONE".format(group=group)

def run_group_fasta(group,gwf,path,AFR):
    infos = build_group_dict(group)
    mrna_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/annotations/cds_fasta"
    gwf.target_from_template(
            name = "retrieve_{group}_fasta".format(group=group),
            template = retrieve_single_ortho_fa(group,infos,mrna_path,path,AFR)
        )
    return LOG_PATH+"/ortho_fasta_{group}.DONE".format(group=group)

def run_group_align(group,path,gwf):
    gwf.target_from_template(
            name = "ortho_align_{group}".format(group=group),
            template = align_ortho(path,group)
        )
    return LOG_PATH+"/ortho_align_{group}.DONE".format(group=group)
def run_group_paml(mode,group,gwf):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/7_dnds/codeml/{mode}/{group}".format(mode=mode,group=group)
    gwf.target_from_template(
            name = "paml_{mode}_{group}".format(mode=mode,group=group),
            template = paml(path,mode,group)
        )
    return LOG_PATH+"/paml/paml_{mode}_{group}.DONE".format(mode=mode,group=group)
    
