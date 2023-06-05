from workflow_templates import *
from workflow_dicts import *
def read_palindrome_tsv(line):
    info = line.strip("\n").split("\t")
    ste = info[0]
    sp = info[1]
    start = int(info[2])
    end = int(info[3])
    cluster = "{sp}_{ste}_{start}".format(sp=sp,ste=ste,start=start)
    return ste,sp,start,end,cluster
def num_2_OF(sp,num):
    OF = of_fasta_dict[sp]+"_"+str(num)
    return OF

def run_reform_gff(gwf,ste,sp,chrom):
    gff_in = sp_gff_path_dict[sp]+"/"+chrom+"/braker/braker.gff3"
    gwf.target_from_template(
        name = "reform_gff_{sp}_{ste}_{chrom}".format(sp=sp,ste=ste,chrom=chrom),
        template = format_gff(ste,sp,chrom,gff_in)
    )
    return LOG_PATH+"/{sp}_{ste}_{chrom}_gff_reform.DONE".format(sp=sp,ste=ste,chrom=chrom)
def run_reform_pep(gwf,ste,sp,chrom):
    if sp in ["MIM","SARA","BI"]:
        pep_in = sp_gff_path_dict[sp]+"/"+chrom+"/fasta/braker_protein_{sp}_{chrom}.fasta".format(sp=sp,chrom=chrom)
    else:
        pep_in = sp_gff_path_dict[sp]+"/{sp}_braker.protein.fasta".format(sp=sp)
    gwf.target_from_template(
        name = "reform_pep_{sp}_{ste}_{chrom}".format(sp=sp,ste=ste,chrom=chrom),
        template = format_pep_fa(ste,sp,chrom,pep_in)
    )
    return LOG_PATH+"/{sp}_{ste}_{chrom}_pep_fa_reform.DONE".format(sp=sp,ste=ste,chrom=chrom)

def run_cat_gff(gwf,ste,sp,gff_list):     
    gwf.target_from_template(
        name = "cat_gff_{sp}_{ste}".format(sp=sp,ste=ste),
        template = cat_group_gff(sp,ste,gff_list)
    )
    return LOG_PATH+"/cat_gff_{sp}_{ste}.DONE".format(sp=sp,ste=ste) 
def run_cat_pep(gwf,ste,sp,pep_list):     
    gwf.target_from_template(
        name = "cat_pep_{sp}_{ste}".format(sp=sp,ste=ste),
        template = cat_group_pep(sp,ste,pep_list)
    )
    return LOG_PATH+"/cat_pep_{sp}_{ste}.DONE".format(sp=sp,ste=ste)
def run_genespace_ste(gwf,ste):
    gwf.target_from_template(
        name = "genespace_{ste}".format(ste=ste),
        template = genespace_ste(ste)
    )
def run_draw_genespace_ste(gwf,ste):
    gwf.target_from_template(
        name = "draw_genespace_{ste}".format(ste=ste),
        template = draw_genespace_ste(ste)
    )
def run_find_fasta_cluster(gwf,line):
    ste,sp,start,end,cluster = read_palindrome_tsv(line)
    logs = [] 
    for i in range(start,end+1):
        OF = num_2_OF(sp,i)
        gwf.target_from_template(
            name = "find_fasta_{cluster}_{OF}".format(cluster=cluster,OF=OF),
            template = find_OF_fasta_sp(ste,OF,sp,cluster)
        )
        logs.append(LOG_PATH+"/getfa_{ste}_{sp}_{OF}.DONE".format(ste=ste,sp=sp,OF=OF))
    return sp,cluster,logs
def run_combine_fasta_cluster(gwf,sp,cluster,logs):
    gwf.target_from_template(
            name = "combine_fasta_{cluster}".format(cluster=cluster),
            template = combine_cluster_fasta(logs,sp,cluster)
        )
    return LOG_PATH+"/combine_fa_{cluster}.DONE".format(cluster=cluster)
def run_diamond_build_db(gwf,sp,cluster):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/region_seq/cross_compare/db"
    fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/region_seq/{sp}/{cluster}/{cluster}.fa".format(sp=sp,cluster=cluster)
    gwf.target_from_template(
            name = "diamond_db_{cluster}".format(cluster=cluster),
            template = diamond_build_db(path,fasta,cluster)
        )
    return LOG_PATH+"/diamond_DB_{cluster}.DONE".format(cluster=cluster) 
