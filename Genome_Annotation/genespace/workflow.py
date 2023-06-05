from gwf import *
from workflow_targets import *
gwf = Workflow()
for ste in chrom_group_dict:
    if ste not in ["Ste_1"]:
        ste_gff_logs = []
        ste_pep_logs = []
        for sp in chrom_group_dict[ste]:
            gff_logs = []
            pep_logs = []
            gff_list = []
            pep_list = []
            for chrom in chrom_group_dict[ste][sp]:
                gff_log = run_reform_gff(gwf,ste,sp,chrom)
                pep_log = run_reform_pep(gwf,ste,sp,chrom)
                gff_logs.append(gff_log)
                pep_logs.append(pep_log)
                gff_list.append("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/{ste}/rawGenomes/{sp}/{sp}/annotation/{sp}_{chrom}_gene.gff".format(ste=ste,sp=sp,chrom=chrom))
                pep_list.append("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/{ste}/rawGenomes/{sp}/{sp}/annotation/{sp}_{chrom}_pep.fa".format(ste=ste,sp=sp,chrom=chrom))
            gwf.target_from_template(
                name = "gff_format_{sp}_{ste}".format(sp=sp,ste=ste),
                template = logs_sum(gff_logs,LOG_PATH+"/gff_{sp}_{ste}.DONE".format(sp=sp,ste=ste))
                )
            gwf.target_from_template(
                name = "pep_format_{sp}_{ste}".format(sp=sp,ste=ste),
                template = logs_sum(pep_logs,LOG_PATH+"/pep_{sp}_{ste}.DONE".format(sp=sp,ste=ste))
                )
            ste_gff_log = run_cat_gff(gwf,ste,sp,gff_list)
            ste_pep_log = run_cat_pep(gwf,ste,sp,pep_list)
            ste_gff_logs.append(ste_gff_log)
            ste_pep_logs.append(ste_pep_log)
        gwf.target_from_template(
            name = "gff_format_{ste}".format(ste=ste),
            template = logs_sum(ste_gff_logs,LOG_PATH+"/gff_{ste}.DONE".format(ste=ste))
            )
        gwf.target_from_template(
            name = "pep_format_{ste}".format(ste=ste),
            template = logs_sum(ste_pep_logs,LOG_PATH+"/pep_{ste}.DONE".format(ste=ste))
            )
        run_genespace_ste(gwf,ste)
        run_draw_genespace_ste(gwf,ste)

palindrome_line = open("parlindrome_region.tsv")
for line in palindrome_line:
    sp,cluster,logs = run_find_fasta_cluster(gwf,line)
    run_combine_fasta_cluster(gwf,sp,cluster,logs)
    run_diamond_build_db(gwf,sp,cluster)
