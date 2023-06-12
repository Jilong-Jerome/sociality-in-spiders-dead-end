from workflow_targets import *
from gwf import *
gwf = Workflow()
for sp in ["BI","SARA","DUM","TENT","MIM","LIN","PAC"]:
    for ind in sp2ind_dict[sp]:
        logs = []
        for chrom in sp2chrom_dict[sp]:
            log=run_psmc_fq_ind(gwf,sp,ind,chrom)
            logs.append(log)
        gwf.target_from_template(
            name = "consensus_fq_{sp}_{ind}".format(sp=sp,ind=ind),
            template = logs_sum(logs,LOG_PATH+"/con_fq_{sp}_{ind}.DONE".format(sp=sp,ind=ind))
        )
        logs = []
        files = []
        path =  PSMC_PATH+"/{sp}/{ind}".format(sp=sp,ind=ind)
        for chrom in sp2Xchrom_dict[sp]:
            logs.append(LOG_PATH+"/consensus_fq_{name}_{chrom}.DONE".format(name=ind,chrom=chrom))
            files.append(path+"/{ind}_{chrom}.fq".format(ind=ind,chrom=chrom))
        combine_name = "consensus_fq_{ind}_X.fq".format(ind=ind)
        run_cat_list(gwf,files,logs,combine_name,path)
        run_fq2psmcfa(gwf,path,combine_name,combine_name.replace(".fq",".psmcfa"))
        combine_name = "consensus_fq_{ind}_X".format(ind=ind)
        psmcfa_in = "consensus_fq_{ind}_X.psmcfa".format(ind=ind)
        psmc_out = "consensus_fq_{ind}_X.psmc".format(ind=ind)
        run_psmc(gwf,path,combine_name,psmcfa_in,psmc_out)
        logs = []
        files = []
        for chrom in sp2Achrom_dict[sp]:
            logs.append(LOG_PATH+"/consensus_fq_{name}_{chrom}.DONE".format(name=ind,chrom=chrom))
            files.append(path+"/{ind}_{chrom}.fq".format(ind=ind,chrom=chrom))
        combine_name = "consensus_fq_{ind}_A.fq".format(ind=ind)
        run_cat_list(gwf,files,logs,combine_name,path)
        run_fq2psmcfa(gwf,path,combine_name,combine_name.replace(".fq",".psmcfa"))
        combine_name = "consensus_fq_{ind}_A".format(ind=ind)
        psmcfa_in = "consensus_fq_{ind}_A.psmcfa".format(ind=ind)
        psmc_out = "consensus_fq_{ind}_A.psmc".format(ind=ind)
        run_psmc(gwf,path,combine_name,psmcfa_in,psmc_out)
        if sp == "PAC" and ind == "PAC_127":
            for chrom in sp2Achrom_dict[sp]:
                combine_name = "consensus_fq_{ind}_{chrom}.fq".format(ind = ind,chrom=chrom)
                logs = [LOG_PATH+"/consensus_fq_{name}_{chrom}.DONE".format(name=ind,chrom=chrom)]
                files = [path+"/{ind}_{chrom}.fq".format(ind=ind,chrom=chrom)]
                run_cat_list(gwf,files,logs,combine_name,path)
                run_fq2psmcfa(gwf,path,combine_name,combine_name.replace(".fq",".psmcfa"))
                combine_name = "consensus_fq_{ind}_{chrom}".format(ind=ind,chrom=chrom)
                psmcfa_in = "consensus_fq_{ind}_{chrom}.psmcfa".format(ind=ind,chrom=chrom)
                psmc_out = "consensus_fq_{ind}_{chrom}.psmc".format(ind=ind,chrom=chrom)
                run_psmc(gwf,path,combine_name,psmcfa_in,psmc_out)
