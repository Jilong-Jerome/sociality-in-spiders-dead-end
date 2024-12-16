from workflow_templates import *

def run_n_bootstrap_psmc(gwf,ROOT_PATH,sp,chrom_type,n,mu,g):
    data_path = ROOT_PATH+"/data/spiders"
    gwf.target_from_template(
    name = "split_psmcfa_{sp}_{chrom_type}".format(sp=sp,chrom_type=chrom_type),
    template = split_psmcfa(data_path,sp,chrom_type)
    )
    for i in range(n):
        path = ROOT_PATH + "/steps/{sp}/{chrom_type}/round_{i}".format(sp=sp,chrom_type=chrom_type,i=i+1)
        gwf.target_from_template(
        name = "round_{i}_{sp}_{chrom_type}".format(i=i+1,sp=sp,chrom_type=chrom_type),
        template = run_psmc_bootstrap(path,data_path,sp,chrom_type,i+1,mu,g)
        )
