from gwf import *

LOG_PATH = "/home/jilong/EcoGenetics/people/jilong/psmc/logs"

def split_psmcfa(data_path,sp,chrom_type):
    inputs = ["{data_path}/{sp}_{chrom_type}.psmcfa".format(data_path=data_path,sp=sp,chrom_type=chrom_type)]
    outputs = [LOG_PATH + "/{sp}_{chrom_type}_split_psmcfa.DONE".format(sp=sp,chrom_type=chrom_type)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate psmc
    echo jobinfo $SLURM_JOBID
    echo "spliting psmcfa from chromsome level to smaller contigs"
    splitfa {data_path}/{sp}_{chrom_type}.psmcfa 100000 > {data_path}/{sp}_{chrom_type}_split.psmcfa
    echo done > {log}
    """.format(data_path=data_path,sp=sp,chrom_type=chrom_type,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def run_psmc_bootstrap(path,data_path,sp,chrom_type,i,mu,g):
    inputs = [LOG_PATH + "/{sp}_{chrom_type}_split_psmcfa.DONE".format(sp=sp,chrom_type=chrom_type)]
    outputs = [LOG_PATH+"/{sp}_iterations/{sp}_{chrom_type}_{i}.DONE".format(sp=sp,chrom_type=chrom_type,i=i)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate psmc
    echo jobinfo $SLURM_JOBID
    echo "run bootstrapping psmc"
    mkdir -p {path}
    cd {path}
    cp {data_path}/{sp}_{chrom_type}_split.psmcfa .
    psmc -b -p "4+25*2+4+6" -o {sp}_{chrom_type}_round_{i}.psmc {sp}_{chrom_type}_split.psmcfa
    psmc_plot.pl -R -u {mu} -g {g} {sp}_{chrom_type}_round_{i} {sp}_{chrom_type}_round_{i}.psmc
    sed -i 's/^/{sp}\t{chrom_type}\t{i}\t/' {sp}_{chrom_type}_round_{i}.0.txt
    mkdir -p {log}
    rmdir {log}
    rm {sp}_{chrom_type}_split.psmcfa
    echo done > {log}
    """.format(log = outputs[0],path=path,data_path=data_path,sp=sp,chrom_type=chrom_type,i=i,mu=mu,g=g)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
