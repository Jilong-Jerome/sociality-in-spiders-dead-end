from gwf import *
from workflow_dicts import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/logs"
OUT_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace"
def logs_sum(logs,sum_log):
    inputs = logs
    outputs = [sum_log]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    echo "finished" > {sum_log}
    echo date
""".format(sum_log=sum_log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def format_gff(ste,sp,chrom,gff_in):
    inputs = []
    outputs = [LOG_PATH+"/{sp}_{ste}_{chrom}_gff_reform.DONE".format(sp=sp,chrom=chrom,ste=ste)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    path = OUT_PATH+"/{ste}/rawGenomes/{sp}/{sp}/annotation".format(ste=ste,sp=sp)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ortho
    echo jobinfo $SLURM_JOBID
    echo "start genespace reformating gff file"
    mkdir -p {path}
    cd {path}
    cat {gff_in} | grep "gene" > {sp}_{chrom}_gene.gff3
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/genespace/gff_reformat_sp.py {sp} {sp}_{chrom}_gene.gff3 {sp}_{chrom}_gene.gff
    rm {sp}_{chrom}_gene.gff3
    echo done > {log}
    """.format(path=path,gff_in=gff_in,sp=sp,chrom=chrom,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def format_pep_fa(ste,sp,chrom,pep_in):
    inputs = []
    outputs = [LOG_PATH+"/{sp}_{ste}_{chrom}_pep_fa_reform.DONE".format(sp=sp,ste=ste,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    path = OUT_PATH+"/{ste}/rawGenomes/{sp}/{sp}/annotation".format(ste=ste,sp=sp)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    echo "start genespace reformating pep fasta file"
    mkdir -p {path}
    cd {path}
    cat {pep_in} | seqkit grep -r -p {chrom} -p ".t1" > {sp}_{chrom}_pep.preparing
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/genespace/pep_reformat_sp.py {sp} {sp}_{chrom}_pep.preparing {sp}_{chrom}_pep.fa
    rm {sp}_{chrom}_pep.preparing
    echo done > {log}
    """.format(path=path,pep_in=pep_in,sp=sp,chrom=chrom,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cat_group_gff(sp,ste,gff_list):
    inputs = [LOG_PATH+"/gff_{sp}_{ste}.DONE".format(sp=sp,ste=ste)]
    outputs = [LOG_PATH+"/cat_gff_{sp}_{ste}.DONE".format(sp=sp,ste=ste)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    cat_string = ""
    rm_string = ""
    for gff in gff_list:
        cat_string = cat_string + gff + " "
        rm_string = rm_string + "rm " + gff + " \n"
    path = OUT_PATH+"/{ste}/rawGenomes/{sp}/{sp}/annotation".format(ste=ste,sp=sp)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    echo "start concatnating  reformated gff  file"
    mkdir -p {path}
    cd {path}
    cat {cat_string} > {sp}_{ste}_gene.gff
    {rm_string}
    echo done > {log}
    """.format(path=path,cat_string=cat_string,rm_string=rm_string,sp=sp,ste=ste,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def cat_group_pep(sp,ste,pep_list):
    inputs = [LOG_PATH+"/pep_{sp}_{ste}.DONE".format(sp=sp,ste=ste)]
    outputs = [LOG_PATH+"/cat_pep_{sp}_{ste}.DONE".format(sp=sp,ste=ste)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    cat_string = ""
    rm_string = ""
    for pep in pep_list:
        cat_string = cat_string + pep + " "
        rm_string = rm_string + "rm " + pep + " \n"
    path = OUT_PATH+"/{ste}/rawGenomes/{sp}/{sp}/annotation".format(ste=ste,sp=sp)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    echo "start concatnating  reformated pep  file"
    mkdir -p {path}
    cd {path}
    cat {cat_string} > {sp}_{ste}_pep.fa
    {rm_string}
    echo done > {log}
    """.format(path=path,cat_string=cat_string,rm_string=rm_string,sp=sp,ste=ste,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
    
def genespace_ste(ste):
#srun --account spider2 --mem 16g -t 12:00:00 Rscript run_genespace.R
    inputs = [LOG_PATH+"/gff_{ste}.DONE".format(ste=ste),
             LOG_PATH+"/pep_{ste}.DONE".format(ste=ste)]
    outputs = [LOG_PATH+"/genespace_{ste}.DONE".format(ste=ste)]
    options = {
              'cores': 8,
              'memory': '24g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    path = OUT_PATH+"/{ste}".format(ste=ste)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ortho
    echo jobinfo $SLURM_JOBID
    echo "start concatnating  reformated pep  file"
    mkdir -p {path}
    cd {path}
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/genespace/run_genespace.R {ste}
    echo done > {log}
    """.format(path=path,ste=ste,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def draw_genespace_ste(ste):
#srun --account spider2 --mem 16g -t 12:00:00 Rscript run_genespace.R
    inputs = [LOG_PATH+"/genespace_{ste}.DONE".format(ste=ste)]
    outputs = [LOG_PATH+"/draw_genespace_{ste}.DONE".format(ste=ste)]
    options = {
              'cores': 8,
              'memory': '24g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    path = OUT_PATH+"/{ste}".format(ste=ste)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ortho
    echo jobinfo $SLURM_JOBID
    echo "start concatnating  reformated pep  file"
    mkdir -p {path}
    cd {path}
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/genespace/draw_genespace.R {ste}
    echo done > {log}
    """.format(path=path,ste=ste,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
    
def find_OF_fasta_sp(ste,OF,sp,cluster):
    sp_id = "Species{num}".format(num=of_fasta_dict[sp])
    inputs = [LOG_PATH+"/genespace_{ste}.DONE".format(ste=ste)]
    outputs = [LOG_PATH+"/getfa_{ste}_{sp}_{OF}.DONE".format(ste=ste,sp=sp,OF=OF)]
    options = {
              'cores': 1,
              'memory': '8g',
              'walltime':'4:00:00',
              'account':"spider2"
              }
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/region_seq/{sp}/{cluster}/genes".format(sp=sp,cluster=cluster)
    fasta_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/{ste}/orthofinder/Results_*/WorkingDirectory".format(ste=ste)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    echo "start retreiveing fasta sequence from genespace orthofinder results"
    mkdir -p {path}
    cd {path}
    cat {fasta_path}/{sp_id}.fa | seqkit grep -n -r -p {OF} > {OF}.fa 
    echo done > {log}
    """.format(fasta_path=fasta_path,path=path,sp_id=sp_id,OF=OF,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def combine_cluster_fasta(pre_logs,sp,cluster):
    inputs = pre_logs
    outputs = [LOG_PATH+"/combine_fa_{cluster}.DONE".format(cluster=cluster)]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/region_seq/{sp}/{cluster}".format(sp=sp,cluster=cluster)
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    echo "start retreiveing fasta sequence from genespace orthofinder results"
    mkdir -p {path}
    cd {path}
    cat ./genes/*.fa > {cluster}.fa 
    echo done > {log}
    """.format(path=path,sp=sp,cluster=cluster,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def diamond_build_db(path,fasta,cluster):
    inputs = [LOG_PATH+"/combine_fa_{cluster}.DONE".format(cluster=cluster)]
    outputs = [
              LOG_PATH+"/diamond_DB_{cluster}.DONE".format(cluster=cluster)]
    options = {
              'cores':4,
              'memory':'8g',
              'walltime':'4:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate braker
    echo "start diamond makedb"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    diamond makedb --in {fasta} -d {out} 
    echo "finished diamond makedb" > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,out=cluster,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

 
