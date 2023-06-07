from gwf import *
from workflow_dicts import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/boostrap/logs"
def random_sample(path, target, n, out):
    inputs = [target]
    outputs = [LOG_PATH+"/sample_{out}.DONE".format(out=out),
              path+"/{out}_id.txt".format(out=out)]
    options = {
              'cores':1,
              'memory':'100m',
              'walltime':'00:05:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate paml
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    shuf -n {n} {target} |cut -f1 >  {out}_id.txt 
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(out=out,path=path,target=target,log=outputs[0],n=n)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def retrive_pop_fa(pop_cds,idname,sample_id,outname,path):
    inputs = [LOG_PATH+"/sample_{out}.DONE".format(out=idname)]
    outputs = [LOG_PATH+"/grep_fa_{out}.DONE".format(out=outname)]
    options = {
              'cores':1,
              'memory':'1g',
              'walltime':'01:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    seqkit grep -f {sample_id} {pop_cds} -o {outname}_split.fa
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/population_dnds/concat_fasta.py {outname}_split.fa {outname}.fa
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(log=outputs[0],path=path,sample_id=sample_id,pop_cds=pop_cds,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def fa2phy(path,prefix):
    inputs = [LOG_PATH+"/grep_fa_{out}_pop1.DONE".format(out=prefix),
              LOG_PATH+"/grep_fa_{out}_pop2.DONE".format(out=prefix)] 
    outputs = [LOG_PATH+"/fa2phy_{out}.DONE".format(out=prefix)]
    options = {
              'cores':1,
              'memory':'1g',
              'walltime':'01:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    cat {prefix}_pop1.fa {prefix}_pop2.fa > {prefix}_NT.fa
    conda activate clustalo
    trimal -in {prefix}_NT.fa -phylip3.2 -out {prefix}_NT.phy
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(log=outputs[0],path=path,prefix=prefix)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def run_paml(path,prefix):
    inputs = [LOG_PATH+"/fa2phy_{out}.DONE".format(out=prefix)]
    outputs = [LOG_PATH+"/paml_{out}.DONE".format(out=prefix)]
    options = {
              'cores':1,
              'memory':'5g',
              'walltime':'04:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    echo "({prefix}_pop1,{prefix}_pop2)#" > tree.txt
    ln -sf ../{prefix}_NT.phy seq.txt
    cp /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/codeml/DUM_auto/codeml.ctl .
    conda activate paml
    codeml
    tail -n 1 results.txt > {prefix}_raw.tab
    date
    jobinfo $SLURM_JOBID
    echo done > {log}
    """.format(log=outputs[0],path=path,prefix=prefix)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


gwf = Workflow()
for i in range(500):
    n = 500
    for sp in ["SARA","DUM","MIM"]:
        out = "{sp}_{n}_{j}".format(sp=sp,n=n,j=i+1)
        target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/population_dnds/{sp}_auto_single_id.tsv".format(sp=sp)
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/boostrap/{sp}/{out}".format(sp=sp,out=out)
        gwf.target_from_template(
            name = "sample_{out}".format(out=out),
            template = random_sample(path, target, n, out)
        )
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/boostrap/{sp}/{out}".format(sp=sp,out=out)
        pop_cds = cds_dicts[sp][0]
        sample_id = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/boostrap/{sp}/{out}/{out}_id.txt".format(sp=sp,out=out)
        outname = "{out}_pop1".format(out=out)
        gwf.target_from_template(
            name = "grep_fa_{out}".format(out=outname),
            template = retrive_pop_fa(pop_cds,out,sample_id,outname,path)
        )
        pop_cds = cds_dicts[sp][1]
        outname = "{out}_pop2".format(out=out)
        gwf.target_from_template(
            name = "grep_fa_{out}".format(out=outname),
            template = retrive_pop_fa(pop_cds,out,sample_id,outname,path)
        )
        gwf.target_from_template(
            name = "fa2phy_{out}".format(out=out),
            template = fa2phy(path,out)
        )
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/boostrap/{sp}/{out}/codeml".format(sp=sp,out=out)
        gwf.target_from_template(
            name = "paml_{out}".format(out=out),
            template = run_paml(path,out)
        )
