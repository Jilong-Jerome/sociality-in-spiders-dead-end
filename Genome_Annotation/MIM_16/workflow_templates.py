from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/augustus_manual/logs" 
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
def blat(path,db,query,outname,log):
    """STAR genome index generate"""
    inputs = []
    outputs = [LOG_PATH + "/"+log]
    options = {
               'cores': 4,
               'memory': '16g',
               'walltime':"48:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate blat
    echo jobinfo $SLURM_JOBID
    echo "start blat"
    date
    mkdir -p {path}
    cd {path}
    blat -noHead {db} {query} {outname}_est.psl    
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,db=db,query=query,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def filter_blat(path,outname,log):
    """filter blat alignment"""
    inputs = [LOG_PATH+"/{outname}_blat.DONE".format(outname=outname)]
    outputs = [LOG_PATH + "/"+log]
    options = {
               'cores': 1,
               'memory': '4g',
               'walltime':"4:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate BRAKER
    echo jobinfo $SLURM_JOBID
    echo "start fileter "
    date
    mkdir -p {path}
    cd {path}
    cat {outname}_est.psl | filterPSL.pl --best > {outname}_est.f.psl
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def augustus(path,fasta,cfg,hints,outname,log):
    """run augustus"""
    inputs = []
    outputs = [LOG_PATH + "/"+log]
    options = {
               'cores': 4,
               'memory': '16g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate BRAKER
    echo jobinfo $SLURM_JOBID
    echo "start fileter "
    date
    mkdir -p {path}
    cd {path}
    augustus --species=MIM_braker_MIM16_part_1 {fasta} --extrinsicCfgFile={cfg} --hintsfile={hints} > {outname}.gff
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,fasta=fasta,cfg=cfg,hints=hints,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
