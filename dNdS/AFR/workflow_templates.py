from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/trinity/logs" 

def trinity_assembly(path,read1,read2,outname):
    inputs = [read1,read2]
    outputs = [LOG_PATH+"{outname}.DONE".format(outname=outname)]
    options = {
              'cores':32,
              'memory':'240g',
              'walltime':'96:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate trinity
    echo "start trinity transcriptome assembly"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    Trinity --seqType fq --trimmomatic --max_memory 250G --left {read1} --right {read2} --CPU 32 --output "{outname}"
    echo "finished trinity assembling" > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(read1=read1,read2=read2,outname=outname,path=path,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
