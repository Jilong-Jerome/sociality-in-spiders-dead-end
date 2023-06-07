from gwf import *
def diamond_build_db(path,fasta,out):
    inputs = [fasta]
    outputs = [
              path+"/{out}.log".format(out=out)]
    options = {
              'cores':8,
              'memory':'16g',
              'walltime':'12:00:00',
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
    diamond makedb --in {fasta} -d {out} > {out}.log.tmp
    mv {out}.log.tmp {out}.log
    echo "finished diamond makedb"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,out=out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def diamond_blastp(path,db,query,out):
    inputs = [fasta]
    outputs = [
              path+"/{out}.log".format(out=out)]
    options = {
              'cores':24,
              'memory':'24g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate braker
    echo "start diamond blastp"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    diamond blastp -p 24 -k 1000 -q {query} -d {db} -o {out}.tsv --ultra-sensitive > {out}.log.tmp
    mv {out}.log.tmp {out}.log
    echo "finished diamond blastp"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,query=query,db=db,out=out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def diamond_blastp_1(path,db,query,out):
    inputs = [fasta]
    outputs = [
              path+"/{out}.log".format(out=out)]
    options = {
              'cores':24,
              'memory':'24g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate braker
    echo "start diamond blastp"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    diamond blastp -p 24 -k 1 -q {query} -d {db} -o {out}.tsv --ultra-sensitive > {out}.log.tmp
    mv {out}.log.tmp {out}.log
    echo "finished diamond blastp"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,query=query,db=db,out=out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def diamond_blastx_n(path,n,db,query,out):
    inputs = [fasta]
    outputs = [
              path+"/{out}.log".format(out=out)]
    options = {
              'cores':24,
              'memory':'24g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate braker
    echo "start diamond blastp"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    diamond blastx -p 24 -k {n} -q {query} -d {db} -o {out}.tsv --ultra-sensitive > {out}.log.tmp
    mv {out}.log.tmp {out}.log
    echo "finished diamond blastp"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,query=query,db=db,out=out,n=n)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def count_by_number(in_aln,outname):
    inputs = [in_aln]
    outputs = [outname]
    options = {
              'cores':1,
              'memory':'100g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate gross
    echo "start counting n"
    echo jobinfo $SLURM_JOBID
    date
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/diamond/count_by_repeat_times.R {in_aln} {outname}.tmp
    mv {outname}.tmp {outname}
    date
    jobinfo $SLURM_JOBID
    """.format(in_aln=in_aln,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


gwf = Workflow()
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/single_ortho_sequences"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/single_ortho_sequences/single_copy_pep.fa"
out = "stegodyphus_single_copy"
gwf.target_from_template(
    name="diamond_makedb_{out}".format(out=out),
    template = diamond_build_db(
        path=path,
        fasta=fasta,
        out = out
        )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/single_ortho_sequences"
query = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/trinity/africanus/africanums_trinity.Trinity.fasta"
db = "stegodyphus_single_copy"
out = "africanus_hit"
n = 1
gwf.target_from_template(
    name="diamond_blasx_{out}".format(out=out),
    template = diamond_blastx_n(
        path=path,
        query = query,
        n = n,
        db = db,
        out = out
        )
    )
