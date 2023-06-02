from gwf import *
def braker_rna_mode(path,genome,bam,species):
    inputs = [bam]
    outputs = [
              path+"/braker_{species}.log".format(species=species)]
    options = {
              'cores':24,
              'memory':'200g',
              'walltime':'96:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate BRAKER
    echo "start braker pipe"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    braker.pl --species={species} --genome={genome} --bam={bam} --softmasking on --cores=24 --gff3 > {path}/braker_{species}.log.tmp
    mv {path}/braker_{species}.log.tmp {path}/braker_{species}.log
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,bam=bam,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def braker_rna_mode_hard(path,genome,bam,species):
    inputs = [bam]
    outputs = [
              path+"/braker_{species}.log".format(species=species)]
    options = {
              'cores':24,
              'memory':'200g',
              'walltime':'96:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate BRAKER
    echo "start braker pipe"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    braker.pl --species={species} --genome={genome} --bam={bam} --cores=24 --gff3 > {path}/braker_{species}.log.tmp
    mv {path}/braker_{species}.log.tmp {path}/braker_{species}.log
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,bam=bam,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def braker_check_evidence(gtf,hint,path,log):
    inputs = [gtf,hint]
    outputs = [
              path+"/"+log]
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate BRAKER
    echo "start braker pipe"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    /home/jilong/software/BRAKER/scripts/predictionAnalysis/selectSupportedSubsets.py {gtf} {hint} --fullSupport full --anySupport any --noSupport no > {path}/{log}.tmp
    mv {path}/{log}.tmp {path}/{log}
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,gtf=gtf,log=log,hint=hint)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



gwf = Workflow()
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine/DUM_hifi_hic_scaffolded_trim.fa.masked"
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam"
species = "dumicola_braker_rna"
gwf.target_from_template(
    name="DUM_braker_rna_mode",
    template = braker_rna_mode(
        path=path,
        genome=genome,
        bam=bam,
        species=species
        )
    )

gtf = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/augustus.gtf"
hint = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/hintsfile.gff"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/rna_mode/braker/dumicola_braker_rna/support"
log = "prediction_analysis.log"
gwf.target_from_template(
    name="DUM_braker_rna_mode_check_evidence",
    template = braker_check_evidence(
        gtf=gtf,
        hint=hint,
        path=path,
        log=log
        )
    )
