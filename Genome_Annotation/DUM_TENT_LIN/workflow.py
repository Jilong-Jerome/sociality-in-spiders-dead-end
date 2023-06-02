from gwf import *
def braker_combine_mode(path,genome,protein,bam,species):
    inputs = [protein,bam]
    outputs = [
              path+"/braker_{species}.log".format(species=species)]
    options = {
              'cores':12,
              'memory':'100g',
              'walltime':'15:00:00',
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
    braker.pl --species={species} --genome={genome} --prot_seq={protein} --bam {bam} --etpmode --softmasking --cores=12 --gff3 > {path}/braker_{species}.log.tmp
    mv {path}/braker_{species}.log.tmp {path}/braker_{species}.log
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,protein=protein,bam=bam,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def braker_combine_mode_full(path,genome,protein,bam,species):
    inputs = [protein,bam]
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
    braker.pl --species={species} --genome={genome} --prot_seq={protein} --bam {bam} --etpmode --softmasking --cores=24 --gff3 > {path}/braker_{species}.log.tmp
    mv {path}/braker_{species}.log.tmp {path}/braker_{species}.log
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,protein=protein,bam=bam,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def seqkit_grep(path,genome,chrom_id,out):
    inputs = [path+"/"+genome]
    outputs = [
              path+"/"+out
              ]
    options = {
              'cores':1,
              'memory':'4g',
              'walltime':'1:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    echo "start get soft masked chromosome"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    seqkit grep -p HiC_scaffold_{chrom_id} {genome} -o {out}.tmp
    mv {out}.tmp {out}
    echo "finished seqkit"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,chrom_id=chrom_id,out=out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



gwf = Workflow()

for chrom_id in range(16):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine"
    genome = "DUM_hifi_hic_scaffolded_trim.fa.masked"
    out = "DUM_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    task_name = "DUM_get_fa_{chrom_id}".format(chrom_id=chrom_id)
    gwf.target_from_template(
        name=task_name,
        template = seqkit_grep(path,genome,chrom_id,out)
        )


#DUM
for chrom_id in range(16):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/etp_mode/HiC_scaffold_{chrom_id}".format(chrom_id=chrom_id)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine/DUM_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
    bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam"
    species = "dumicola_braker_etp_{chrom_id}".format(chrom_id=chrom_id)
    if chrom_id != 11:
        gwf.target_from_template(
            name="DUM_braker_etp_mode_{chrom_id}".format(chrom_id=chrom_id),
            template = braker_combine_mode(
                path=path,
                genome=genome,
                protein=protein,
                bam = bam,
                species=species
                )
            )
#TENT
for chrom_id in range(16):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase/combine"
    genome = "TENT_hifi_hic_scaffolded_trim.fa.masked"
    out = "TENT_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    task_name = "TENT_get_fa_{chrom_id}".format(chrom_id=chrom_id)
    gwf.target_from_template(
        name=task_name,
        template = seqkit_grep(path,genome,chrom_id,out)
        )
for chrom_id in range(16):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/tentoriicola/etp_mode/HiC_scaffold_{chrom_id}".format(chrom_id=chrom_id)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase/combine/TENT_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
    bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/TENT_RNA_STAR_sort.bam"
    species = "TENT_braker_etp_{chrom_id}".format(chrom_id=chrom_id)
    gwf.target_from_template(
        name="TENT_braker_etp_mode_{chrom_id}".format(chrom_id=chrom_id),
        template = braker_combine_mode(
            path=path,
            genome=genome,
            protein=protein,
            bam = bam,
            species=species
            )
        )
# LIN
for chrom_id in range(25):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase/combine"
    genome = "LIN_hifi_hic_scaffolded_trim.fa.masked"
    out = "LIN_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    task_name = "LIN_get_fa_{chrom_id}".format(chrom_id=chrom_id)
    gwf.target_from_template(
        name=task_name,
        template = seqkit_grep(path,genome,chrom_id,out)
        )
for chrom_id in range(25):
    chrom_id = chrom_id + 1
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/lineatus/etp_mode/HiC_scaffold_{chrom_id}".format(chrom_id=chrom_id)
    genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase/combine/LIN_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(chrom_id=chrom_id)
    protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
    bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/LIN_RNA_STAR_sort.bam"
    species = "LIN_braker_etp_{chrom_id}".format(chrom_id=chrom_id)
    gwf.target_from_template(
        name="LIN_braker_etp_mode_{chrom_id}".format(chrom_id=chrom_id),
        template = braker_combine_mode(
            path=path,
            genome=genome,
            protein=protein,
            bam = bam,
            species=species
            )
        )

