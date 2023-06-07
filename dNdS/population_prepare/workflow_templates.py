from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/population_dnds/logs"

def bam2vcf(path,bam,ref,name):
    inputs = [bam]
    outputs = [
    LOG_PATH+"/bam2vcf_{name}.DONE".format(name = name)
    ]
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'24:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate platypus
    echo "start vcf 2 bam"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    platypus callVariants --bamFiles={bam} --refFile={ref} --output={name}.vcf
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,name=name,bam=bam,ref=ref,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def vcf_noindel(path,vcfin,vcfout,name):
    inputs = [
    LOG_PATH+"/bam2vcf_{name}.DONE".format(name = name)
    ]
    outputs = [
    LOG_PATH+"/noindel_{name}.DONE".format(name = name)
    ]
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate vcftools
    echo "remove indels in the vcf files"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    vcftools --vcf {vcfin} --remove-indels --recode --out {vcfout}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,name=name,vcfin=vcfin,vcfout=vcfout,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def vcf2ref(path,vcfin,refin,refout,name):
    inputs = [
    LOG_PATH+"/noindel_{name}.DONE".format(name = name)
    ]
    outputs = [
    LOG_PATH+"/altref_{name}.DONE".format(name = name)
    ]
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "make alternative reference fasta from vcfs"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    bgzip -c {vcfin} > {name}.noindel.vcf.gz
    bcftools index {name}.noindel.vcf.gz
    bcftools consensus -f {refin} --sample {name} {name}.noindel.vcf.gz > {refout}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(name=name,path=path,refin=refin,vcfin=vcfin,refout=refout,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def agat_retrieve_coding_fasta(path,refin,gffin,faout,name):
    inputs = [
    LOG_PATH+"/altref_{name}.DONE".format(name = name)
    ]
    outputs = [
    LOG_PATH+"/getcds_{name}.DONE".format(name = name)
    ]
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate agat
    echo "Retrive cds from refernce based on annotations"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    agat_sp_extract_sequences.pl -g {gffin} -f {refin} -t cds -o {faout}.temp 
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/pep_reformat_sp.py {sp} {faout}.temp {faout}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(name=name,path=path,gffin=gffin,refin=refin,faout=faout,sp=name,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def concat_fasta(path,fain,faout,name):
    inputs = [
    fain
    ]
    outputs = [
    LOG_PATH+"/catcds_{name}.DONE".format(name = name)
    ]
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "concat cds from refernce based on annotations"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/population_dnds/concat_fasta.py {fain} {faout}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(name=name,path=path,fain=fain,faout=faout,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
