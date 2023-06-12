from gwf import *

LOG_PATH="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/PSMC/logs"

def logs_sum(logs,sum_log):
    inputs = logs
    outputs = [sum_log]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    spec = """
    echo "finished" > {sum_log}
    echo date
""".format(sum_log=sum_log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def get_consensus_fq_chr(path,bam_in,chrom,fa_ref,name):
    inputs = []
    outputs = [LOG_PATH+"/consensus_fq_{name}_{chrom}.DONE".format(name=name,chrom=chrom)]
    options = {
              'cores': 1,
              'memory': '128g',
              'walltime':'72:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo jobinfo $SLURM_JOBID
    echo "start preparing psmcfa by finding consensus fq file"
    mkdir -p {path}
    cd {path}
    # general 
    /home/jilong/software/samtools-1.12/samtools mpileup -Q 30 -q 30 -u -v -f {ref} -r {chrom} {bam_in} |  bcftools call -c -V indels |  vcfutils.pl vcf2fq -d 5 -D 80 -Q 30 > {name}_{chrom}.fq
    # pacificus
    # /home/jilong/software/samtools-1.12/samtools mpileup -Q 20 -q 30 -u -v -f {ref} -r {chrom} {bam_in} |  bcftools call -c -V indels |  vcfutils.pl vcf2fq -d 5 -Q 20 > {name}_{chrom}.fq
    echo done > {log}
    """.format(name=name,path=path,ref=fa_ref,chrom=chrom,bam_in=bam_in,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cat_files(file_list,logs,combine_name,path):
    file_string = ""
    for filename in file_list:
        file_string = file_string + " {filename}".format(filename=filename)
    inputs = logs
    outputs = [LOG_PATH+"/cat_{combine_name}.DONE".format(combine_name=combine_name)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start concatenating a list of files"
    mkdir -p {path}
    cd {path}
    cat {file_string} > {combine_name}
    echo done > {log}
    """.format(path=path,file_string=file_string,combine_name=combine_name,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def fq2psmcfa(path,fq_file,combine_name):
    inputs = [LOG_PATH+"/cat_{combine_name}.DONE".format(combine_name=fq_file)]
    outputs = [LOG_PATH+"/fq2psmcfa_{combine_name}.DONE".format(combine_name=combine_name)]
    options = {
              'cores': 1,
              'memory': '2g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate psmc
    echo jobinfo $SLURM_JOBID
    echo "start converting fq to psmcfa"
    mkdir -p {path}
    cd {path}
    fq2psmcfa {fq_file} > {combine_name}
    echo done > {log}
    """.format(path=path,fq_file=fq_file,combine_name=combine_name,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def psmc(path,psmcfa,psmc_out):
    inputs = [LOG_PATH+"/fq2psmcfa_{combine_name}.DONE".format(combine_name=psmcfa)]
    outputs = [LOG_PATH+"/psmc_{combine_name}.DONE".format(combine_name=psmc_out)]
    options = {
              'cores': 1,
              'memory': '4g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate psmc
    echo jobinfo $SLURM_JOBID
    echo "start running psmc"
    mkdir -p {path}
    cd {path}
    psmc -p "4+25*2+4+6" -o {psmc_out} {psmcfa}
    echo done > {log}
    """.format(path=path,psmcfa=psmcfa,psmc_out=psmc_out,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

