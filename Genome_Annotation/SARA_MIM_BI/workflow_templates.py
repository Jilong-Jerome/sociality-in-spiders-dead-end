from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/logs" 
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
def STAR_index(path,ref,log):
    """STAR genome index generate"""
    inputs = [ref]
    outputs = [LOG_PATH + "/"+log]
    options = {
               'cores': 18,
               'memory': '64g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate STAR
    echo jobinfo $SLURM_JOBID
    echo "start STAR indexing"
    date
    mkdir -p {path}
    cd {path}
    STAR --runThreadN 18 --runMode genomeGenerate --genomeDir {path} --genomeFastaFiles {ref}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,ref=ref,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def STAR_mapping(align_path,index_path,sp,fq1,fq2,RG):
    """STAR RNA alignment"""
    if sp == "MIM16":
        new_RG = "MIM16_{RG}".format(RG=RG)
    else:
        new_RG = RG
    inputs = [LOG_PATH+"/star_index_{sp}.DONE".format(sp=sp)]
    outputs = [LOG_PATH+"/STAR_align_"+new_RG+".DONE"]
    options = {
               'cores': 16,
               'memory': '48g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate STAR
    echo jobinfo $SLURM_JOBID
    echo "start star align {RG}"
    date
    mkdir -p {align_path}/{RG}
    STAR --genomeDir {index_path} --runThreadN 16 --readFilesCommand zcat --readFilesIn {fq1} {fq2} --outSAMattrRGline ID:{RG} --outFileNamePrefix {align_path}/{RG}/{RG} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(align_path=align_path,index_path=index_path,fq1=fq1,fq2=fq2,RG=RG,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def merge_bams(path,outname,files,log_list):
    """samtools merge multiple bam files into one"""
    inputs = log_list
    outputs = [LOG_PATH+"/STAR_merge_"+outname+".DONE"]
    options = {
               'cores': 16,
               'memory': '48g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start merging bam files with samtools"
    date
    mkdir -p {path}
    /home/jilong/software/samtools-1.12/samtools merge -@ 16 -O BAM {path}/{outname}.bam {files}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,outname=outname,files=files,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def sort_bam(path,bam_in,bam_out):
    """samtools merge multiple bam files into one"""
    inputs = [LOG_PATH+"/STAR_merge_"+bam_in+".DONE"]
    outputs = [LOG_PATH+"/STAR_merge_"+bam_out+"_sort.DONE"]
    options = {
               'cores': 16,
               'memory': '16g',
               'walltime':"12:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start merging bam files with samtools"
    date
    mkdir -p {path}
    cd {path}
    /home/jilong/software/samtools-1.12/samtools sort -@ 16 -O BAM -o {bam_out}.bam {bam_in}.bam
    /home/jilong/software/samtools-1.12/samtools index -@ 16 -b {bam_out}.bam {bam_out}.bai
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,bam_in=bam_in,bam_out=bam_out,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def repeat_model(genome,out,path):
    inputs = [genome]
    outputs = [LOG_PATH+"/repeat_model_{out}.DONE".format(out=out)]
    options = {
              'cores':24,
              'memory':'64g',
              'walltime':'72:00:00',
              'account':'spider2'
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start RepeatModeler"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    BuildDatabase -name {out} -engine ncbi {genome}
    RepeatModeler -pa 24 -engine ncbi -database {out}
    echo finished > {log}
    echo "finished RepeatModeler"
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,out=out,genome = genome,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def repeat_masker_lib(path,fasta,lib1,lib2,lib,sp):
    inputs = [LOG_PATH+"/repeat_model_{sp}.DONE".format(sp=sp)]
    outputs = [LOG_PATH+"/repeat_mask_{sp}.DONE".format(sp=sp)]
    options = {
              'cores':24,
              'memory':'64g',
              'walltime':'72:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start RepeatMasker"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    cat {lib1} {lib2} > {lib}
    RepeatMasker -e ncbi -pa 24 -xsmall -dir {path} -lib {lib} {fasta}
    echo finished > {log}
    date
    echo jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,lib1=lib1,lib2=lib2,lib=lib,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def seqkit_grep(sp,path,genome,chrom_id,out):
    inputs = [LOG_PATH+"/repeat_mask_{sp}.DONE".format(sp=sp)]
    outputs = [
              LOG_PATH+"/seqkit_grep_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
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
    seqkit grep -p HiC_scaffold_{chrom_id} {genome} -o {out}
    echo finished > {log}
    echo "finished seqkit"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,chrom_id=chrom_id,out=out,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def braker_combine_mode(path,sp,aug_sp,chrom_id,genome,protein,bam):
    inputs = [
        LOG_PATH+"/seqkit_grep_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id) 
    ]
    outputs = [
              LOG_PATH+"/BRAKER_ETP_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
    ]
    options = {
              'cores':12,
              'memory':'72g',
              'walltime':'48:00:00',
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
    braker.pl --species={aug_sp} --genome={genome} --prot_seq={protein} --bam {bam} --etpmode --softmasking --cores=12 --gff3 > {path}/braker_{aug_sp}.log.tmp
    mv {path}/braker_{aug_sp}.log.tmp {path}/braker_{aug_sp}.log
    echo finished > {log}
    echo "finished braker"
    date
    jobinfo $SLURM_JOBID
    """.format(aug_sp=aug_sp,log=outputs[0],path=path,genome=genome,protein=protein,bam=bam,sp=sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def braker_protein_mode(path,sp,aug_sp,chrom_id,genome,protein,bam):
    inputs = [
        LOG_PATH+"/seqkit_grep_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id) 
    ]
    outputs = [
              LOG_PATH+"/BRAKER_protein_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
    ]
    options = {
              'cores':12,
              'memory':'72g',
              'walltime':'48:00:00',
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
    braker.pl --species={aug_sp} --genome={genome} --prot_seq={protein} --epmode --cores 12 --gff3 --softmasking > {path}/braker_{aug_sp}.log.tmp
    mv {path}/braker_{aug_sp}.log.tmp {path}/braker_{aug_sp}.log
    echo finished > {log}
    echo "finished braker"
    date
    jobinfo $SLURM_JOBID
    """.format(aug_sp=aug_sp,log=outputs[0],path=path,genome=genome,protein=protein,sp=sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
    #braker.pl --species={aug_sp} --genome={genome} --prot_seq={protein} --bam {bam} --softmasking --cores=12 --trainFromGth --prg=gth --ALIGNMENT_TOOL_PATH=/home/jilong/software/gth-1.7.1-Linux_x86_64-64bit/bin/gth --gff3 > {path}/braker_{aug_sp}.log.tmp

def agat_gff2fasta(path,gff,genome,fasta,sp,chrom_id):
    inputs = [LOG_PATH+"/BRAKER_ETP_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)]
    outputs = [
              LOG_PATH+"/gff2fasta_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
              ]
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate agat
    echo "start finding best match protein"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    agat_sp_extract_sequences.pl -g {gff} -f {genome} -o {fasta}.tmp -t gene
    mv {fasta}.tmp {fasta}
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,gff=gff,fasta=fasta,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def agat_gff2fasta_prot(path,gff,genome,fasta,sp,chrom_id):
    inputs = [LOG_PATH+"/BRAKER_ETP_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)]
    outputs = [
              LOG_PATH+"/gff2prot_fasta_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
              ]
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate agat
    echo "start finding best match protein"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    agat_sp_extract_sequences.pl -g {gff} -f {genome} -o {fasta}.tmp -p
    mv {fasta}.tmp {fasta}
    date
    echo finished > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,genome=genome,gff=gff,fasta=fasta,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cat_braker_sp(sp,path,string_gff,string_n,string_p):
    inputs = [LOG_PATH+"/gff2fa_{sp}.DONE".format(sp=sp)]
    outputs = [LOG_PATH+"/braker_result_combine_{sp}_.DONE".format(sp=sp)]
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo "start finding best match protein"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    cat {string_gff} > {sp}_braker.gff3
    cat {string_n} > {sp}_braker.fasta
    cat {string_p} > {sp}_braker.protein.fasta
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/full_annotation/rename_naive.py {sp}_braker.fasta {sp}_braker.rename.fasta 1 0 
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/full_annotation/rename_naive.py {sp}_braker.protein.fasta {sp}_braker.protein.rename.fasta 2 0 
    date
    echo finished > {log}
    jobinfo $SLURM_JOBID
    """.format(path=path,string_gff=string_gff,string_n=string_n,string_p=string_p,sp=sp,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def run_busco_protein(sp,fasta,path,outname):
    inputs=[LOG_PATH+"/braker_result_combine_{sp}_.DONE".format(sp=sp)]
    outputs=[LOG_PATH+"/braker_result_BUSCO_{sp}_{outname}.DONE".format(sp=sp,outname=outname)]
    options={
        'cores':16,
        'memory':'32g',
        'walltime':'12:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate busco
    echo "start busco for checking completeness"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    busco -i {fasta} -o {outname} -l arthropoda_odb10 -m proteins -c 16 -f
    echo finished > {log}
    """.format(fasta=fasta,path=path,outname=outname,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def braker_rna_mode(path,sp,aug_sp,chrom_id,genome,protein,bam):
    inputs = [
        LOG_PATH+"/seqkit_grep_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id) 
    ]
    outputs = [
              LOG_PATH+"/BRAKER_rna_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id)
    ]
    options = {
              'cores':12,
              'memory':'72g',
              'walltime':'48:00:00',
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
    braker.pl --species={aug_sp} --genome={genome} --bam {bam} --softmasking --cores=12 --gff3 > {path}/braker_{aug_sp}.log.tmp
    mv {path}/braker_{aug_sp}.log.tmp {path}/braker_{aug_sp}.log
    echo finished > {log}
    echo "finished braker"
    date
    jobinfo $SLURM_JOBID
    """.format(aug_sp=aug_sp,log=outputs[0],path=path,genome=genome,bam=bam,sp=sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
