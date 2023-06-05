from gwf import *
def STAR_index(path,ref,log):
    """STAR genome index generate"""
    inputs = [ref]
    outputs = [path+"/"+log]
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
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir {path} --genomeFastaFiles {ref} > {path}/{log}.tmp 
    mv {path}/{log}.tmp {path}/{log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,ref=ref,log=log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def STAR_index_fast(path,ref,log):
    """STAR genome index generate"""
    inputs = [ref]
    outputs = [path+"/"+log]
    options = {
               'cores': 16,
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
    STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {path} --genomeFastaFiles {ref} > {path}/{log}.tmp 
    mv {path}/{log}.tmp {path}/{log}
    date
    jobinfo $SLURM_JOBID 
    """.format(path=path,ref=ref,log=log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def STAR_mapping(align_path,index_path,index_log,fq1,fq2,RG):
    """STAR RNA alignment"""
    inputs = [fq1,fq2,index_log]
    outputs = ["{align_path}/{RG}/{RG}.log".format(align_path=align_path,RG=RG)]
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
    STAR --genomeDir {index_path} --runThreadN 16 --readFilesCommand zcat --readFilesIn {fq1} {fq2} --outSAMattrRGline ID:{RG} --outFileNamePrefix {align_path}/{RG}/{RG} --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx > {align_path}/{RG}/{RG}.log.tmp
    mv {align_path}/{RG}/{RG}.log.tmp {align_path}/{RG}/{RG}.log
    date
    jobinfo $SLURM_JOBID
    """.format(align_path=align_path,index_path=index_path,fq1=fq1,fq2=fq2,RG=RG)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def merge_bams(path,outname,logs,files):
    """samtools merge multiple bam files into one"""
    files_log = open(logs)
    log_list = []
    for file_id in files_log:
        log_list.append(file_id.strip("\n"))
    inputs = log_list
    outputs = [path+"/"+outname+".bam"]
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
    /home/jilong/software/samtools-1.12/samtools merge -@ 16 -b {files} -O BAM {path}/{outname}.bam.tmp
    mv {path}/{outname}.bam.tmp {path}/{outname}.bam
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,outname=outname,files=files)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def sort_bam(path,bam_in,bam_out,log):
    """samtools merge multiple bam files into one"""
    inputs = [bam_in]
    outputs = [path+"/"+log]
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
    /home/jilong/software/samtools-1.12/samtools sort -@ 16 -O BAM -o {bam_out} {bam_in} > {path}/{log}.tmp
    mv {path}/{log}.tmp {path}/{log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,bam_in=bam_in,bam_out=bam_out,log=log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def index_bam(bam,bai):
    """samtools merge multiple bam files into one"""
    inputs = [bam]
    outputs = [bai]
    options = {
               'cores': 8,
               'memory': '8g',
               'walltime':"12:00:00",
               'account':"spider2"
    } 
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo jobinfo $SLURM_JOBID
    echo "start indexing bam files with samtools"
    date
    /home/jilong/software/samtools-1.12/samtools index -@ 16 -b {bam} {bai}.tmp
    mv {bai}.tmp {bai}
    date
    jobinfo $SLURM_JOBID
    """.format(bam=bam,bai=bai)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

gwf = Workflow()

##build genome index for DUM
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/STAR_index"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/DUM_hifi_hic_scaffolded_trim.fa"
log = "DUM_star_index.log"
species = "DUM"
gwf.target_from_template(
        name="STAR_index_{sp}".format(sp=species),
        template = STAR_index(path,ref,log)
    )
reads_list = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/90-524847844/dumicola_samples.txt")
for line in reads_list:
    ind = line.strip("\n").split("\t")[0]
    fq1 = line.strip("\n").split("\t")[2]
    fq2 = line.strip("\n").split("\t")[3]
    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/STAR_index"
    index_log = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/STAR_index/DUM_star_index.log"
    RG = ind.replace("-","_")
    gwf.target_from_template(
        name="STAR_align_{RG}".format(RG=RG),
        template = STAR_mapping(align_path,index_path,index_log,fq1,fq2,RG)
    )

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
outname = "DUM_RNA_STAR"
logs = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/DUM_STAR_align_logs.txt"
files = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/DUM_STAR_align_bams.txt"
gwf.target_from_template(
    name="samtools_merge_bams_DUM",
    template = merge_bams(path,outname,logs,files)
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR.bam" 
bam_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam" 
log = "sort_DUM_RNA_bam.log"
gwf.target_from_template(
    name="samtools_sort_DUM",
    template = sort_bam(path,bam_in,bam_out,log)
    )

bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam" 
bai = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/DUM_RNA_STAR_sort.bam.bai" 
gwf.target_from_template(
    name="samtools_index_DUM",
    template = index_bam(bam,bai)
    )

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/STAR_index"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/TENT_hifi_hic_scaffolded_trim.fa"
log = "TENT_star_index.log"
species = "TENT"
gwf.target_from_template(
        name="STAR_index_{sp}".format(sp=species),
        template = STAR_index(path,ref,log)
    )
reads_list = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/TENT_30_RNA_meta.txt")
for line in reads_list:
    fq1 = line.strip("\n").split("\t")[0]
    fq2 = line.strip("\n").split("\t")[1]
    ind = fq1.split("/")[-2]
    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/STAR_index"
    index_log = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/STAR_index/TENT_star_index.log"
    RG = ind.replace("-","_")
    gwf.target_from_template(
        name="STAR_align_{RG}".format(RG=RG),
        template = STAR_mapping(align_path,index_path,index_log,fq1,fq2,RG)
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
outname = "TENT_RNA_STAR"
logs = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/TENT_STAR_align_logs.txt"
files = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/TENT_STAR_align_bams.txt"
gwf.target_from_template(
    name="samtools_merge_bams_TENT",
    template = merge_bams(path,outname,logs,files)
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/TENT_RNA_STAR.bam" 
bam_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/TENT_RNA_STAR_sort.bam" 
log = "sort_TENT_RNA_bam.log"
gwf.target_from_template(
    name="samtools_sort_TENT",
    template = sort_bam(path,bam_in,bam_out,log)
    )
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/TENT_RNA_STAR_sort.bam" 
bai = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/TENT_RNA_STAR_sort.bam.bai" 
gwf.target_from_template(
    name="samtools_index_TENT",
    template = index_bam(bam,bai)
    )
## SARA

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/STAR_index"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
log = "SARA_star_index.log"
species = "SARA"
gwf.target_from_template(
        name="STAR_index_{sp}".format(sp=species),
        template = STAR_index(path,ref,log)
    )
#reads_list = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/SARA_18_RNA_meta.txt")
#for line in reads_list:
#    fq1 = line.strip("\n").split("\t")[0]
#    fq2 = line.strip("\n").split("\t")[1]
#    ind = fq1.split("/")[-2]
#    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
#    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/STAR_index"
#    index_log = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/STAR_index/SARA_star_index.log"
#    RG = ind.replace("-","_")
#    gwf.target_from_template(
#        name="STAR_align_{RG}".format(RG=RG),
#        template = STAR_mapping(align_path,index_path,index_log,fq1,fq2,RG)
#    )
#path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
#outname = "SARA_RNA_STAR"
#logs = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/SARA_STAR_align_logs.txt"
#files = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/SARA_STAR_align_bams.txt"
#gwf.target_from_template(
#    name="samtools_merge_bams_SARA",
#    template = merge_bams(path,outname,logs,files)
#    )
#path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
#bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/SARA_RNA_STAR.bam" 
#bam_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/SARA_RNA_STAR_sort.bam" 
#log = "sort_SARA_RNA_bam.log"
#gwf.target_from_template(
#    name="samtools_sort_SARA",
#    template = sort_bam(path,bam_in,bam_out,log)
#    )
#bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/SARA_RNA_STAR_sort.bam" 
#bai = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/SARA_RNA_STAR_sort.bam.bai" 
#gwf.target_from_template(
#    name="samtools_index_SARA",
#    template = index_bam(bam,bai)
#    )


## LIN

path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/STAR_index"
ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/LIN_hifi_hic_scaffolded_trim.fa"
log = "LIN_star_index.log"
species = "LIN"
gwf.target_from_template(
        name="STAR_index_{sp}".format(sp=species),
        template = STAR_index(path,ref,log)
    )

reads_list = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/LIN_27_RNA_meta.txt")
for line in reads_list:
    fq1 = line.strip("\n").split("\t")[0]
    fq2 = line.strip("\n").split("\t")[1]
    ind = fq1.split("/")[-2]
    align_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
    index_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/STAR_index"
    index_log = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/STAR_index/LIN_star_index.log"
    RG = ind.replace("-","_")
    gwf.target_from_template(
        name="STAR_align_{RG}".format(RG=RG),
        template = STAR_mapping(align_path,index_path,index_log,fq1,fq2,RG)
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
outname = "LIN_RNA_STAR"
logs = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/LIN_STAR_align_logs.txt"
files = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/RNA_seq/star_map/LIN_STAR_align_bams.txt"
gwf.target_from_template(
    name="samtools_merge_bams_LIN",
    template = merge_bams(path,outname,logs,files)
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR"
bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/LIN_RNA_STAR.bam" 
bam_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/LIN_RNA_STAR_sort.bam" 
log = "sort_LIN_RNA_bam.log"
gwf.target_from_template(
    name="samtools_sort_LIN",
    template = sort_bam(path,bam_in,bam_out,log)
    )
bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/LIN_RNA_STAR_sort.bam" 
bai = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/RNA_analysis/STAR/LIN_RNA_STAR_sort.bam.bai" 
gwf.target_from_template(
    name="samtools_index_LIN",
    template = index_bam(bam,bai)
    )


