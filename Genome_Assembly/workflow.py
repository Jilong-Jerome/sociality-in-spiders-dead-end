from gwf import *

def hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out):
    inputs = [hifi_fq,hic_R1,hic_R2]
    outputs = [hifi_out+".hic.p_ctg.gfa"]
    options = {
              'cores':32,
              'memory':'256g',
              'walltime':'72:00:00',
              'account':"spider2"}
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate hifiasm
    echo "hifiasem on {hifi_reads}"
    echo jobinfo $SLURM_JOBID
    date
    hifiasm -o {hifi_out}.tmp -t32 --h1 {hic_1} --h2 {hic_2} {hifi_reads}
    mv {hifi_out}.tmp.hic.p_ctg.gfa {hifi_out}.hic.p_ctg.gfa
    date
    echo "hifiasm finish"
    """.format(hifi_reads=hifi_fq,hic_1=hic_R1,hic_2=hic_R2,hifi_out=hifi_out)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def gfa_fa(gfa):
    inputs = [gfa+".bp.p_ctg.gfa"]
    outputs = [gfa+".p_ctg.fasta"]
    options = {
              'cores':1,
              'memory':'1g',
              'walltime':'4:00:00',
              'account':"spider2"}
    code ="""
    "/^S/{print '>' $2 '\n' $3}"
    """
    spec="""
    awk {code} {gfa} > {fasta}.tmp
    mv {fasta}.tmp {fasta}
    """.format(gfa=inputs[0],fasta=outputs[0],code=code)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


gwf = Workflow()
hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/SARA/SARA_hifi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/SARA_hifi_hic"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_2_finish/SARA_hic_R1.fastq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_2_finish/SARA_hic_R2.fastq"
gwf.target_from_template(
    name = "SARA_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)


hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/DUM/DUM_HiFi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/DUM/DUM_hifi"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sdum_OmniC_lib2/211104_X600_FCHH2WHCCX2_L4_CHKPE85221100091_2.fq"
gwf.target_from_template(
    name = "DUM_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)


hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/TENT/TENT_HiFi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/TENT/TENT_hifi"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_1_finish/211203_X517_FCHH3HCCCX2_L4_CHKPE85221110117_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_1_finish/211203_X517_FCHH3HCCCX2_L4_CHKPE85221110117_2.fq"
gwf.target_from_template(
    name = "TENT_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)

hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/LIN/LIN_HiFi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/LIN/LIN_hifi"
hic_R1  ="/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_3_finish/211203_X517_FCHH3HCCCX2_L6_CHKPE85221110119_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/Sample_3_finish/211203_X517_FCHH3HCCCX2_L6_CHKPE85221110119_2.fq"
gwf.target_from_template(
    name = "LIN_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)

hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/MIM/MIM_hifi.fastq"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/MIM/MIM_hifi"
hic_R1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/BI_MIM/Mim/Mim_L1_1.fq"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/BI_MIM/Mim/Mim_L1_2.fq"
gwf.target_from_template(
    name = "MIM_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)

hifi_fq = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Pacbio_Hifi/BICO/S_Bi265_hifi_merged.fastq.gz"
hifi_out = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/hifi_asm/BI/BI_hifi"
hic_R1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/BI_MIM/B4/B4_L1_1.fq.gz"
hic_R2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/Hi-C/BI_MIM/B4/B4_L1_2.fq.gz"
gwf.target_from_template(
    name = "BI_hifi_hic_asm",
    template = hifi_asm(hifi_fq,hic_R1,hic_R2,hifi_out)
)
