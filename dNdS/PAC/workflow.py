from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/consensus/logs"
def samtools_consensus_step(path,bam_in,fasta_out,id):
    inputs = [bam_in]
    outputs = [LOG_PATH + "bam2fa_{id}.DONE".format(id=id)]
    options = {
        "cores": 8,
        "memory": "32g",
        "walltime": "24:00:00",
        "account": "spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh

    echo jobinfo $SLURM_JOBID
    echo "Starting samtools consensus generation"

    mkdir -p {path}
    cd {path}
    conda activate samtools117
    samtools consensus -a --show-del yes --show-ins no -f fasta -o {fasta_out} {bam_in} -@ 8
    echo done > {log}
    """.format(log=outputs[0], bam_in=bam_in, path=path, fasta_out=fasta_out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def samtools_consensus_for_psmc(path,bam_in,fasta_out,id):
    inputs = [bam_in]
    outputs = [LOG_PATH + "bam2fa_{id}.DONE".format(id=id)]
    options = {
        "cores": 8,
        "memory": "32g",
        "walltime": "24:00:00",
        "account": "spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh

    echo jobinfo $SLURM_JOBID
    echo "Starting samtools consensus generation"

    mkdir -p {path}
    cd {path}
    conda activate samtools117
    samtools consensus -a --show-del yes --show-ins yes --min-MQ 30 --min-BQ 30  -d 10 -f fasta -o {fasta_out} {bam_in} -@ 8
    echo done > {log}
    """.format(log=outputs[0], bam_in=bam_in, path=path, fasta_out=fasta_out)
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


gwf = Workflow()

path =  "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/consensus/fasta"
id = "PAC_128"
bam_in = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/pacificus_consensus/pac128.bam"
fasta_out = id+"_consensus.fa"
gwf.target_from_template(
        name = "bam2fa_{id}".format(id=id),
        template = samtools_consensus_step(path,bam_in,fasta_out,id)
    )
id = "PAC_128_psmc"
fasta_out = id+"_consensus.fa"
gwf.target_from_template(
        name = "bam2fa_{id}".format(id=id),
        template = samtools_consensus_for_psmc(path,bam_in,fasta_out,id)
    )
