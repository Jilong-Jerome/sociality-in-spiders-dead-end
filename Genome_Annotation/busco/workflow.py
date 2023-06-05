from gwf import *

#### BUSCO run for AUGUSTUS gene model training

## Functions
def run_busco_protein(maker_fasta,path,outname):
    inputs=[maker_fasta]
    outputs=[path+"/"+outname+".log"]
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
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m proteins -c 16 -f > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
gwf = Workflow()
#LIN
maker_fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/lineatus/etp_mode/LIN_braker.protein.rename.fasta"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/lineatus/etp_mode/braker_busco/LIN"
outname = "LIN_braker_etp_mode"
gwf.target_from_template(
            name="busco_{name}".format(name=outname),
            template = run_busco_protein(maker_fasta,path,outname)
        )

#MIM
maker_fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/BRAKER_annotation/peptide/MIM.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/mimosarum/etp_mode/braker_busco/MIM"
outname = "MIM_braker_etp_mode"
gwf.target_from_template(
            name="busco_{name}".format(name=outname),
            template = run_busco_protein(maker_fasta,path,outname)
        )

#BI
maker_fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/BRAKER_annotation/peptide/BI.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/bicolor/etp_mode/braker_busco/BI"
outname = "BI_braker_etp_mode"
gwf.target_from_template(
            name="busco_{name}".format(name=outname),
            template = run_busco_protein(maker_fasta,path,outname)
        )

#SARA
maker_fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/BRAKER_annotation/peptide/SARA.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/sarasinorum/etp_mode/braker_busco/SARA"
outname = "SARA_braker_etp_mode"
gwf.target_from_template(
            name="busco_{name}".format(name=outname),
            template = run_busco_protein(maker_fasta,path,outname)
        )


