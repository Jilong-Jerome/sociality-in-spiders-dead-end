from gwf import *

#### BUSCO run for AUGUSTUS gene model training

## Functions
##Put together training sequences
def put_train_together(path,noseq_gff,fasta,fai,name):
    inputs=[noseq_gff,fasta,fai]
    outputs=[path+"/sbatch.log"]
    options={
        'cores':1,
        'memory':'4g',
        'walltime':'12:00:00',
        'account':'spider2'
    }
    spec="""
    echo jobinfo $SLURM_JOBID
    sbatch gather_train.sh {noseq_gff} {fai} {fasta} {path} {name}> {path}/sbatch.log 
    """.format(noseq_gff=noseq_gff,fai=fai,fasta=fasta,path=path,name=name)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def run_busco(maker_fasta,path,outname):
    inputs=[maker_fasta]
    outputs=[path+"/"+outname+".log"]
    options={
        'cores':32,
        'memory':'16g',
        'walltime':'12:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo "start busco for augustus training"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m genome -c 32 --long -f --augustus_species fly  --augustus_parameters='--progress=true' > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def run_busco_test(maker_fasta,path,outname):
    inputs=[maker_fasta]
    outputs=[path+"/"+outname+".log"]
    options={
        'cores':32,
        'memory':'16g',
        'walltime':'12:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo "start busco for augustus training"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m genome -c 32 --long -f --augustus_species parasteatoda --augustus_parameters='--progress=true' > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def run_busco_round(maker_fasta,path,outname,sp):
    inputs=[maker_fasta]
    outputs=[path+"/"+outname+".log"]
    options={
        'cores':32,
        'memory':'16g',
        'walltime':'12:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo "start busco for augustus training"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m genome -c 32 --long -f --augustus_species {sp}  --augustus_parameters='--progress=true' > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname,sp=sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def run_busco_check(maker_fasta,path,outname,sp):
    inputs=[maker_fasta]
    outputs=[path+"/"+outname+".log"]
    options={
        'cores':32,
        'memory':'64g',
        'walltime':'24:00:00',
        'account':'spider2'
    }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base
    echo "start busco for checking completeness"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m transcriptome -c 32 -f --augustus_species {sp}  --augustus_parameters='--progress=true' > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname,sp=sp)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def run_busco_trinity(maker_fasta,path,outname):
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
    busco -i {maker_fasta} -o {outname} -l arthropoda_odb10 -m transcriptome -c 16 -f > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    """.format(maker_fasta=maker_fasta,path=path,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
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
#TENT
maker_fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/tentoriicola/etp_mode/TENT_braker.rename.fasta"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/tentoriicola/etp_mode/braker_busco/TENT"
outname = "TENT_braker_etp_mode"
gwf.target_from_template(
            name="busco_{name}".format(name=outname),
            template = run_busco_trinity(maker_fasta,path,outname)
        )
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


