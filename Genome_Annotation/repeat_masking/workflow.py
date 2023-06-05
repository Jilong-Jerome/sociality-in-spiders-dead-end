from gwf import *

def combine_fq_gzip(filelist,outpath,outname):
    """Combine fastq file in gzipped formmat from a list of files
    filelist: a list of files with absolute directory path
    outpath: folder storing the output combined fastq file
    outname: New name given to the new combined file
    """
    inputs=filelist
    outputs = ["{outpath}/{outname}.fq".format(outpath=outpath,outname=outname)]
    options = {
              'cores':1,
              'memory':'8g',
              'walltime':'24:00:00',
              'account':'spider2'
              }
    file_string = ""
    for f in filelist:
        file_string = file_string+"\t"+f
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start combining fastq files (gzipped)"
    echo jobinfo $SLURM_JOBID
    date
    zcat {file_string} > {outpath}/{outname}.fq.tmp
    mv {outpath}/{outname}.fq.tmp {outpath}/{outname}.fq
    echo "finished combining fastq file"
    date
    jobinfo $SLURM_JOBID
    """.format(file_string=file_string,outpath=outpath,outname=outname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def trinity_assembly(path,read1,read2,outname):
    inputs = [read1,read2]
    outputs = ["{path}/{outname}.log".format(path=path,outname=outname)]
    options = {
              'cores':16,
              'memory':'256g',
              'walltime':'72:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start trinity transcriptome assembly"
    echo jobinfo $SLURM_JOBID
    date
    Trinity --seqType fq --max_memory 240G --left {read1} --right {read2} --CPU 16 --output "{path}/{outname}" > {path}/{outname}.log.tmp
    mv {path}/{outname}.log.tmp {path}/{outname}.log
    echo "finished trinity assembling"
    date
    jobinfo $SLURM_JOBID
    """.format(read1=read1,read2=read2,outname=outname,path=path)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def maker_pipe(path,species):
    inputs = [path+"/maker_bopts.ctl",
              path+"/maker_exe.ctl",
              path+"/maker_opts.ctl"]
    outputs = [
              path+"/maker_{species}.log".format(species=species)]
    options = {
              'cores':32,
              'memory':'100g',
              'walltime':'167:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start maker pipe"
    echo jobinfo $SLURM_JOBID
    date
    cd {path}
    maker > {path}/maker_{species}.log.tmp
    mv {path}/maker_{species}.log.tmp {path}/maker_{species}.log
    echo "finished maker"
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def repeat_model(genome,out,path):
    inputs = [genome]
    outputs = [path+"/{out}.log".format(out=out)]
    options = {
              'cores':16,
              'memory':'128g',
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
    BuildDatabase -name {out} -engine ncbi {genome}
    RepeatModeler -pa 16 -engine ncbi -database {out} 2>&1 | tee {log}.tmp
    mv {log}.tmp {log}
    echo "finished RepeatModeler"
    date
    jobinfo $SLURM_JOBID 
    """.format(out=out,genome = genome,log = outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def repeat_masker_default(path,fasta,species):
    inputs = [fasta]
    outputs = [path+"/RM_{species}.log".format(species=species)]
    options = {
              'cores':16,
              'memory':'128g',
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
    RepeatMasker -e ncbi -pa 16 -dir {path} --species Eukaryota {fasta} > {path}/RM_{species}.log.tmp
    mv {path}/RM_{species}.log.tmp {path}/RM_{species}.log
    date
    echo jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def repeat_masker_simple(path,fasta,species):
    inputs = [fasta]
    outputs = [path+"/RM_{species}.log".format(species=species)]
    options = {
              'cores':16,
              'memory':'128g',
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
    RepeatMasker -e ncbi -pa 16 -dir {path} --species simple {fasta} > {path}/RM_{species}.log.tmp
    mv {path}/RM_{species}.log.tmp {path}/RM_{species}.log
    date
    echo jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,species=species)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def repeat_masker_lib(path,fasta,lib,species):
    inputs = [fasta]
    outputs = [path+"/RM_{species}.log".format(species=species)]
    options = {
              'cores':16,
              'memory':'128g',
              'walltime':'120:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start RepeatMasker"
    echo jobinfo $SLURM_JOBID
    date
    RepeatMasker -e ncbi -pa 16 -xsmall -dir {path} -lib {lib} {fasta} > {path}/RM_{species}.log.tmp
    mv {path}/RM_{species}.log.tmp {path}/RM_{species}.log
    date
    echo jobinfo $SLURM_JOBID
    """.format(path=path,fasta=fasta,species=species,lib=lib)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def out2gff3(path,out):
    inputs = [path+"/"+out]
    outputs = [path+"/"+out+".gff3"]
    options = {
              'cores':4,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start Making gff3"
    echo jobinfo $SLURM_JOBID
    rmOutToGFF3custom -o {infile} > {output}.tmp
    mv {output}.tmp {output}
    echo "create GFF3 finished"
    """.format(infile=inputs[0],output=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
     
def isolate_complex_gff(path,out):
    inputs = [path+"/"+out+".gff3"]
    outputs = [path+"/"+out+".complex.reformat.gff3"] 
    options = {
              'cores':1,
              'memory':'8g',
              'walltime':'6:00:00',
              'account':'spider2'
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start isolating complex repeats from gff3"
    echo jobinfo $SLURM_JOBID
    grep -v -e "Satellite" -e "Low_complexity" -e "Simple_repeat" {infile} > {output}.tmp 
    bash /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/annotate/reform_maker_gff3.sh {output}.tmp {output}.reform.tmp 
    mv {output}.reform.tmp {output}
    rm {output}.tmp
    echo "complex GFF3 finished"
    """.format(infile=inputs[0],output=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def bedtools_maskfasta(path,softmasked,complex_gff,new_name):
    inputs = [softmasked]
    outputs = [path+"/"+new_name+".split.masked.fa"] 
    options = {
              'cores':1,
              'memory':'8g',
              'walltime':'12:00:00',
              'account':'spider2'
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate annotate
    echo "start hard masking complex repeats based on soft masked assembly"
    echo jobinfo $SLURM_JOBID
    bedtools maskfasta -fi {softmasked} -bed {complex_gff} -fo {output}.tmp
    mv {output}.tmp {output}
    echo "complex hard masked"
    """.format(softmasked=softmasked,complex_gff=complex_gff,output=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

gwf = Workflow()


path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/maker/dumicola"
species = "dumicola"
genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/DUM_hifi_hic_scaffolded_trim.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_modeler"
outname = "dumicola"
gwf.target_from_template(
    name="repeat_modeler_{}".format(outname),
    template = repeat_model(
              genome = genome,
              out = outname,
              path = path
              )
    )

genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_modeler"
outname = "Sarasinorum"
gwf.target_from_template(
    name="repeat_modeler_{}".format(outname),
    template = repeat_model(
              genome = genome,
              out = outname,
              path = path
              )
    )

genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/TENT_hifi_hic_scaffolded_trim.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_modeler"
outname = "Tentorricola"
gwf.target_from_template(
    name="repeat_modeler_{}".format(outname),
    template = repeat_model(
              genome = genome,
              out = outname,
              path = path
              )
    )

genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/LIN_hifi_hic_scaffolded_trim.fa"
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_modeler"
outname = "Lineatus"
gwf.target_from_template(
    name="repeat_modeler_{}".format(outname),
    template = repeat_model(
              genome = genome,
              out = outname,
              path = path
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/DUM/DUM_hifi_hic_scaffolded_trim.fa"
species = "dumicola"
gwf.target_from_template(
    name="repeat_masker_default_{}".format(species),
    template = repeat_masker_default(
              path = path,
              fasta = fasta,
              species = species
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/annotate/Dumicola-families.fa"
species = "dumicola_selflib"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/repbase_arthropod"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/repbase/repbase_arthropoda.fa"
species = "dumicola_repbase_arthropod"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine/Dumicola_Arthropoda.fa"
species = "dumicola_combine"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/TENT/TENT_hifi_hic_scaffolded_trim.fa"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/annotate/Tentorricola-families.fa"
species = "tentoriicola_selflib"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase/combine"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase/combine/Tentorricola_Arthropoda.fa"
species = "tentorricola_combine"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/SARA/RMdatabase"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/SARA/SARA_hifi_hic_scaffolded_trim.fa"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/annotate/Sarasinorum-families.fa"
species = "sarasinorum_selflib"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/SARA/RMdatabase/combine"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/SARA/RMdatabase/combine/Sarasinorum_Arthropoda.fa"
species = "sarasinorum_combine"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase"
fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/LIN/LIN_hifi_hic_scaffolded_trim.fa"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/annotate/Lineatus-families.fa"
species = "lineatus_selflib"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase/combine"
lib = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase/combine/Lineatus_Arthropoda.fa"
species = "lineatus_combine"
gwf.target_from_template(
    name="repeat_masker_lib_{}".format(species),
    template = repeat_masker_lib(
              path = path,
              fasta = fasta,
              species = species,
              lib = lib
              )
    )
path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/DUM/RMdatabase/combine"
out = "DUM_hifi_hic_scaffolded_trim.fa.out"
species = "dumicola"
gwf.target_from_template(
    name="create_gff3_{}".format(species),
    template = out2gff3(
              path = path,
              out = out
              )
)
gwf.target_from_template(
    name="isolate_complex_gff_{}".format(species),
    template = isolate_complex_gff(
              path = path,
              out = out
              )
)
path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/TENT/RMdatabase"
out = "TENT_hifi_hic_scaffolded_trim.fa.out"
species = "tentorricola"
gwf.target_from_template(
    name="create_gff3_{}".format(species),
    template = out2gff3(
              path = path,
              out = out
              )
)
gwf.target_from_template(
    name="isolate_complex_gff_{}".format(species),
    template = isolate_complex_gff(
              path = path,
              out = out
              )
)
path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/SARA/RMdatabase"
out = "SARA_hifi_hic_scaffolded_trim.fa.out"
species = "sarasinorum"
gwf.target_from_template(
    name="create_gff3_{}".format(species),
    template = out2gff3(
              path = path,
              out = out
              )
)
gwf.target_from_template(
    name="isolate_complex_gff_{}".format(species),
    template = isolate_complex_gff(
              path = path,
              out = out
              )
)
path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/LIN/RMdatabase"
out = "LIN_hifi_hic_scaffolded_trim.fa.out"
species = "lineatus"
gwf.target_from_template(
    name="create_gff3_{}".format(species),
    template = out2gff3(
              path = path,
              out = out
              )
)
gwf.target_from_template(
    name="isolate_complex_gff_{}".format(species),
    template = isolate_complex_gff(
              path = path,
              out = out
              )
)

