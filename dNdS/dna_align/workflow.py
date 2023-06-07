from gwf import *
#####WORKFLOW to ALIGN and CHECK BAM COVERAGE
####FUNCTIONS
## faidx
def samtool_faidx(fasta,fai):
    """samtools for creating fasta index file"""
    inputs = [fasta]
    outputs = [fai]
    options = {
              'cores': 1,
              'memory': '8g',
              'walltime':'2:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools for indexing assembly fasta"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/samtools-1.12/samtools faidx {ifa} -o {faidx}.tmp
    mv {faidx}.tmp {faidx}
    echo "finished"
    date
    jobinfo $SLURM_JOBID
    """.format(ifa = inputs[0],faidx=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
## bwa2mem indexing
def bwa2_index(ref_genome,path,basename):
    """ Template for indexing a genome with bwa2-mem"""
    inputs = [ref_genome,
              ref_genome+".fai"]
    outputs = ["{path}/{basename}_bwa2index.log".format(path=path,basename=basename)] 
    options = {
               'cores': 12,
               'memory': '256g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo jobinfo $SLURM_JOBID
    echo "start indexing bwa2"
    date
    mkdir -p {path}
    /home/jilong/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index -p {path}/{basename} {ref_genome} > {path}/{basename}_bwa2index.log.tmp
    mv {path}/{basename}_bwa2index.log.tmp {path}/{basename}_bwa2index.log
    echo "finish indexing bwa2"
    date
    jobinfo $SLURM_JOBID
    """.format(ref_genome=ref_genome,path=path,basename=basename)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

## merge fq
def combine_fq(combine_list,path,outname):
    """ Combine reads from several Lane to a single file"""
    inputs = combine_list
    catfiles=""
    for i in inputs:
        catfiles=catfiles+str(i)+" "
    outputs = ['{path}/{outname}.fq.gz'.format(path=path,outname=outname)]
    options = {
              'cores': 1,
              'memory': "2g",
              'walltime':"08:00:00",
              'account':"spider2"
              }
    spec = """
    echo jobinfo $SLURM_JOBID
    echo "start combining fq"
    date
    cat {catfiles} > {combined}.tmp
    mv {combined}.tmp {combined}
    echo "finished combination"
    date
    jobinfo $SLURM_JOBID
    """.format(catfiles=catfiles,combined=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

##  bwa2 align
def bwa2_align(index_prefix,index_log,path,indname,read1,read2):
    """ Template for aligning paired read to a genome with bwa2-mem"""
    inputs = [ index_log,
               read1,
               read2]
    outputs = ['{path}/{indname}.bam'.format(path=path,indname=indname)]
    options = {
               'cores': 24,
               'memory': '80g',
               'walltime':"24:00:00",
               'account':"spider2"
    }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo jobinfo $SLURM_JOBID
    echo "start aligning bwa-mem2"
    date
    mkdir -p {path}
    /home/jilong/software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 24 {index_prefix} {read1} {read2} | /home/jilong/software/samtools-1.12/samtools view -Sb > {out}.tmp
    mv {out}.tmp {out}
    echo "finish aligning with bwa-mem2"
    date
    jobinfo $SLURM_JOBID
    """.format(index_prefix=index_prefix,path=path,read1=read1,read2=read2,out=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

## bam filters
def picard_addRG(path,indname):
    """GATK4.0 preprocessing, read in bam file need RGs"""
    inputs = [path+'/{}.bam'.format(indname)]
    outputs = [path+'/{}_RG.bam'.format(indname)]
    options = {
              'cores': 8,
              'memory': '16g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start gatk markduplicatesSpark"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/java/jre1.8.0_291/bin/java -jar /home/jilong/software/picard.jar AddOrReplaceReadGroups --I {ibam} --O {obam}.tmp -RGID {sample} -RGPU unknown -RGSM {sample} -RGPL illumina -RGLB lib0
    mv {obam}.tmp {obam}
    echo "finished gatk markduplicatesSpark"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obam=outputs[0],sample=indname)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_fixmate(path,indname):
    """Samtools fixtmate for remove duplicates"""
    inputs = [path+'/{}_RG.bam'.format(indname)]
    outputs = [path+'/{}_fixmate.bam'.format(indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools fixmate, remove secondary and unmapped reads"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/samtools-1.12/samtools fixmate -rm {ibam} {obam}.tmp -@ 16
    mv {obam}.tmp {obam}
    echo "finished samtools fixmate"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obam=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_sort(path,indname):
    """Samtools sort bam file"""
    inputs = [path+'/{}_fixmate.bam'.format(indname)]
    outputs = [path+'/{}_sorted.bam'.format(indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools sort bam"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/samtools-1.12/samtools sort {ibam} -o {obam}.tmp -@ 16
    mv {obam}.tmp {obam}
    echo "finished samtools sort"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obam=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_markdup(path,indname):
    """Samtools markdup for remove duplicates"""
    inputs = [path+'/{}_sorted.bam'.format(indname)]
    outputs = [path+'/{}_dup.bam'.format(indname),
               path+'/{}_dup.stat'.format(indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools markdup to remove duplicate reads"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/samtools-1.12/samtools markdup -r -f {stat}.tmp -s  {ibam} {obam}.tmp -@ 16
    mv {obam}.tmp {obam}
    mv {stat}.tmp {stat}
    echo "finished samtools markdup"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obam=outputs[0],stat=outputs[1])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def samtool_mapped(path,indname):
    """Samtools take only the mapped reads"""
    inputs = [path+'/{}_dup.bam'.format(indname)]
    outputs = [path+'/{}_final.bam'.format(indname)]
    options = {
              'cores': 16,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start samtools to remove unmapped reads"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/samtools-1.12/samtools view -@ 16 -bq 60  -f 0x2 -F 0x4  {ibam} > {obam}.tmp
    mv {obam}.tmp {obam}
    echo "finished samtools get mapped reads"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obam=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def picard_bai(bam):
    """GATK4.0 preprocessing, index bam files"""
    inputs = [bam]
    outputs = ['{}.bai'.format(bam)]
    options = {
              'cores': 12,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec = """    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate spider
    echo "start picar build bam index"
    echo jobinfo $SLURM_JOBID
    date
    /home/jilong/software/java/jre1.8.0_291/bin/java -jar /home/jilong/software/picard.jar BuildBamIndex --INPUT {ibam} --OUTPUT {obai}.tmp
    mv {obai}.tmp {obai}
    echo "finished picard index bam"
    date
    jobinfo $SLURM_JOBID
    """.format(ibam=inputs[0],obai=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def get_depth_of_scaffold(bam,scaffold,path,ind):
    obam = path+"/"+ind+"_"+scaffold+".bam"
    inputs = [bam]
    outputs = [path+"/"+ind+"_"+scaffold+".depth"]
    options = {
              'cores': 1,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bamtools
    echo "start depth calling"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    samtools view -b {bam} {scaffold} > {obam}.tmp
    mv {obam}.tmp {obam}
    bamtools coverage -in {obam} -out {out}.tmp
    mv {out}.tmp {out}
    rm {obam}
    date
    """.format(path=path,bam=bam,obam=obam,scaffold=scaffold,out=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


def depth_of_scaffold(bam,ind,scaffold,region,path):
    inputs = [bam]
    outputs = [path+"/"+ind+"_"+scaffold+".depth"]
    options = {
              'cores': 1,
              'memory': '32g',
              'walltime':'12:00:00',
              'account':"spider2"
              }
    spec="""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bamtools
    echo "start depth calling"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    samtools coverage {bam} -r {region} > {out}.tmp
    mv {out}.tmp {out}
    date
    """.format(path=path,bam=bam,out=outputs[0],region=region)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

#Workflow
gwf = Workflow()


###prepare meta data
reference_dict={}
refs = open("genome_reference.txt")
for ref in refs:
    species = ref.strip("\n").split("\t")[0]
    genome = ref.strip("\n").split("\t")[1]
    reference_dict[species]=genome
species_reads = {}
for species in ["SARA","LIN"]:
    reads = open("{}_reads.txt".format(species))
    reads_dict={}
    for read in reads:
        read_pair = read.strip("\n").split("/")[-1].split(".")[0][-1]
        family = read.strip("\n").split("/")[-3]
        ind = read.strip("\n").split("/")[-2].split("-")[-1]
        if ind == "F":
            suffix = "female"
        elif ind == "M":
            suffix = "male"
        else:
            suffix = "offspring"
        newname = species+"_"+family+"_"+ind+"_"+suffix
        if newname in reads_dict:
            if read_pair == "1":
                reads_dict[newname][0].append(read.strip("\n"))
            if read_pair == "2":
                reads_dict[newname][1].append(read.strip("\n"))
        else:
            reads_dict[newname]=[[],[]]
            if read_pair == "1":
                reads_dict[newname][0].append(read.strip("\n"))
            if read_pair == "2":
                reads_dict[newname][1].append(read.strip("\n"))
    species_reads[species]=reads_dict   
ind_name_list=[]
inds = open("ind_list.txt")
for ind in inds:
    indname = ind.split(".")[0][0:-3]
    ind_name_list.append(indname)
## merging fq
path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations/ind_fq_merged"
for species in ["SARA","LIN"]:
    for ind in species_reads[species]:    
        R1_combine = species_reads[species][ind][0]
        R2_combine = species_reads[species][ind][1]
        gwf.target_from_template(
            name = "merge_R1_{}".format(ind),
            template = combine_fq(
                combine_list = sorted(R1_combine),
                path = path,
                outname = ind+"_R1"
            )
        )
        gwf.target_from_template(
            name = "merge_R2_{}".format(ind),
            template = combine_fq(
                combine_list = sorted(R2_combine),
                path = path,
                outname = ind+"_R2"
            )
        )
### Indexing of the genome
## Create ".fai" file
for species in ["BI","DUM","TENT","SARA","LIN","PAC"]:
    fasta = reference_dict[species]
    fai = fasta+".fai"
    #print("Samtools indexing of {species}".format(species=species))
    #print("Input fasta file:{fasta}".format(fasta=fasta))
    #print("Output fai file:{fai}".format(fai=fai))
    #print("\n")
    gwf.target_from_template(
        name="fasta_index_{}".format(species),
        template = samtool_faidx(
               fasta=fasta,
               fai=fai
               )
    )
## Create bwa index
for species in ["BI","DUM","TENT","SARA","LIN","MIM","PAC"]:
    fasta = reference_dict[species]
    basename = species
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_index/{}".format(species) 
    gwf.target_from_template(
        name="bwa2_index_{}".format(species),
        template = bwa2_index(
               ref_genome = fasta,
               path = path,
               basename = basename
               )
    )
### bwa mem 2 align
species_trans={"D":"DUM",
               "T":"TENT",
               "M":"MIM",
               "SARA":"SARA",
               "LIN":"LIN",
               "BI":"BI",
               "PAC":"PAC"}
for ind in ind_name_list:
    species = ind.split("_")[0]
    species_code = species_trans[species]
    if species_code != "NONE":
        index_prefix = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_index/{code}/{code}".format(code = species_code)
        index_log = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_index/{code}/{code}_bwa2index.log".format(code = species_code)
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/{code}/".format(code = species_code)
        indname = ind
        read1 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations/ind_fq_merged/{ind}_R1.fq.gz".format(ind=ind)
        read2 = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/mutations/ind_fq_merged/{ind}_R2.fq.gz".format(ind=ind)
###bam filter
        gwf.target_from_template(
            name = "align_bwa2_{}".format(indname),
            template = bwa2_align(
                index_prefix = index_prefix,
                index_log = index_log,
                path = path,
                indname = indname,
                read1 = read1,
                read2 = read2
                    )
        )
        gwf.target_from_template(
            name = "picard_addRG_{}".format(indname),
            template = picard_addRG(path,indname)
        )
        gwf.target_from_template(
            name = "samtools_fixmate_{}".format(indname),
            template = samtool_fixmate(path,indname)
        )
        gwf.target_from_template(
            name = "samtools_sort_{}".format(indname),
            template = samtool_sort(path,indname)
        )
        gwf.target_from_template(
            name = "samtools_markdup_{}".format(indname),
            template = samtool_markdup(path,indname)
        )
        gwf.target_from_template(
            name = "samtools_mapped_MQ_{}".format(indname),
            template = samtool_mapped(path,indname)
        )
        bam = path+"/{}_final.bam".format(indname)
        gwf.target_from_template(
            name = "picard_bai_{}".format(indname),
            template = picard_bai(bam)
        )
## bam coverage call
#for ind in ind_name_list:
#    #print(ind)
#    species = ind.split("_")[0]
#    species_code = species_trans[species]
#    sex = ind.split("_")[-1].strip("\n")
#    #print(sex)
#    if sex == "male" and species_code != "MIM":
#        IDs_file = open("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/repeat_masker/{species_code}/RMdatabase/combine/id.txt".format(species_code=species_code))
#        for ID in IDs_file:
#            scaffold = ID.strip("\n")
#            #print(scaffold)
#            path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/depth/{species_code}/{ind}".format(species_code=species_code,ind=ind)
#            bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/{species_code}/{ind}_final.bam".format(species_code=species_code,ind=ind)
#            gwf.target_from_template(
#                name = "call_depth_{ind}_{scaffold}".format(ind=ind,scaffold=scaffold),
#                template = get_depth_of_scaffold(bam,scaffold,path,ind)
#            )

region_dict={}
for code in ["DUM","TENT","LIN","SARA"]:
    region_dict[code]=[]
    regions = open("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{code}/{code}_hifi_hic_scaffolded_trim.fa.fai".format(code=code))
    for line in regions:
        chrom = line.split("\t")[0]
        end = line.split("\t")[1]
        region = chrom+":1-"+end
        region_dict[code].append(region)
for ind in ind_name_list:
    #print(ind)
    species = ind.split("_")[0]
    species_code = species_trans[species]
    sex = ind.split("_")[-1].strip("\n")
    #print(sex)
    if sex == "male" and species_code != "MIM":
        path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/samtools_coverage/{species_code}/{ind}".format(ind=ind,species_code=species_code)
        bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/{species_code}/{ind}_final.bam".format(species_code=species_code,ind=ind)
        regions = region_dict[species_code]
        for region in regions:
            scaffold = region.split(":")[0]
            gwf.target_from_template(
                name = "samtools_coverage_{ind}_{scaffold}".format(ind=ind,scaffold=scaffold),         
                template = depth_of_scaffold(bam,ind,scaffold,region,path)
            )
for ind in ind_name_list:
    #print(ind)
    species = ind.split("_")[0]
    species_code = species_trans[species]
    sex = ind.split("_")[-1].strip("\n")
    #print(sex)
    if sex == "female" and species_code != "MIM":
        path="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/samtools_coverage/{species_code}/{ind}".format(ind=ind,species_code=species_code)
        bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/chromosome_X/bwa2_align/{species_code}/{ind}_final.bam".format(species_code=species_code,ind=ind)
        regions = region_dict[species_code]
        for region in regions:
            scaffold = region.split(":")[0]
            gwf.target_from_template(
                name = "samtools_coverage_{ind}_{scaffold}".format(ind=ind,scaffold=scaffold),         
                template = depth_of_scaffold(bam,ind,scaffold,region,path)
            )



