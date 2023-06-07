from workflow_templates import *
from gwf import *

def run_bam2vcf(gwf,path,bam,ref,name):
    gwf.target_from_template(
        name = "bam2vcf_{name}".format(name=name),
        template = bam2vcf(path,bam,ref,name)
    )
    return LOG_PATH+"/bam2vcf_{name}.DONE".format(name = name)

def run_noindel(gwf,path,name):
    vcfin = path+"/{name}.vcf".format(name=name)
    vcfout = path+"/{name}_noindel".format(name=name)
    gwf.target_from_template(
        name = "noindel_{name}".format(name=name),
        template = vcf_noindel(path,vcfin,vcfout,name)
    )
    return LOG_PATH+"/noindel_{name}.DONE".format(name = name)

def run_altref(gwf,path,vcfin,refin,name):
    refout = path+"/{name}.fa".format(name=name)
    gwf.target_from_template(
        name = "vcf2ref_{name}".format(name=name),
        template = vcf2ref(path,vcfin,refin,refout,name)
    )
    return LOG_PATH+"/altref_{name}.DONE".format(name = name)

def run_agat(gwf,path,refin,gffin,name):
    faout = "{name}_cds.fa".format(name=name)
    gwf.target_from_template(
        name = "getcds_{name}".format(name=name),
        template = agat_retrieve_coding_fasta(path,refin,gffin,faout,name)
    )
    return LOG_PATH+"/getcds_{name}.DONE".format(name = name)

def run_concat(gwf,path,fain,faout,name):
    gwf.target_from_template(
        name = "catcds_{name}".format(name=name),
        template = concat_fasta(path,fain,faout,name)
    )
    return LOG_PATH+"/catcds_{name}.DONE".format(name = name)
