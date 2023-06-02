from workflow_templates import *

def run_blat(gwf,path,db,query,outname):
    log = "{}_blat.DONE".format(outname)
    gwf.target_from_template(
            name=outname,
            template = blat(path,db,query,outname,log)
            )
    return LOG_PATH+"/"+log

def run_filter_blat(gwf,path,outname):
    log = "{}_filter_blat.DONE".format(outname)
    gwf.target_from_template(
            name=outname+"_filter",
            template = filter_blat(path,outname,log)
            )
    return LOG_PATH+"/"+log

def run_augustus(gwf,path,fasta,cfg,hints,outname):
    log = "{}_augustus.DONE".format(outname)
    gwf.target_from_template(
            name=outname+"_augustus",
            template = augustus(path,fasta,cfg,hints,outname,log)
            )
    return LOG_PATH+"/"+log
