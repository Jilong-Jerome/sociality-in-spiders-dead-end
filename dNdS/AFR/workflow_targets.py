from gwf import *
from workflow_templates import *
def run_trinity_single(path,read1,read2,outname,gwf):
    gwf.target_from_template(
    name="trinity_{outname}".format(outname=outname),
    template = trinity_assembly(
              path = path,
              read1=read1,
              read2=read2,
              outname = outname
              )
    )
    return LOG_PATH+"/"+outname+".DONE"
