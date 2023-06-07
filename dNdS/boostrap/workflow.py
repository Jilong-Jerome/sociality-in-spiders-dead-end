from gwf import *
LOG_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/logs"
ALN_PATH = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/single_ortho"

def random_sample(path, target, n, out):
    inputs = [target]
    outputs = [LOG_PATH+"/sample_{out}.DONE".format(out=out),
              path+"/{out}_id.txt".format(out=out)]
    options = {
              'cores':1,
              'memory':'100m',
              'walltime':'00:05:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate paml
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    shuf -n {n} {target} |cut -f1 >  {out}_id.txt 
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(out=out,path=path,target=target,log=outputs[0],n=n)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def cat_og_phy_filter(path,id_list,out):
    #cat_string  = "cat "
    cat_string  = ""
    ids = open(path+"/"+id_list+"_id.txt")
    for og_id in ids:
        phy_file = ALN_PATH+"/{og_id}_NT_fix.fasta".format(og_id=og_id.strip("\n"))
        cat_string = cat_string + phy_file + " "
    #print(cat_string)
    inputs = [LOG_PATH+"/sample_{id_list}.DONE".format(id_list=id_list),
              path+"/{out}_id.txt".format(out=out.replace(".fa",""))]
    outputs = [LOG_PATH+"/cat_{id_list}.DONE".format(id_list=id_list)]
    options = {
              'cores':1,
              'memory':'4g',
              'walltime':'02:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate goalign
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    # {cat_string} > {out}.tmp
    goalign concat -i {cat_string} > {out}.cat
    # rm {out}.tmp 
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/align_filter.py {out}.cat {out}.region
    echo raw_region_done
    conda activate gwf
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/find_region.R {out}.region {out}.region.all
    echo region_summary_done
    conda activate gwf
    Rscript /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/filter_region.R {out}.region.all {out}.region.filtered
    echo region_retrieve_done
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/combine_OG/align_filter_region_fa.py {out}.cat {out}.region.filtered.tsv {out}.filtered.phy
    echo retreive_region_align_done
    conda activate goalign
    goalign concat -i {out}.filtered.phy -p > {out}.filtered.concat.phy
    conda activate clustalo
    trimal -in {out}.filtered.concat.phy -phylip3.2 -out {out}.filtered.concat.paml.phy
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(out=out,path=path,cat_string=cat_string,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def paml_filter(alignment,path,codeml_path,mode,id_list,reform_script):
    inputs = [LOG_PATH+"/cat_{id_list}.DONE".format(id_list=id_list)]
    outputs = [
              LOG_PATH+"/paml_{mode}_{id_list}.DONE".format(mode=mode,id_list=id_list)
              ]
    options = {
              'cores':1,
              'memory':'5g',
              'walltime':'04:00:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate prank
    echo "start prank aligning single ortholog groups fasta"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    cp {codeml_path}/codeml_{mode}.ctl codeml.ctl
    cp {codeml_path}/tree_{mode}.txt tree.txt
    #prank -convert -d={alignment} -f=phylips -o=seq.txt -keep
    #mv seq.txt.phy seq.txt
    ln -f -s {group}.fa.filtered.concat.paml.phy seq.txt
    conda activate paml
    codeml     
    echo {group} > {group}.name
    paste {group}.name rst1 > {group}.res
    rm {group}.name
    tail -n 10 results.txt > results.raw
    conda activate ete3
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/combine_OG/boostrap/paml2tab.py results.raw {group}.tab {group}
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/boostrap/pair_social.py {group}.tab {group}_social.tab
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,alignment=alignment,reform_script=reform_script,codeml_path=codeml_path,mode=mode,group=id_list,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)


gwf = Workflow()
for i in range(500):
    n = 500
    out = "auto_{n}_{j}".format(n=n,j=i+1)
    target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/auto_og_pass.txt"
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}".format(out=out)
    gwf.target_from_template(
            name = "sample_{out}".format(out=out),
            template = random_sample(path, target, n, out)
        )
    gwf.target_from_template(
        name = "cat_{out}".format(out=out),
        template = cat_og_phy_filter(path,out,out+".fa")
    )
    codeml_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/codeml"
    mode = "ofree"
    reform_script = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/reform_results.py"
    alignment = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}/{out}.fa".format(out=out)
    gwf.target_from_template(
            name = "paml_{mode}_{id_list}".format(mode=mode,id_list=out),
            template = paml_filter(alignment,path,codeml_path,mode,out,reform_script)
        )
for i in range(1):
    n = 2293
    out = "auto_all"
    target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/auto_og_pass.txt"
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}".format(out=out)
    gwf.target_from_template(
            name = "sample_{out}".format(out=out),
            template = random_sample(path, target, n, out)
        )
    gwf.target_from_template(
        name = "cat_{out}".format(out=out),
        template = cat_og_phy_filter(path,out,out+".fa")
    )
    codeml_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/codeml"
    mode = "ofree"
    reform_script = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/reform_results.py"
    alignment = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}/{out}.fa".format(out=out)
    gwf.target_from_template(
            name = "paml_{mode}_{id_list}".format(mode=mode,id_list=out),
            template = paml_filter(alignment,path,codeml_path,mode,out,reform_script)
        )
for i in range(1):
    n = 347
    out = "auto_x"
    target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/x_og_pass.txt"
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/x/{out}".format(out=out)
    gwf.target_from_template(
            name = "sample_{out}".format(out=out),
            template = random_sample(path, target, n, out)
        )
    gwf.target_from_template(
        name = "cat_{out}".format(out=out),
        template = cat_og_phy_filter(path,out,out+".fa")
    )
    codeml_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/codeml"
    mode = "ofree"
    reform_script = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/reform_results.py"
    alignment = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}/{out}.fa".format(out=out)
    gwf.target_from_template(
            name = "paml_{mode}_{id_list}".format(mode=mode,id_list=out),
            template = paml_filter(alignment,path,codeml_path,mode,out,reform_script)
        )

for i in range(100):
    n = 100
    out = "x_{n}_{j}".format(n=n,j=i+1)
    target = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/x_og_pass.txt"
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/x/{out}".format(out=out)
    gwf.target_from_template(
            name = "sample_{out}".format(out=out),
            template = random_sample(path, target, n, out)
        )
    gwf.target_from_template(
        name = "cat_{out}".format(out=out),
        template = cat_og_phy_filter(path,out,out+".fa")
    )
    codeml_path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/codeml"
    mode = "ofree"
    reform_script = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/7_sp/reform_results.py"
    alignment = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/boostrap/auto/{out}/{out}.fa".format(out=out)
    gwf.target_from_template(
            name = "paml_{mode}_{id_list}".format(mode=mode,id_list=out),
            template = paml_filter(alignment,path,codeml_path,mode,out,reform_script)
        )
