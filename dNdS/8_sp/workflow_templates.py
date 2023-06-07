from gwf import *
LOG_PATH="/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/8_dnds/logs"
def logs_sum(logs,sum_log):
    inputs = logs
    outputs = [sum_log]
    options = {
              'cores': 1,
              'memory': '1g',
              'walltime':'1:00:00',
              'account':"spider2"
              }
    spec = """
    echo "finished" > {sum_log}
    echo date
""".format(sum_log=sum_log)
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def retrieve_ortholog_fasta(path,pac_cds,pac_name,group):
    inputs = []
    outputs = [
              LOG_PATH+"/ortho_fasta_{group}.DONE".format(group=group)
              ]
    options = {
              'cores':1,
              'memory':'100m',
              'walltime':'00:30:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate seqkit
    mkdir -p {path}
    cd {path}
    seqkit grep -p {pac_name} {pac_cds} > PAC_{group}.fa
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/rename_fasta.py PAC_{group}.fa PAC_{group}_rename.fa PAC
    cat /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/7_dnds/prank_ortho/{group}.rename.fasta PAC_{group}_rename.fa > {group}_unalign.fasta
    rm PAC_{group}.fa
    rm PAC_{group}_rename.fa
    echo done > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,pac_name=pac_name,pac_cds=pac_cds,group=group,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def align_ortho(path,group):
    inputs = [LOG_PATH+"/ortho_fasta_{group}.DONE".format(group=group)]
    outputs = [
              LOG_PATH+"/ortho_align_{group}.DONE".format(group=group)
              ]
    options = {
              'cores':1,
              'memory':'4g',
              'walltime':'02:00:00',
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
    conda activate macse
    macse -prog alignSequences -seq {group}_unalign.fasta -out_AA {group}_AA.fasta -out_NT {group}_NT.fasta -fs 10 -stop 10
    conda activate biopython
    python /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/PUB1_GENOME/dnds/8_sp/fix_fasta.py {group}_NT.fasta {group}_NT_fix.fasta 
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,group=group,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def paml(path,mode,group):
    inputs = [LOG_PATH+"/ortho_prank_{group}.DONE".format(group=group)]
    outputs = [
              LOG_PATH+"/paml/paml_{mode}_{group}.DONE".format(mode=mode,group=group)
              ]
    options = {
              'cores':1,
              'memory':'100m',
              'walltime':'00:30:00',
              'account':'spider2'
              }
    spec = """
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate paml
    echo "start prank aligning single ortholog groups fasta"
    echo jobinfo $SLURM_JOBID
    date
    mkdir -p {path}
    cd {path}
    cp ../../codeml_{mode}.ctl codeml.ctl
    ln -s -f /home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/dnds/7_dnds/prank_ortho/{group}_NT.phy seq.txt
    codeml     
    echo {group} > {group}.name
    paste {group}.name rst1 > {group}.res
    rm {group}.name
    cat results.txt |tail -n 26|head -n 12|grep "\S"| tr -s ' ' > results.tab
    python {reform_script} results.tab {group} {group}.tab
    rm results.tab
    echo finished > {log}
    date
    jobinfo $SLURM_JOBID
    """.format(path=path,mode=mode,group=group,log=outputs[0])
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



