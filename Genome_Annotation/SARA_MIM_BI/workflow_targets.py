from workflow_templates import *
sara_dict = {}
david_ind = open("SARA_david_name.txt")
for line in david_ind:
    ind = line.strip("\n")
    fq1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/Sarasinorum_RNA_David/Clean/{ind}/{ind}_1.fq.gz".format(ind=ind)
    fq2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/Sarasinorum_RNA_David/Clean/{ind}/{ind}_2.fq.gz".format(ind=ind)
    sara_dict[ind] = [fq1,fq2]
jesper_ind = open("SARA_jesper_name.txt")
for line in jesper_ind:
    ind = line.strip("\n")
    fq1 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/S_sarasinorum/{ind}/{ind}_R1_001.fastq.gz".format(ind=ind)
    fq2 = "/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/S_sarasinorum/{ind}/{ind}_R2_001.fastq.gz".format(ind=ind)
    sara_dict[ind] = [fq1,fq2]
mim_dict = {}
mim_info = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/MIM_33_RNA_meta.txt")
for line in mim_info:
    line_info = line.strip("\n")
    ind = line_info.split("\t")[0].split("/")[-2]
    fq1 = line_info.split("\t")[0]
    fq2 = line_info.split("\t")[1]
    mim_dict[ind] = [fq1,fq2]
bi_dict = {}
bi_info = open("/home/jilong/spider2/faststorage/social_spiders_2020/data/BACKUP/RNAseq_2021/BI_18_RNA_meta.txt")
for line in bi_info:
    line_info = line.strip("\n")
    ind = line_info.split("\t")[0].split("/")[-2]
    fq1 = line_info.split("\t")[0]
    fq2 = line_info.split("\t")[1]
    bi_dict[ind] = [fq1,fq2]
sp_fq_dict = {
    "SARA":sara_dict,
    "MIM":mim_dict,
    "BI":bi_dict,
    "MIM16":mim_dict
}
def build_STAR_index(path,ref,sp,log,gwf):
    gwf.target_from_template(
        name = "STAR_index_{sp}".format(sp=sp),
        template = STAR_index(path,ref,log)
        )
    return LOG_PATH+"/{log}".format(log=log)

def STAR_align_dict(sp_dict,sp,align_path,index_path,gwf):
    log_list = []
    for ind in sp_dict:
        fq1 = sp_dict[ind][0]
        fq2 = sp_dict[ind][1]
        indname = ind.replace("-","_")
        if sp == "MIM16":
            task_name = "STAR_mapping_MIM16_{ind}".format(ind=indname)
        else:
            task_name = "STAR_mapping_{ind}".format(ind=indname)
        gwf.target_from_template(
            name = task_name,
            template = STAR_mapping(align_path,index_path,sp,fq1,fq2,ind)
            )
        log_list.append(LOG_PATH+"/STAR_align_"+ind+".DONE")
    return log_list

def merge_STAR_bam(path,sp,filepath,filenames,logs,gwf):
    file_string = ""
    for filename in filenames:
        file_dir = filepath+"/"+filename+"/{filename}Aligned.sortedByCoord.out.bam".format(filename=filename)
        file_string = file_string+file_dir+" "
    outname = "STAR_merged_{sp}".format(sp=sp)
    gwf.target_from_template(
            name = "STAR_merge_bam_{sp}".format(sp=sp),
            template = merge_bams(path,outname,file_string,logs)
            )
    log = LOG_PATH+"/STAR_merge_"+outname+".DONE"
    return log

def sort_STAR_bam(path,bam_in,bam_out,sp,gwf):
    gwf.target_from_template(
            name = "STAR_sort_bam_{sp}".format(sp=sp),
            template = sort_bam(path,bam_in,bam_out)
            )
    log = LOG_PATH+"/STAR_merge_"+bam_out+"_sort.DONE"
    return log

def model_sp(genome,sp,path,out,gwf):
    gwf.target_from_template(
            name = "repeat_model_{sp}".format(sp=sp),
            template = repeat_model(genome,out,path)
            )
    log = LOG_PATH+"/repeat_model_{out}.DONE".format(out=out)
    return log
def repeat_sp(genome,lib1,lib2,lib,path,sp,gwf):
    gwf.target_from_template(
            name = "repeat_mask_{sp}".format(sp=sp),
            template = repeat_masker_lib(path,genome,lib1,lib2,lib,sp)
            )
    log = [LOG_PATH+"/repeat_mask_{sp}.DONE".format(sp=sp)]
    return log

def grep_chrom_fa(sp,numbers,gwf):
    log_list = []
    for chrom_id in range(numbers):
        chrom_id = chrom_id + 1
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/split".format(sp=sp)
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/mask/{sp}_hifi_hic_scaffolded_trim.fa.masked".format(sp=sp)
        out = "{sp}_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        task_name = "{sp}_get_fa_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        gwf.target_from_template(
            name=task_name,
            template = seqkit_grep(sp,path,genome,chrom_id,out)
            )
        log_list.append(LOG_PATH+"/seqkit_grep_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id))
    return log_list

def braker_sp_chrom(sp,numbers,gwf):
    log_list = []
    for chrom_id in range(numbers):
        chrom_id = chrom_id + 1
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/split/{sp}_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged/STAR_merged_{sp}_sort.bam".format(sp=sp) 
        protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
        task_name = "braker_{sp}_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        aug_sp = "{sp}_braker_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        gwf.target_from_template(
            name=task_name,
            template = braker_combine_mode(path,sp,aug_sp,chrom_id,genome,protein,bam)
            )
        log_list.append(LOG_PATH+"/BRAKER_ETP_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id))
    return log_list 

def braker_sp_chrom_list(sp,lists,gwf):
    log_list = []
    for chrom_id in lists:
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/repeat_masker/split/{sp}_hifi_hic_scaffolded_trim.fa.masked.{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        bam = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_merged/STAR_merged_{sp}_sort.bam".format(sp=sp) 
        protein = "/home/jilong/spider2/faststorage/social_spiders_2020/data/public_data/orthodb/combined_dumicola_arthopoda_orthoddb.fa"
        task_name = "braker_{sp}_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        aug_sp = "{sp}_braker_{chrom_id}".format(sp=sp,chrom_id=chrom_id)
        gwf.target_from_template(
            name=task_name,
            template = braker_combine_mode(path,sp,aug_sp,chrom_id,genome,protein,bam)
            )
        log_list.append(LOG_PATH+"/BRAKER_ETP_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id))
    return log_list 
def agat_gff2fa(sp,chrom_numbers,gwf):
    log_list = []
    for chrom_id in chrom_numbers:
        path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}/fasta".format(sp=sp,chrom_id = chrom_id)
        gff = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}/braker/braker.gff3".format(sp=sp,chrom_id=chrom_id)
        genome = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
        fasta = "braker_{sp}_HiC_scaffold_{chrom_id}.fasta".format(sp=sp,chrom_id=chrom_id)
        gwf.target_from_template(
            name = "gff2fa_{sp}_{chrom_id}".format(sp=sp,chrom_id=chrom_id),
            template = agat_gff2fasta(path,gff,genome,fasta,sp,chrom_id)
        )
        log_list.append(LOG_PATH+"/gff2fasta_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id))
        fasta = "braker_protein_{sp}_HiC_scaffold_{chrom_id}.fasta".format(sp=sp,chrom_id=chrom_id)
        gwf.target_from_template(
            name = "gff2fa_prot_{sp}_{chrom_id}".format(sp=sp,chrom_id=chrom_id),
            template = agat_gff2fasta_prot(path,gff,genome,fasta,sp,chrom_id)
        )
        log_list.append(LOG_PATH+"/gff2prot_fasta_{sp}_{chrom}_.DONE".format(sp=sp,chrom=chrom_id))
    return log_list

def combine_results_braker_fa(sp,chrom_numbers,gwf):
    string_gff = ""
    string_n = ""
    string_p = ""
    for chrom_id in chrom_numbers:
        gff_file = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}/braker/braker.gff3".format(sp=sp,chrom_id=chrom_id)
        n_file = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}/fasta/braker_{sp}_HiC_scaffold_{chrom_id}.fasta".format(sp=sp,chrom_id=chrom_id)
        p_file = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/HiC_scaffold_{chrom_id}/fasta/braker_protein_{sp}_HiC_scaffold_{chrom_id}.fasta".format(sp=sp,chrom_id=chrom_id)
        string_gff = string_gff + " " + gff_file
        string_n = string_n + " " + n_file
        string_p = string_p + " " + p_file
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker".format(sp=sp)
    gwf.target_from_template(
        name = "cat_BRAKER_{sp}".format(sp=sp),
        template = cat_braker_sp(sp,path,string_gff,string_n,string_p)
    )

def busco_sp(sp,gwf):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/busco".format(sp=sp)
    fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/braker/{sp}_braker.protein.rename.fasta".format(sp=sp)
    outname = "braker_{sp}_busco".format(sp=sp)
    gwf.target_from_template(
        name = "BUSCO_BRAKER_{sp}".format(sp=sp),
        template = run_busco_protein(sp,fasta,path,outname)
    )

def busco_per_sp(sp,gwf):
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/BUSCO/{sp}/busco".format(sp=sp)
    fasta = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/PUB1_GENOME/BUSCO/{sp}/fasta/{sp}_all.fa".format(sp=sp)
    outname = "{sp}_busco".format(sp=sp)
    gwf.target_from_template(
        name = "BUSCO_annotation_{sp}".format(sp=sp),
        template = run_busco_protein(sp,fasta,path,outname)
    )
def build_STAR_index_sp(sp,gwf):
    #sp = "SARA"
    path = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/{sp}/star_index".format(sp=sp)
    ref = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/3D_dna/final_3d/{sp}/{sp}_hifi_hic_scaffolded_trim.fa".format(sp=sp)
    log = "star_index_{sp}.DONE".format(sp=sp)
    build_STAR_index(path,ref,sp,log,gwf)
