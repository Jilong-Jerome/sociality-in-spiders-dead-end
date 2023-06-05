sp_gff_path_dict = {
    "MIM":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM/braker",
    "SARA":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/SARA/braker",
    "BI":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/BI/braker",
    "DUM":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/etp_mode",
    "TENT":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/tentoriicola/etp_mode",
    "LIN":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/lineatus/etp_mode"
}

sp_pep_path_dict = {
    "MIM":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/MIM/braker",
    "SARA":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/SARA/braker",
    "BI":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/full_annotation/BI/braker",
    "DUM":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/dumicola/etp_mode",
    "TENT":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/tentoriicola/etp_mode",
    "LIN":"/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/annotate/braker/lineatus/etp_mode"
}

of_fasta_dict = {
    "BI":"0",
    "DUM":"1",
    "LIN":"2",
    "MIM":"3",
    "SARA":"4",
    "TENT":"5"
}
chrom_group_file = open("ste_chr_map.tsv")
chrom_group_dict = {}
for line in chrom_group_file:
    line_info = line.strip("\n").split("\t")
    ste_name = line_info[0]
    ste_dict = {}
    for i in range(len(line_info)-1):
        scaffolds = line_info[i+1].split(";")
        chrom_list = []
        for scaffold in scaffolds:
            sp = scaffold.split(".")[0]
            chrom = scaffold.split(".")[1]
            chrom_list.append(chrom)
        ste_dict[sp]=chrom_list
    chrom_group_dict[ste_name]=ste_dict
        
