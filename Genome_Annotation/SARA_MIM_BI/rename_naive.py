import sys
in_file = sys.argv[1]
out_file = sys.argv[2]
chrom_col = int(sys.argv[3])
gene_col = int(sys.argv[4])
fasta = open(in_file)
out = open(out_file,"w")
for line in fasta:
    if line[0] == ">":
        chrom = line.split(" ")[chrom_col].split("=")[1]
        gene_id = line.split(" ")[gene_col].split("_")[-1]
        new_name = chrom+"_"+gene_id
        out.write(">{new_name}\n".format(new_name=new_name))
    else:
        out.write(line.replace("*",""))

