import sys

sp = sys.argv[1]
gff_in = sys.argv[2]
gff_out = sys.argv[3]

def read_gff_gene(gene_in):
    info = gene_in.strip("\n").split("\t")
    chrom = sp+"."+"HiC_"+info[0].split("_")[-1]
    method = info[1]
    note = info[2]
    start = info[3]
    end = info[4]
    strand = info[6]
    gene_id = info[8].split("=")[-1].strip(";")
    new_id = chrom+"_"+gene_id
    out_string = "{chrom}\t{method}\t{note}\t{start}\t{end}\t\t{strand}\t\tlocus={new_id}\n".format(chrom=chrom,method=method,note=note,start=start,end=end,strand=strand,new_id=new_id)
    return out_string

infile = open(gff_in)
out = open(gff_out,"w")
for line in infile:
    outstring = read_gff_gene(line)
    out.write(outstring)
