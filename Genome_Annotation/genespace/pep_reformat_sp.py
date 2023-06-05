import sys

sp = sys.argv[1]
pep_in = sys.argv[2]
pep_out = sys.argv[3]

def read_pep_line(line):
    if line[0] == ">":
        line_info = line.split(" ")
        gene = line_info[1].split("=")[1]
        chrom = line_info[2].split("=")[1].split("_")[2]
        new_string = ">{sp}.HiC_{chrom}_{gene}\n".format(sp=sp,chrom=chrom,gene=gene)
    else:
        new_string = line.replace("*","")
    return new_string

infile = open(pep_in)
out = open(pep_out,"w")
for line in infile:
    outstring = read_pep_line(line)
    out.write(outstring)
