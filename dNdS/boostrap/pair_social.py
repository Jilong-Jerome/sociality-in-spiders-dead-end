import sys
filename = sys.argv[1]
outname = sys.argv[2]

# Read the data and store it in a dictionary
data = {}
with open(filename, "r") as f:
    next(f)  # Skip the header
    for line in f:
        nodename, dN, dS, id = line.strip().split('\t')
        data[nodename] = {"dN": float(dN), "dS": float(dS),"id": id}
# Define the contrast pairs
contrast_pairs = [("AFR", "MIM"), ("PAC", "SARA"), ("TENT", "DUM")]

# Compute dN and dS ratios for each pair and store the results in a list
results = []
for nodename1, nodename2 in contrast_pairs:
    if nodename1 in data and nodename2 in data:
        dN_ratio = data[nodename2]["dN"] / data[nodename1]["dN"]
        dS_ratio = data[nodename1]["dS"] / data[nodename2]["dS"]
        dN_shift = data[nodename2]["dN"] - data[nodename1]["dN"]
        dS_shift = data[nodename2]["dS"] - data[nodename1]["dS"]
        dNdS_ratio = (data[nodename2]["dN"]/data[nodename2]["dS"])/(data[nodename1]["dN"]/data[nodename1]["dS"])
        dNdS_social = data[nodename2]["dN"]/data[nodename2]["dS"]
        dNdS_sub = data[nodename1]["dN"]/data[nodename1]["dS"]
        dN_social = data[nodename2]["dN"]
        dS_social = data[nodename2]["dS"]
        dN_sub = data[nodename1]["dN"]
        dS_sub = data[nodename1]["dS"]
        id = data[nodename1]["id"]
        results.append((id,nodename2, nodename1, dN_ratio, dS_ratio, dNdS_ratio,dN_shift,dS_shift,dNdS_social,dNdS_sub,dN_social,dS_social,dN_sub,dS_sub))

# Output the results as a tab-separated table
with open(outname, 'w') as outfile:
    # Write the header
    outfile.write("id\tsocial\tsubsocial\tdN_ratio\tdS_ratio\tdNdS_ratio\tdN_shift\tdS_shift\tdNdS_social\tdNdS_sub\tdN_social\tdS_social\tdN_sub\tdS_sub\n")

    # Write the table content
    for id, nodename2, nodename1, dN_ratio, dS_ratio, dNdS_ratio,dN_shift,dS_shift,dNdS_social,dNdS_sub,dN_social,dS_social,dN_sub,dS_sub in results:
        outfile.write(f"{id}\t{nodename2}\t{nodename1}\t{dN_ratio}\t{dS_ratio}\t{dNdS_ratio}\t{dN_shift}\t{dS_shift}\t{dNdS_social}\t{dNdS_sub}\t{dN_social}\t{dS_social}\t{dN_sub}\t{dS_sub}\n")

