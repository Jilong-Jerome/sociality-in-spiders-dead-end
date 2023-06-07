from Bio import SeqIO
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

sequences = list(SeqIO.parse(input_file, "fasta"))

# Concatenate all the sequences into a single long sequence
concat_seq = sum(sequences, start=sequences[0].__class__(""))

# Update the ID and description of the concatenated sequence
concat_seq.id = output_file.split(".")[0]
#concat_seq.description = "Concatenated sequences from " + input_file


# Write the concatenated sequence to a new FASTA file
with open(output_file, "w") as out_handle:
    SeqIO.write(concat_seq, out_handle, "fasta")

