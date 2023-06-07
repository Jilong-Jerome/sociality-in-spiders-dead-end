import sys
from Bio import SeqIO

def rename_fasta(input_file, output_file, new_id):
    with open(input_file, "r") as fasta_file:
        seq_record = SeqIO.read(fasta_file, "fasta")
        seq_record.id = new_id
        seq_record.description = ""

        with open(output_file, "w") as output:
            SeqIO.write(seq_record, output, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python rename_fasta.py <input_fasta> <output_fasta> <new_id>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    new_id = sys.argv[3]

    rename_fasta(input_file, output_file, new_id)

