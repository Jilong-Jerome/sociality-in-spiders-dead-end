import sys
def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequences = {}
        current_seq = ''
        current_header = ''
        
        for line in file:
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = current_seq
                
                current_header = line.strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        
        if current_header:
            sequences[current_header] = current_seq
    
    return sequences


def replace_invalid_nucleotides(sequences):
    valid_nucleotides = 'ATCG'
    
    for header, sequence in sequences.items():
        replaced_seq = ''.join([nuc if nuc in valid_nucleotides else '-' for nuc in sequence])
        sequences[header] = replaced_seq
    
    return sequences


def write_fasta(file_path, sequences):
    with open(file_path, 'w') as file:
        for header, sequence in sequences.items():
            file.write(header + '\n')
            file.write(sequence + '\n')


if __name__ == '__main__':
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    # Read the input FASTA file
    sequences = read_fasta(input_fasta)

    # Replace invalid nucleotides
    replaced_sequences = replace_invalid_nucleotides(sequences)

    # Write the output FASTA file
    write_fasta(output_fasta, replaced_sequences)

