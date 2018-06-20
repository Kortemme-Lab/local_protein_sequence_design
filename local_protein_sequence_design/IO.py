def sequence_to_fasta_file(fasta_path, title, sequence):
    '''Save a sequence into a fasta file.'''
    with open(fasta_path, 'w') as f:
        f.write('>' + title + '\n')
        
        start = 0
        end = min(80, len(sequence))
        f.write(sequence[start:end] + '\n')
        start += 80

        while start < len(sequence):
            end = min(start + 80, len(sequence))
            f.write(sequence[start:end] + '\n')
            start += 80

def fasta_file_to_sequences(fasta_file):
    '''Read a list of sequences from the fasta file.
    Return the titles and sequences'''
    titles = []
    sequences = []
  
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith('>'):
                titles.append(lines[i][1:])
                sequences.append('')
                i += 1

                while i < len(lines):
                    if lines[i].startswith('>'):
                        break
                    sequences[-1] = ''.join([sequences[-1], lines[i].strip()])
                    i += 1

    return titles, sequences

