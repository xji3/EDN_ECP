# This file splits cdna (e.g. dna sequence from amino acid alignment)
# sequence alignment into 3 separate files by codon site position
# Xiang Ji
# xji3@ncsu.edu

import argparse, os

def main(args):
    alignment = args.alignment
    split_files = [alignment.replace('.fasta', '_CS_' + str(i)+'.fasta') for i in range(1, 4)]

    first_line = True
    with open(alignment, 'r') as f:
        for line in f:
            if line[0] == '>':
                if first_line:
                    for split_file in split_files:
                        with open(split_file, 'w+') as g:
                            g.write(line)
                    first_line = False
                else:
                    for split_file in split_files:
                        with open(split_file, 'a') as g:
                            g.write(line)
            else:
                seq = line.replace('\n', '')
                assert(len(seq) %3 == 0)
                for codon_site in range(1, 4):
                    site_seq = ''.join([seq[3*i + codon_site - 1] for i in range(len(seq)/3)])
                    with open(split_files[codon_site - 1], 'a+') as g:
                        g.write(site_seq + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', dest = 'alignment', required = True, help = 'Input alignment file')
    
    main(parser.parse_args())
                
        
