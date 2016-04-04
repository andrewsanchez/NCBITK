#!/usr/bin/env python

import os
import sys

just_fastas = sys.argv[1] # *_just_fastas
summary = sys.argv[2] # assembly_summary.txt
df = pd.read_csv(summary, delimiter='\t', index_col=0)

total_fastas = []
correctly_named = []
misnamed = []

for root, dirs, files in os.walk('tmp/genbank_bacteria_fastas/'):
    for dir in dirs:
      genus = dir.split('_')[0]
    for name in files:
        total_fastas.append(name)
        abs_path = os.path.join(root, name)
        fasta = open(abs_path, 'r')
        line = fasta.readline()
        species = name.split('_')[3:4]
        species = '_'.join(species)
        if genus in line:
            correctly_named.append(name)
            fasta.close()
        else:
            file = os.path.join(root, name)
            misnamed.append(file)
            fasta.close()

print('Total number of files:  {}'.format(len(total_fastas)))
print('Correctly named files:  {}'.format(len(correctly_named)))
print('Misnamed files:  {}'.format(len(misnamed)))
print('Misnamed files:')
print('\n')

for f in misnamed:
    with open(f, 'r') as fasta:
        line = fasta.readline()
        print(f)
        print(line)
