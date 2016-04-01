#!/usr/bin/env python

import pandas as pd
import os
import sys

just_fastas = sys.argv[1] # *_just_fastas
summary = sys.argv[2] # assembly_summary.txt
df = pd.read_csv(summary, delimiter='\t', index_col=0)

fna_files = []
correctly_named = []
incorrectly_named = []

for f in os.listdir(just_fastas):
    fna_files.append(f)
    abs_path = just_fastas + f
    fasta = open(abs_path, 'r')
    line = fasta.readline()
    organism = f.split('_')[2:3]
    organism = '_'.join(organism)
    if organism in line:
        correctly_named.append(f)
    else:
        incorrectly_named.append(f)
    fasta.close()

print('Total number of files:  {}'.format(len(fna_files)))
print('Correctly named files:  {}'.format(len(correctly_named)))
print('Incorrectly named files:  {}'.format(len(incorrectly_named)))
print('\n')

for f in incorrectly_named:
    abs_path = just_fastas + f
    fasta = open(abs_path, 'r')
    line = fasta.readline()
    print(line)

    #organism = f.split('_')[2:3]
    #organism = '_'.join(organism)
    #print(organism + '<- Filename')
    print(f)

    id = (f.split('_')[0:2])
    id = ('_'.join(id))
    for row in df.index:
        if id == row:
            print(df.get_value(id, 'organism_name'))

print('\n')
