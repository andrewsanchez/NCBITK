#!/usr/bin/env python

import pandas as pd
import os
import sys
import subprocess

just_fastas = sys.argv[1]
summary = 'assembly_summary.txt'
df = pd.read_csv(summary, delimiter='\t', index_col=0)

fna_files = []
correctly_named = []
incorrectly_named = []

for f in os.listdir('just_fastas'):
    fna_files.append(f)
    abs_path = 'just_fastas' + f
    fasta = open(abs_path, 'r')
    line = fasta.readline()
    if organism in line:
        correctly_named.append(f)
    else:
        incorrectly_named.append(f)
# print('Correctly named files:  {}'.format(len(correctly_named)))
# print('Incorrectly named files:  {}'.format(len(incorrectly_named)))
# print('Total number of files:  {}'.format(len(fna_files)))

for f in incorrectly_named:
    id = (f.split('_')[0:2])
    id = ('_'.join(id))
    organism = f.split('_')[2:3]
    organism = '_'.join(organism)
    abs_path = 'just_fastas' + f
    fasta = open(abs_path, 'r')
    line = fasta.readline()
    print(line)
    for row in df.index:
        if id == row:
            print(df.get_value(id, 'organism_name'))
    
#     print(f)
    abs_path = 'just_fastas' + f
    fna = open(abs_path, 'r')
#     print(fna.readline())
