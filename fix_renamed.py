#!/usr/bin/env python

import os, sys, re
import pandas as pd

summary = sys.argv[1] # location of assembly_summary.txt
df = pd.read_csv(summary, delimiter='\t', index_col=0)

funky_files = open('fix_renamed.txt', 'r')

for line in funky_files:
	id = (line.split('_')[0:2])
	id = ('_'.join(id))
	print(line).strip('\n')
	print(id)
	print(df.loc[id][6:8][:][0])
	print(df.loc[id][6:8][:][1])
	print('\n')
