#!/usr/bin/env python3

import pandas as pd
import os
import sys
import subprocess

just_fastas = sys.argv[1] # folder whose contents you want to rename
summary = sys.argv[2] # location of assembly_summary.txt
df = pd.read_csv(summary, delimiter='\t', index_col=0)

df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.infraspecific_name.replace({'[ =\-\;]': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)

for f in os.listdir(just_fastas):
    if f.startswith('GCA'):
        id = (f.split('_')[0:2])
        id = ('_'.join(id))

    for row in df.index:
        if id == row:
            org_name = df.get_value(id, 'organism_name')
            strain = df.get_value(id, 'infraspecific_name').strip('strain=')
            assembly_level  = df.get_value(id, 'assembly_level')
            newname = '{}_{}{}_{}.fna'.format(id, org_name, strain, assembly_level)
            old = just_fastas+f
            new = just_fastas+newname
            os.rename(old, new)

# subprocess.run(['gunzip', '-f', '/home/truthling/MGGen/Acinetobacter_nosocomialis_fastas/*.fna.gz'])
