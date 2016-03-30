#!/usr/bin/env python

import sys
import subprocess
import os
import pandas as pd

organism = sys.argv[1]
local_mirror = sys.argv[2]
ftp_directory = 'bacteria/' + organism + '/latest_assembly_versions/'
just_fastas = local_mirror + '_fastas/'

if os.path.isdir(organism):
    subprocess.run(['rsync',
                    '-iPrLtm',
                    '--exclude=**/unplaced_scaffolds/**',
                    '-f=+ *.fna.gz',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/' + ftp_directory,
                    '--log-file=log.txt',
                    local_mirror])

else:
    os.mkdir(organism)
    subprocess.run(['rsync',
                    '-iPrLtm',
                    '--exclude=**/unplaced_scaffolds/**',
                    '-f=+ *.fna.gz',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/' + ftp_directory,
                    '--log-file=log.txt',
                    local_mirror])

# copy just fastas from local NCBI mirror to just_fastas directory
if os.path.isdir(just_fastas):
    subprocess.run(['sudo', 'find', organism, '-type', 'f',
                    '-exec', 'cp',
                    '-t', just_fastas,
                    '-- {}', '+'])

else:
    os.mkdir(just_fastas)
    subprocess.run(['find', organism, '-type', 'f',
                '-exec', 'cp',
                '-t', just_fastas,
                '-- {}', '+'])

# Rename files in just_fastas directory
summary = sys.argv[3]
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
            newname = '{}_{}{}_{}.fna.gz'.format(id, org_name, strain, assembly_level)
            old = just_fastas+f
            new = just_fastas+newname
            os.rename(old, new)

# Decompress files in just_fastas directory
subprocess.run(['pigz', '-d', just_fastas+'*.fna.gz'])
