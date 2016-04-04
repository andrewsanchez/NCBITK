#!/usr/bin/env python

import sys
import subprocess
import os
import pandas as pd
from ftplib import FTP

ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')

summary = 'assembly_summary.txt' # location of assembly_summary.txt
df = pd.read_csv(summary, delimiter='\t', index_col=0)

df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.infraspecific_name.replace({'[ =\-\;]': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)

local_mirror = 'genbank_bacteria/'
all_fastas = local_mirror.strip('/') + '_fastas/'

if not os.path.isdir(local_mirror):
    os.mkdir(local_mirror)

if not os.path.isdir(all_fastas):
    os.mkdir(all_fastas)

dirs = ftp.nlst()
for organism in dirs[0:100]:
    print(str(dirs.index(organism))+ ' out of ' + str(len(dirs)))

    just_fastas = all_fastas + organism + '/'
    subprocess.run(['rsync',
                    '-iPrLt',
                    '-f=+ GCA*fna.gz',
                    '--exclude=/unplaced_scaffolds/**',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/' + organism + '/latest_assembly_versions/',
                    '--log-file=log.txt',
                    local_mirror + organism])
    if os.path.isdir(just_fastas):
        organism = local_mirror + organism
        subprocess.run(['sudo', 'find', organism, '-type', 'f',
                        '-exec', 'cp',
                        '-t', just_fastas,
                        '-- {}', '+'])

        subprocess.run(['sudo', 'find', just_fastas, '-name', '*.gz',
                        '-exec', 'pigz',
                        '-d',
                        '-- {}', '+'])
    else:
        os.mkdir(just_fastas)
        organism = local_mirror + organism
        subprocess.run(['find', organism, '-type', 'f',
                    '-exec', 'cp',
                    '-t', just_fastas,
                    '-- {}', '+'])

        subprocess.run(['sudo', 'find', just_fastas, '-name', '*.gz',
                        '-exec', 'pigz',
                        '-d',
                        '-- {}', '+'])

    for f in os.listdir(just_fastas):
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
