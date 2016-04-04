#!/usr/bin/env python3

import os, sys, stat
import subprocess
import pandas as pd
from ftplib import FTP

ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')

local_mirror = sys.argv[1]
all_fastas = local_mirror.strip('/') + '_fastas/'
summary = sys.argv[2] # location of assembly_summary.txt

df = pd.read_csv(summary, delimiter='\t', index_col=0)

df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.infraspecific_name.replace({'[ =\-\;]': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)

if not os.path.isdir(local_mirror):
    os.mkdir(local_mirror)

if not os.path.isdir(all_fastas):
    os.mkdir(all_fastas)

dirs = ftp.nlst()
for organism in dirs[0:3]: # sync with any number of folders with dirs[n:n2]
    #print(str(dirs.index(organism))+ ' out of ' + str(len(dirs)))

    single_organism = all_fastas + organism + '/'
    subprocess.call(['rsync',
                    '-iPrLtm',
                    '-f=+ GCA*fna.gz',
                    '--exclude=/unplaced_scaffolds/**',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/' + organism + '/latest_assembly_versions/',
                    '--log-file='+local_mirror+'log.txt',
                    local_mirror + organism])

    for root, dirs, files in os.walk(local_mirror+organism):
        for dir in [os.path.join(root, d) for d in dirs]:
            os.chmod(dir, stat.S_IRWXU)
        for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, stat.S_IRWXU)

    if os.path.isdir(single_organism):
        organism = local_mirror + organism
        subprocess.call(['find', organism, '-type', 'f',
                        '-exec', 'cp',
                        '-t', single_organism,
                        '-- {}', '+'])

        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz',
                        '-d',
                        '-- {}', '+'])
    else:
        os.mkdir(single_organism)
        organism = local_mirror + organism
        subprocess.call(['find', organism, '-type', 'f',
                    '-exec', 'cp',
                    '-t', single_organism,
                    '-- {}', '+'])

        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz',
                        '-d',
                        '-- {}', '+'])

    for f in os.listdir(single_organism):
        id = (f.split('_')[0:2])
        id = ('_'.join(id))

        if id in df.index:
            org_name = df.get_value(id, 'organism_name')
            strain = df.get_value(id, 'infraspecific_name').strip('strain=')
            assembly_level  = df.get_value(id, 'assembly_level')
            newname = '{}_{}{}_{}.fna'.format(id, org_name, strain, assembly_level)
            old = single_organism+f
            new = single_organism+newname
            os.rename(old, new)
