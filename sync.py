#!/usr/bin/env python3

import os, sys, stat, subprocess
import pandas as pd
from ftplib import FTP

# get current version of assembly_summary.txt
#if os.path.isfile('assembly_summary.txt'):
#   os.remove('assembly_summary.txt')
#   subprocess.call(['wget',
#                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])
#else:
#   subprocess.call(['wget',
#                   'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])

ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')

local_mirror = sys.argv[1] # make sure to include a trailing slash
all_fastas = local_mirror.strip('/') + '_fastas/'

# remove duplicate strings during renaming
def rmDuplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

# clean up assembly_summary.txt
df = pd.read_csv('assembly_summary.txt', delimiter='\t', index_col=0)
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].isnull())].fillna('NA'))
df.update(df['infraspecific_name'][(df['infraspecific_name'].isnull()) & (df['isolate'].notnull())].fillna(df['isolate']))
#df.organism_name.replace({' = ': '_'}, regex=True, inplace=True)
#df.infraspecific_name.replace({' = ': '_'}, regex=True, inplace=True)
df.organism_name.replace({' ': '_'}, regex=True, inplace=True)
df.organism_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.infraspecific_name.replace({'[ =\-\;]': '_'}, regex=True, inplace=True)
df.infraspecific_name.replace({'[\W]': ''}, regex=True, inplace=True)
df.assembly_level.replace({' ': '_'}, regex=True, inplace=True)

# make directories if necessary
if not os.path.isdir(local_mirror):
    os.mkdir(local_mirror)

if not os.path.isdir(all_fastas):
    os.mkdir(all_fastas)

dirs = ftp.nlst()
for organism in dirs[0:2]: # sync with any number of folders with dirs[m:n]
    single_organism = all_fastas + organism + '/'
    subprocess.call(['rsync',
                    '--ignore-existing',
                    '-irLtm',
                    '-f=+ GCA*fna.gz',
                    '--exclude=/unplaced_scaffolds/**',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/' + organism + '/latest_assembly_versions/',
                    '--log-file='+local_mirror+'log.txt',
                    local_mirror + organism])

    # change permissions so that files can be renamed
    for root, dirs, files in os.walk(local_mirror+organism):
        for dir in [os.path.join(root, d) for d in dirs]:
            os.chmod(dir, stat.S_IRWXU)
        for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, stat.S_IRWXU)

    # copy files to different directory
    if os.path.isdir(single_organism):
        organism = local_mirror + organism
        subprocess.call(['find', organism, '-type', 'f',
                        '-exec', 'cp',
                        '-nt', single_organism,
                        '-- {}', '+'])

        # decompress files
        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz', '-dfk', '-- {}', '+'])
    else:
        os.mkdir(single_organism)
        organism = local_mirror + organism
        subprocess.call(['find', organism, '-type', 'f',
                    '-exec', 'cp',
                    '-nt', single_organism,
                    '-- {}', '+'])

        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz', '-dfk', '-- {}', '+'])

    for f in os.listdir(single_organism):
        id = (f.split('_')[0:2])
        id = ('_'.join(id))
        print(f)

        if id in df.index:
            org_name = df.get_value(id, 'organism_name')
            strain = df.get_value(id, 'infraspecific_name')
            assembly_level  = df.get_value(id, 'assembly_level')
            newname = '{}_{}{}_{}.fna'.format(id, org_name, strain, assembly_level)
            newname = f.split('_')
            newname = rmDuplicates(newname)
            print('after rmDuplicates:  ' + newname)
            newname = '_'.join(newname)
            newname = newname.replace('subsp_', '')
            newname = newname.replace('str_', '')
            newname = newname.replace('strain_', '')
            print(newname)
            old = single_organism+f
            print(old)
            print(new)
            new = single_organism+newname
            os.rename(old, new)
