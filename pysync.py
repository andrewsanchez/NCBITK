#!/usr/bin/env python

import sys
import subprocess
import os
from ftplib import FTP

ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')

local_mirror = 'genbank_bacteria/'
all_fastas = local_mirror + '_fastas'
if not os.path.isdir(local_mirror):
    os.mkdir(local_mirror)

if not os.path.isdir(all_fastas):
    os.mkdir(all_fastas)

dirs = ftp.nlst()
for organism in dirs[0:100]:
    print(str(dirs.index(organism))+ ' out of ' + str(len(dirs)))

    just_fastas = local_mirror + organism + '_fastas'
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
