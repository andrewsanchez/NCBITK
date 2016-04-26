#!/usr/bin/env python3

import os, sys, stat, subprocess, rename
import pandas as pd
from ftplib import FTP

local_mirror = sys.argv[1]
all_fastas = local_mirror.strip('/') + '_fastas/'

# get current version of assembly_summary.txt
#if os.path.isfile('assembly_summary.txt'):
#  os.remove('assembly_summary.txt')
#  subprocess.call(['wget',
#                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])
#else:
#  subprocess.call(['wget',
#                  'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])

# connect with ftp.ncbi and get list of directories
ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')
dirs = ftp.nlst()

# make directories if necessary
if not os.path.isdir(local_mirror):
    os.mkdir(local_mirror)

if not os.path.isdir(all_fastas):
    os.mkdir(all_fastas)

for organism in dirs[0:3]: # sync with any number of folders with dirs[m:n]
    single_organism = all_fastas + organism + '/'
    subprocess.call(['rsync',
                    '--ignore-existing',
                    '--chmod=777',
                    '-irLtm',
                    '-f=+ GCA*fna.gz',
                    '--exclude=/unplaced_scaffolds/**',
                    '-f=+ */',
                    '--exclude=*',
                    'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/' + organism + '/latest_assembly_versions/',
                    '--log-file='+local_mirror+'log.txt',
                    local_mirror + '/' + organism])

    # copy files to different directory
    if os.path.isdir(single_organism):
        organism = local_mirror + '/' + organism
        subprocess.call(['find', organism, '-type', 'f',
                        '-exec', 'cp',
                        '-t', single_organism,
                        '-- {}', '+'])

        # decompress files
        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz', '-k', '-fd', '-- {}', '+'])
    else:
        os.mkdir(single_organism)
        organism = local_mirror + '/' + organism
        subprocess.call(['find', organism, '-type', 'f',
                        '-exec', 'cp',
                        '-t', single_organism,
                        '-- {}', '+'])

        subprocess.call(['find', single_organism, '-name', '*.gz',
                        '-exec', 'pigz', '-k', '-fd', '-- {}', '+'])

    rename.rename(single_organism)

