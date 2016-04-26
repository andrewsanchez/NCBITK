#!/usr/bin/env python3

import os, sys, stat, subprocess, rename, argparse
import pandas as pd
from ftplib import FTP

# get current version of assembly_summary.txt
def get_assembly_summary():
    if os.path.isfile('assembly_summary.txt'):
       os.remove('assembly_summary.txt')
       subprocess.call(['wget',
                       'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])
    else:
       subprocess.call(['wget',
                       'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])

# connect with ftp.ncbi and get list of directories
def get_organism_list(source):
    if source:
        dirs = [] 
        try:
            with open('organism_list.txt', 'r') as f:
                for line in f:
                    line = line.strip()
                    dirs.append(line)
            return dirs

        except IOError:
            print('the --input_file flag requires a file "organism_list.txt"' \
                    ' to be in your current working directory.')

    else:
        ftp_site = 'ftp.ncbi.nlm.nih.gov'
        ftp = FTP(ftp_site)
        ftp.login()
        ftp.cwd('genomes/genbank/bacteria')
        dirs = ftp.nlst()

        return dirs

# make directories if necessary
def check_dirs(local_mirror):
    all_fastas = local_mirror.strip('/') + '_fastas/'
    if not os.path.isdir(local_mirror):
        os.mkdir(local_mirror)

    if not os.path.isdir(all_fastas):
        os.mkdir(all_fastas)

def get_fastas(local_mirror, organism_list):
    all_fastas = local_mirror.strip('/') + '_fastas/'
    try:
        dirs = get_organism_list(organism_list)
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
                            '--log-file='+local_mirror+'/log.txt',
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
    except TypeError:
        print('Create the file "organism_list.txt" in the current working directory')

def Main():
    parser = argparse.ArgumentParser(description = 'Sync with NCBI database.')
   #group = parser.add_mutually_exclusive_group()
   #group.add_argument('--input_file')
   #group.add_argument('--from_list')
    parser.add_argument('local_mirror', help = 'Your local directory to save fastas to, e.g. "bacteria"', type=str)
    parser.add_argument('-G', '--no_wget', help = "Don't fetch assembly_summary.txt", action='store_true')
    parser.add_argument('-i', '--input_file', help = 'Input file containing directories to sync with.  ' \
            'Requirements:  must be named organism_list.txt, exist in current working directory' \
            ' and have one directory per line.', action='store_true')
    parser.add_argument('-l', '--from_list', help = 'Comma separated list of directories to be downloaded'\
            'i.e., Escherichia_coli,Acinetobacter_nosocomialis', action='store_true')
    args = parser.parse_args()

    if args.no_wget:
        check_dirs(args.local_mirror)
        get_fastas(args.local_mirror, args.input_file)
    #if args.from_list:
    else:
        get_assembly_summary()
        check_dirs(args.local_mirror)
        get_fastas(args.local_mirror, args.input_file)

if __name__ == '__main__':
    Main()
