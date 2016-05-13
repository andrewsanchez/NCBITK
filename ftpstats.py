#!/usr/bin/env python3

import os, argparse
from ftplib import FTP

def compare_dirs(organisms, local):
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/'+organisms)
    dirs = ftp.nlst()
    countdirs = str(len(dirs))
    localdirs = str(len(os.listdir(local)))
    missingdirs = int(countdirs) - int(localdirs)
    print('Dirs at ftp.ncbi:  '+countdirs)
    print('Dirs in {}:  {}'.format(local, localdirs))
    print('Missing dirs = ' + str(missingdirs))

    for i in dirs:
        if i not in os.listdir(local):
            print(i)

def Main():
    parser = argparse.ArgumentParser(
    description = "Compare the number of directories at ftp.ncbi.nlm.nih.gov/genomes/genbank/[organism]"\
            "with the number of directories in your local directory")
    parser.add_argument("organism", help = "e.g., bacteria or fungi")
    parser.add_argument("local_dir", help = "local directory you want to compare against NCBI's FTP site")
    args = parser.parse_args()

    compare_dirs(args.organism, args.local_dir)

Main()
