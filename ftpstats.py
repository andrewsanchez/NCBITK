#!/usr/bin/env python3

from ftplib import FTP
import os

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

if __name__ == "__main__":
    import sys
    compare_dirs(sys.argv[1], sys.argv[2])
