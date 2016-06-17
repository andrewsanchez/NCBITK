#!/usr/bin/env python3

import re
from ftplib import FTP

ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)
ftp.login()
ftp.cwd('genomes/genbank/bacteria')
dirs = ftp.nlst()

#exclude = re.compile("(.*(from|rna|cds|protein).*)|(.*txt)|()")
def count_dirs():
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()
    dir_count = 0
    for dir in dirs:
        try:
            dir_count += len(ftp.nlst(dir+"/latest_assembly_versions"))
            print(dir_count)
        except:
            continue

def count_files():
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()
    for dir in dirs:
        try:
            for parent in ftp.nlst(dir+"/latest_assembly_versions"):
                id = parent.split("/")[-1]
                id = id+"_genomic.fna.gz"
                gz = re.compile(id)
                for child in ftp.nlst(parent):
                    child = child.split("/")[-1]
                    if gz.match(child):
                        with open("fastas_all.txt", "a") as file:
                            file.write(child)
                            file.write("\n")
        except:
            continue
