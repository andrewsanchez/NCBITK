#!/usr/bin/env python3

import re, sys, os, subprocess
from ftplib import FTP

#exclude = re.compile("(.*(from|rna|cds|protein).*)|(.*txt)|()")
def count_dirs():
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()
    dir_count = 0
    for organism in dirs:
        try:
            dir_count += len(ftp.nlst(organism+"/latest_assembly_versions"))
            with open("dir_count.txt", "a") as file:
                file.write(dir_count)
                file.write("\n")
        except:
            continue

def count_files():
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()
    for organism in dirs:
        ftp = FTP(ftp_site)
        ftp.login()
        ftp.cwd('genomes/genbank/bacteria')
        try:
            for parent in ftp.nlst(organism+"/latest_assembly_versions"):
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

def get_fastas(local_mirror):
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()

    rsync_log = os.path.join(local_mirror, "rsync_log.txt")
    for organism in dirs:
        ftp = FTP(ftp_site)
        ftp.login()
        ftp.cwd('genomes/genbank/bacteria')
        latest = os.path.join(organism, 'latest_assembly_versions')
        for parent in ftp.nlst(latest):
            accession = parent.split("/")[-1]
            fasta = accession+"_genomic.fna.gz"
            organism_dir = os.path.join(local_mirror, organism)
            subprocess.call(['rsync',
                            '--dry-run',
                            '--chmod=ugo=rwX',
                            '--copy-links',
                            '--recursive',
                            '--times',
                            '--prune-empty-dirs',
                            '-f=+ '+fasta,
                            '-f=+ */',
                            '--exclude=*',
                            'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/'+parent,
                            '--log-file='+rsync_log,
                            #'--itemize-changes',
                           #'-f=+ '+organism+'/',
                           #'-f=+ '+organism+'/'+'latest_assembly_versions/',
                           #'-f=+ '+organism+'/'+'latest_assembly_versions/'+accession+'/',
                            organism_dir])

def main(argv):
    if argv == "files":
        count_files()
    if argv == "dirs":
        count_dirs()
    if argv == "rsync":
        get_fastas('scratch/test_dir')

if __name__ == "__main__":
    main(sys.argv[1])
