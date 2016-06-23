#!/usr/bin/env python3

import re, sys, os, subprocess
from ftplib import FTP

def get_dirs():
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd('genomes/genbank/bacteria')
    dirs = ftp.nlst()
    return dirs

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

def get_fastas(local_mirror, dirs='y'):
    if dirs == 'y':
        ftp_site = 'ftp.ncbi.nlm.nih.gov'
        ftp = FTP(ftp_site)
        ftp.login()
        ftp.cwd('genomes/genbank/bacteria')
        dirs = ftp.nlst()
    else:
        pass
    for organism in dirs:
        latest = os.path.join(organism, 'latest_assembly_versions')
        for parent in ftp.nlst(latest):
            accession = parent.split("/")[-1]
            fasta = accession+"_genomic.fna.gz"
            subprocess.call(['rsync',
                            '--info=progress1,flist2,copy1',
                            #'-v',
                            #'-vv',
                            #'--dry-run',
                            #'--itemize-changes',
                            '--chmod=ugo=rwX',
                            '--copy-links',
                            '--recursive',
                            '--times',
                            '--prune-empty-dirs',
                            '-f=+ '+fasta,
                            '-f=+ */',
                           #'-f=+ '+organism+'/',
                           #'-f=+ '+organism+'/'+'latest_assembly_versions/',
                           #'-f=+ '+organism+'/'+'latest_assembly_versions/'+accession+'/',
                            '--exclude=*',
                            '--log-file='+local_mirror+'/log.txt',
                            #'--log-file-format=%n',
                            'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/'+parent,
                            local_mirror + '/' + organism])

def main(argv):
    get_fastas('scratch/test_dir', argv)

if __name__ == "__main__":
    main(sys.argv[1])
