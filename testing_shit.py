import os, rename_fastas, shutil
from subprocess import call
from ftplib import FTP, error_temp
from time import strftime

def ftp_complete_species_list(local_mirror, ftp):
   complete_species_list = ftp.nlst()
    ftp_stats = os.path.join(local_mirror, "ftp_stats.txt")
    time_stamp = strftime("%y/%m/%d_%H:%M")
    with open(ftp_stats, "a") as stats:
        stats.write("{} - {} dirs at genbank/bacteria/\n".format(time_stamp, len(complete_species_list)))
    return complete_species_list


def ftp_login(directory="genomes/genbank/bacteria"):
    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    print("Logging into ftp.ncbi.nlm.nih.gov/{}/".format(directory))
    ftp.login()
    ftp.cwd(directory)
    return ftp

