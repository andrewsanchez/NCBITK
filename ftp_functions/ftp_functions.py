import os
from ftplib import FTP, error_temp

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    print(directory)
    print("Logging into ftp.ncbi.nlm.nih.gov")
    ftp.login()
    ftp.cwd(directory)

    return ftp

def ftp_complete_species_list():

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    ftp = ftp_login()

    print("Getting list of all bacteria directories.")
    print("Estimated wait time:  1 minute.")
    try:
        complete_species_list = ftp.nlst()
    except error_temp:
        sleep(30)
        complete_species_list = ftp.nlst()

    print("All bacteria directories succesfully read from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/")

    return complete_species_list

def check_dirs(genbank_mirror):

    """Make directories to store fastas if they don't already exist."""

    latest_assembly_versions_list = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")

    if not os.path.isdir(genbank_mirror):
        os.mkdir(genbank_mirror)

    if os.path.isfile(latest_assembly_versions_list):
        os.remove(latest_assembly_versions_list)

    info_dir = os.path.join(genbank_mirror, ".info")
    if not os.path.isdir(info_dir):
        os.mkdir(info_dir)
