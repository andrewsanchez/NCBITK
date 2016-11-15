from ftplib import FTP, error_temp

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    print("Logging into ftp.ncbi.nlm.nih.gov/{}/".format(directory))
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

