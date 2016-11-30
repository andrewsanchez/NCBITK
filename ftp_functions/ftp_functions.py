import os
import argparse
from ftplib import FTP, error_temp
from urllib.request import urlretrieve
from urllib.error import URLError

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


def grab_zipped_genome(genbank_mirror, species, genome_id, genome_path, ext=".fna.gz"):

    """
    Download compressed genome from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
    """

    zipped = "{}_genomic{}".format(genome_path, ext)
    zipped_dst = "{}_genomic{}".format(genome_id, ext)
    zipped_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}".format(genome_path, zipped)
    zipped_dst = os.path.join(genbank_mirror, species, zipped_dst)
    urlretrieve(zipped_url, zipped_dst)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("species")
    parser.add_argument("genome_id")
    parser.add_argument("genome_path")
    parser.add_argument("-g", "--grab", action="store_true")
    args = parser.parse_args()

    if args.grab:
        grab_zipped_genome(args.genbank_mirror, args.species, args.genome_id, args.genome_path)

main()

