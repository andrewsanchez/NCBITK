import os
import argparse
from ftplib import FTP, error_temp
from urllib.request import urlretrieve
from urllib.error import URLError

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login(user='anonymous', passwd='aas229@nau.edu')
    ftp.cwd(directory)

    return ftp

def ftp_complete_species_list():

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    ftp = ftp_login()

    print("Getting list of all bacteria directories.")
    print("Estimated wait time:  1 minute.")
    try:
        complete_species_list = [i for i in ftp.nlst() if not i.endswith('txt')]
    except error_temp:
        sleep(30)
        complete_species_list = [i for i in ftp.nlst() if not i.endswith('txt')]

    print("{} bacteria directories succesfully read from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/".format(len(complete_species_list)))

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
    print(species,genome_id,genome_path)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("info", nargs='+')
    parser.add_argument("-g", "--grab", action="store_true")
    args = parser.parse_args()

    if args.grab:
        for grab_genome_args in args.info:
            grab_genome_args = grab_genome_args.split(',')
            species, genome_id, genome_path = grab_genome_args
            grab_zipped_genome(args.genbank_mirror, species, genome_id, genome_path)

if __name__ == "__main__":
    main()
