#!/usr/bin/env python

import os, rename_fastas, shutil, argparse, gzip
import pandas as pd
from re import sub
from urllib.request import urlretrieve
from urllib.error import URLError
from ftplib import FTP, error_temp
from time import strftime, sleep

def get_assembly_summary(genbank_mirror, assembly_summary_location, assembly_summary_url):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    urlretrieve(assembly_summary_url, assembly_summary_location)
    assembly_summary = pd.read_csv(assembly_summary_location, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

def check_dirs(genbank_mirror):

    """Make directories to store fastas if they don't already exist."""

    if not os.path.isdir(genbank_mirror):
        os.mkdir(genbank_mirror)

    info_dir = os.path.join(genbank_mirror, ".info")
    if not os.path.isdir(info_dir):
        os.mkdir(info_dir)

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

def write_latest_assembly_versions(genbank_mirror, species, ftp):

    latest_dir = os.path.join(species, "latest_assembly_versions")
    latest_assembly_versions = [accession.split("/")[-1] for accession in ftp.nlst(latest_dir)]
    dir_structure = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    with open(dir_structure, "a") as f:
        for accession in latest_assembly_versions:
            print("{} in {}".format(accession, species))
            f.write("{},{}\n".format(species, accession))

def get_latest_assembly_versions(genbank_mirror, complete_species_list, ymdt):

    """
    Create DataFrame to represent genbank's directory structure.
    """

    ftp = ftp_login()
    genbank_stats = os.path.join(genbank_mirror, "genbank_stats_{}.txt".format(ymdt))
    ymd = strftime("%y/%m/%d")
    no_latest_msg = "no latest_assembly_versions dir for"

    print("Getting latest assembly versions for {} species.".format(len(complete_species_list)))
    print("This will take several minutes.")

    def log_error():
            with open(genbank_stats, "a") as stats:
                stats.write("{} - {} {}\n".format(species, no_latest_msg, ymd))

    for species in complete_species_list:
        try:
            write_latest_assembly_versions(genbank_mirror, species, ftp)
        except error_temp:
            log_error()
         #  with open(genbank_stats, "a") as stats:
         #      stats.write("{} - {} {}\n".format(species, no_latest_msg, ymd))
        except BrokenPipeError:
            try:
                ftp = ftp_login()
                write_latest_assembly_versions(genbank_mirror, species, ftp)
            except error_temp:
                log_error()
              # with open(genbank_stats, "a") as stats:
              #     stats.write("{} - {} {}\n".format(species, no_latest_msg, ymd))

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None)
    latest_assembly_versions.columns = ["id"]

    return latest_assembly_versions

def check_species_dirs(genbank_mirror, species):

    species_dir = os.path.join(genbank_mirror, species)
    if not os.path.isdir(species_dir):
        os.mkdir(species_dir)

    return species_dir

def grab_zipped_genome(genbank_mirror, species, genome_id, ext=".fna.gz"):

    print(genome_id)
    zipped = "{}_genomic{}".format(genome_id, ext)
    zipped_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}".format(genome_id, zipped)
    zipped_dst = os.path.join(genbank_mirror, species, zipped)
    urlretrieve(zipped_url, zipped_dst)
    print("{} -> {}".format(zipped_url, zipped_dst))
    zipped_src = os.path.join(genbank_mirror, species, zipped)

    return zipped_src

def unzip_genome(genbank_mirror, zipped_src, species, genome_id):

    """
    Decompress genome and remove the compressed genome.
    """

    unzipped = "{}_genomic.fasta".format(genome_id)
    unzipped = os.path.join(genbank_mirror, species, unzipped)
    zipped = gzip.open(zipped_src)
    decoded = zipped.read()
    unzipped = open(unzipped, "wb")
    unzipped.write(decoded)
    unzipped.close()
    zipped.close()
    os.remove(zipped_src)

def mash():
    None

def log_errors(genbank_mirror):
    None

def grab_and_organize_genomes(genbank_mirror, latest_assembly_versions):

    species_directories = set(list(latest_assembly_versions.index))

    for species in species_directories:
        species_dir = check_species_dirs(genbank_mirror, species)
        local_genomes = ["_".join(genome_id.split("_")[:3]) for genome_id in os.listdir(species_dir)]
        latest_genomes = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species].values.tolist()]
        print(species, latest_genomes)

        for genome_id in local_genomes:
            if genome_id not in latest_genomes:
                os.remove(os.path.join(genbank_mirror, species, "{}_genomic.fasta".format(genome_id)))

        for genome_id in latest_genomes:
            if genome_id not in local_genomes:
                try:
                    zipped_src = grab_zipped_genome(genbank_mirror, species, genome_id)
                    unzip_genome(genbank_mirror, zipped_src, species, genome_id)
                except URLError:
                    zipped_src = grab_zipped_genome(genbank_mirror, species, genome_id, ext=".fasta.gz")
                    unzip_genome(genbank_mirror, zipped_src, species, genome_id)
                except URLError:
                    genbank_stats = os.path.join(genbank_mirror, "genbank_stats_{}.txt".format(ymdt))
                    with open(genbank_stats, "a") as stats:
                        stats.write("URLError for {}\n".format(genome_id))

def main():

    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror.strip("/")
#   ymdt = strftime("%y.%m.%d_%H%M")
#   assembly_summary_location = os.path.join(genbank_mirror, "assembly_summary.txt")
#   assembly_summary_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
#   check_dirs(genbank_mirror)
#   complete_species_list = ftp_complete_species_list()
#   assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_location, assembly_summary_url)
#   latest_assembly_versions = get_latest_assembly_versions(genbank_mirror, complete_species_list, ymdt)

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None)
    latest_assembly_versions.columns = ["id"]
    grab_and_organize_genomes(genbank_mirror, latest_assembly_versions)

main()
