#!/usr/bin/env python

import os, shutil, argparse, gzip
import pandas as pd
import rename
from urllib.request import urlretrieve
from ftplib import FTP, error_temp
from urllib.error import URLError
from time import strftime, sleep
from re import sub

def get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")
    urlretrieve(assembly_summary_url, assembly_summary_dst)
    assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

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
    latest_dirs = [accession.split("/")[-1] for accession in ftp.nlst(latest_dir)]
    accession_ids = ["_".join(id.split("_")[:2]) for id in latest_dirs]
    latest_assembly_versions_list = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    dirs_and_ids = zip(latest_dirs, accession_ids)

    with open(latest_assembly_versions_list, "a") as f:
        for item in dirs_and_ids:
            f.write("{},{},{}\n".format(species, item[1], item[0]))

def get_latest_assembly_versions(genbank_mirror, complete_species_list, genbank_stats, ymdt):


    """
    Create DataFrame to represent genbank's directory structure.
    """

    ftp = ftp_login()

    print("Getting latest assembly versions for {} species.".format(len(complete_species_list)))
    print("Estimated wait time:  2 hours.")

    for species in complete_species_list:
        try:
            write_latest_assembly_versions(genbank_mirror, species, ftp)
        except error_temp:
            log_error("latest", species, genbank_stats, ymdt)
        except BrokenPipeError:
            try:
                ftp = ftp_login()
                write_latest_assembly_versions(genbank_mirror, species, ftp)
            except error_temp:
                log_error("latest", species, genbank_stats, ymdt)

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None)
    latest_assembly_versions.columns = ["id", "dir"]

    return latest_assembly_versions

def check_species_dirs(genbank_mirror, species):

    species_dir = os.path.join(genbank_mirror, species)
    if not os.path.isdir(species_dir):
        os.mkdir(species_dir)

    return species_dir

def grab_zipped_genome(genbank_mirror, species, genome_url, genome_id, ext=".fna.gz"):

    """
    Download compressed genome from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
    """

    zipped = "{}_genomic{}".format(genome_url, ext)
    zipped_dst = "{}_genomic.fasta".format(genome_id)
    zipped_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}".format(genome_url, zipped)
    zipped_dst = os.path.join(genbank_mirror, species, zipped)
    urlretrieve(zipped_url, zipped_dst)
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

def check_local_genomes(local_genomes, latest_genomes, species, genome_id, genbank_stats):

    for genome_id in local_genomes:
        if genome_id not in latest_genomes:
            os.remove(os.path.join(genbank_mirror, species, "{}_genomic.fasta".format(genome_id)))
            with open(genbank_stats, "a") as stats:
                stats.write("{} removed\n".format(genome_id))

def sync_latest_genomes(local_genomes, latest_genome_info, species, genome_id, genbank_stats):

    for info in latest_genomes:
        genome_id = info[0]
        genome_url = info[1]
        if genome_id not in local_genomes:
            try:
                zipped_src = grab_zipped_genome(genbank_mirror, species, genome_url)
                unzip_genome(genbank_mirror, zipped_src, species, genome_id)
            except URLError:
                zipped_src = grab_zipped_genome(genbank_mirror, species, genome_url, ext=".fasta.gz")
                unzip_genome(genbank_mirror, zipped_src, species, genome_id)
            except URLError:
                with open(genbank_stats, "a") as stats:
                    stats.write("URLError for {}\n".format(genome_id))
            with open(genbank_stats, "a") as stats:
                stats.write("{} downloaded\n".format(genome_id))

def grab_and_organize_genomes(genbank_mirror, genbank_stats, latest_assembly_versions):

    species_directories = set(list(latest_assembly_versions.index))

    for species in species_directories:
        species_dir = check_species_dirs(genbank_mirror, species)
        local_genomes = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
        latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        latest_genome_dirs = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
        latest_genome_info = zip(latest_genome_ids, latest_genome_dirs)
        print("species", species)
        print("local_genomes", local_genomes)
        print("latest_genomes", latest_genomes)

        check_local_genomes(local_genomes, latest_genomes, species, genome_id, genbank_stats)
        sync_latest_genomes(local_genomes, latest_genomes, latest_genome_dirs, species, genome_id, genbank_stats)

def log_error(msg, species, genbank_stats, ymdt):
    if msg == "latest":
        msg = "no latest_assembly_versions dir for"
        with open(genbank_stats, "a") as stats:
            stats.write("{} - {} {}\n".format(species, msg, ymdt))

def main():

    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    ymdt = strftime("%y.%m.%d_%H%M")
    genbank_stats = os.path.join(genbank_mirror, ".info", "genbank_stats_{}.txt".format(ymdt))
    check_dirs(genbank_mirror)
    complete_species_list = ftp_complete_species_list()
    latest_assembly_versions = get_latest_assembly_versions(genbank_mirror, complete_species_list, genbank_stats, ymdt)
#   latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
#   latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None)
#   latest_assembly_versions.columns = ["id"]
#   grab_and_organize_genomes(genbank_mirror, genbank_stats, latest_assembly_versions)
#   rename(genbank_mirror, assembly_summary)

main()
