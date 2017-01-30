#!/usr/bin/env python

import os, argparse, gzip
import pandas as pd
import tarfile
from urllib.request import urlretrieve
from urllib.error import URLError
from ftplib import FTP, error_temp
from time import strftime, sleep
from glob import glob
from re import sub

def get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")
    urlretrieve(assembly_summary_url, assembly_summary_dst)
    assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

def get_species_names(genbank_mirror, taxdump_url="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"):

    """
    Get names.dmp from the taxonomy dump
    """

    info_dir = os.path.join(genbank_mirror, ".info")
    taxdump = urlretrieve(taxdump_url)
    taxdump_tar = tarfile.open(taxdump[0])
    taxdump_tar.extract('names.dmp', info_dir)
    names_dmp = os.path.join(genbank_mirror, ".info", 'names.dmp')
    names = pd.read_csv(names_dmp, sep='\t|', skiprows=3, index_col=0, header=None, usecols=[0,2,6], engine='python', iterator=True)
    names = pd.concat(chunk[chunk[6] == 'scientific name'] for chunk in names)
    names.drop(6, axis=1, inplace=True)
    names.columns = ['species']
    names.index.name = 'species_taxid'
    names.species.replace({' ': '_'}, regex=True, inplace=True)
    names.species.replace({'/': '_'}, regex=True, inplace=True)
    names.to_csv(names_dmp)

    return names

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

    slurm = os.path.join(genbank_mirror, ".info", "slurm")
    if not os.path.isdir(slurm):
        os.mkdir(slurm)

def ftp_login(directory="genomes/genbank/bacteria", email="aas229@nau.edu"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site, user="anonymous", passwd=email)
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
            genome_path = item[1]
            accession_id = item[0]
            f.write("{},{},{}\n".format(species, genome_path, accession_id))

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

def grab_zipped_genome(genbank_mirror, species, genome_id, genome_path, ext=".fna.gz"):

    """
    Download compressed genome from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
    """

    zipped = "{}_genomic{}".format(genome_path, ext)
    zipped_dst = "{}_genomic{}".format(genome_id, ext)
    zipped_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{}/{}".format(genome_path, zipped)
    zipped_dst = os.path.join(genbank_mirror, species, zipped_dst)
    urlretrieve(zipped_url, zipped_dst)

def unzip_genome(root, f, genome_id):

    """
    Decompress genome and remove the compressed genome.
    """

    zipped_src = os.path.join(root, f)
    zipped = gzip.open(zipped_src)
    decoded = zipped.read()
    unzipped = "{}.fasta".format(genome_id)
    unzipped = os.path.join(root, unzipped)
    unzipped = open(unzipped, "wb")
    unzipped.write(decoded)
    zipped.close()
    unzipped.close()
    os.remove(zipped_src)

def unzip_genbank_mirror(genbank_mirror):

    for root, files, dirs, in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("gz"):
                genome_id = "_".join(f.split("_")[:2])
                unzip_genome(root, zipped_src, genome_id)

def check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids):

    for genome_id in local_genome_ids:
        if genome_id not in latest_genome_ids:
            fasta = glob("{}*".format(genome_id))
            os.remove(os.path.join(genbank_mirror, species, fasta[0]))

def sync_latest_genomes(genbank_mirror, assembly_summary, names):

        # local_genome_ids = get_local_genome_ids(species_dir)
        # if genome_id not in local_genome_ids:
        #     try:
        #         grab_zipped_genome(genbank_mirror, species, genome_id, genome_path)
        #     except URLError:
        #         grab_zipped_genome(genbank_mirror, species, genome_id, genome_path, ext=".fasta.gz")
        #     except URLError:
        #         with open(genbank_stats, "a") as stats:
        #             stats.write("URLError for {}\n".format(genome_id))
        #     with open(genbank_stats, "a") as stats:
        #         stats.write("{} downloaded\n".format(genome_id))
    for genome in assembly_summary.index:
        url = assembly_summary.ftp_path[genome]
        taxid = assembly_summary.species_taxid.loc[genome]
        species = names.species.loc[taxid]
        dst = os.path.join(genbank_mirror, species, genome)
        urlretrieve(url, dst)
        print("Downloaded {} to {}".format(genome, dst))

def get_local_genome_ids(species):

        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]

        return local_genome_ids

def grab_and_organize_genomes(genbank_mirror, genbank_stats, latest_assembly_versions):

    species_directories = os.listdir(genbank_mirror)

    for species in species_directories:
        local_genome_ids = get_local_genome_ids(species_dir)
        latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        latest_genome_paths = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)
        print("species", species)
        print("local_genome_ids", local_genome_ids)
        print("latest_genome_ids", latest_genome_ids)
        print("latest_genome_paths", latest_genome_paths)

        check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids, genbank_stats)
        sync_latest_genomes(genbank_mirror, species, local_genome_ids, ids_and_paths, genbank_stats)

def log_error(msg, species, genbank_stats, ymdt):
    if msg == "latest":
        msg = "no latest_assembly_versions dir for"
        with open(genbank_stats, "a") as stats:
            stats.write("{} - {} {}\n".format(species, msg, ymdt))

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-u", "--unzip", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror

    assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
    unzip_genbank_mirror(genbank_mirror)
    rename(genbank_mirror, assembly_summary)

if __name__ == "__main__":
    main()
