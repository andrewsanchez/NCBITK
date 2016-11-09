#!/usr/bin/env python

import os, rename_fastas, shutil, argparse
import pandas as pd
from urllib.request import urlretrieve
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

def clean_up_files_and_dirs(genbank_mirror, assembly_summary):

    ncbitk_log = os.path.join(genbank_mirror, "ncbitk_log.txt")
    renamed_dir = "{}_renamed".format(genbank_mirror)

    # Remove fastas from genbank_mirror no longer in assembly_summary.txt
    with open(ncbitk_log, "a") as f:
        for root, dirs, files in os.walk(genbank_mirror):
            for name in files:
                if name.startswith("GCA"):
                    accession_id = name.split("_")[0:2]
                    accession_id = "_".join(accession_id)
                    if accession_id not in assembly_summary_df.index:
                        os.remove(os.path.join(root, name))
                        f.write("{} removed, {}\n".format(name, time_stamp))

        # remove empty directories
        for name in dirs:
            name = os.path.join(root, name)
            if not os.listdir(name): # rm empty dirs
                f.write("{} removed, {}\n".format(name, time_stamp))
                os.rmdir(name)

    # Remove fastas from renamed_dir no longer in assembly_summary.txt
    for root, dirs, files in os.walk(renamed_dir):
        for name in files:
            if name.startswith("GCA"):
                accession_id = name.split("_")[0:2]
                accession_id = "_".join(accession_id)
                if accession_id not in assembly_summary_df.index:
                    os.remove(os.path.join(root, name))
                    with open(ncbitk_log, "a") as f:
                        f.write("{} removed, {}\n".format(name, time_stamp))

        # remove empty directories
        for name in dirs:
            name = os.path.join(root, name)
            if not os.listdir(name):
                os.rmdir(name)

    complete_ftp_paths = os.path.join(genbank_mirror, "fasta_list.txt")
    if os.path.isfile(complete_ftp_paths):
        os.remove(complete_ftp_paths)

    ftp_paths = os.path.join(genbank_mirror, "ftp_paths_from_assembly_summary.txt")
    if os.path.isfile(ftp_paths):
        os.remove(ftp_paths)

    species_and_accession_ids = os.path.join(genbank_mirror, "species_and_accession_ids.csv")
    if os.path.isfile(species_and_accession_ids):
        print("Removing {}".format(species_and_accession_ids))
        os.remove(species_and_accession_ids)

    # Remove old filter_file
    filter_files_dir = os.path.join(genbank_mirror, "filter_files")
    for root, dirs, files in os.walk(filter_files_dir):
        for f in files:
            os.remove(os.path.join(root, f))

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
            f.write("{},{}".format(species, accession))

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

    def log_error():
            with open(genbank_stats, "a") as stats:
                stats.write("{} - {} {}\n".format(species, no_latest_msg, ymd))

    latest_assembly_versions = pd.read_csv(species_and_accession_ids)
    latest_assembly_versions.columns = ["species", "id"]

    return latest_assembly_versions


def generate_get_tree_commands(complete_species_list):

    """
    Generate commands for use in slurm array to recreate genbank directory structure
    """

    ftp = ftp_login()
    with open(os.path.join(genbank_mirror, "get_dir_structure_commands.txt"), "a") as f:
        for species in complete_species_list:
            pass

def mk_dir_structure(genbank_mirror, assembly_summary_df):

    """
    Write decompressed FASTA's to renamed_dir and rename files.
    """

    import gzip

    renamed_dir = "{}_renamed".format(genbank_mirror)
    species_and_accession_ids = os.path.join(genbank_mirror, "species_and_accession_ids.csv")

    with open(species_and_accession_ids) as f:
        for line in f:
            complete_id = line.split(",")[1].strip()
            assembly_id = "_".join(complete_id.split("_")[0:2])
            if assembly_id in assembly_summary_df.index:
                species = line.split(",")[0]
                species_dir = os.path.join(renamed_dir, species)
                if not os.path.isdir(species_dir):
                    print("Creating {}.\n".format(species_dir))
                    os.mkdir(species_dir)

                source = os.path.join(genbank_mirror, complete_id, "{}_genomic.fna.gz".format(complete_id))
                destination = os.path.join(species_dir, "{}.fasta".format(complete_id))
                if not os.path.isfile(destination):
                    print("Updating {}\n".format(species_dir))
                    zipped = gzip.open(source)
                    unzipped = open(destination, "wb")
                    decoded = zipped.read()
                    unzipped.write(decoded)
                    print("Successfully unzipped {}".format(destination))
                    zipped.close()
                    unzipped.close()

    rename_fastas.rename(renamed_dir, assembly_summary_df)

def organize_and_get_genomes(latest_assembly_versions):

    for genome in latest_assembly_versions["id"]:
    None

def get_species_from_list(species_and_or_genera):

    """Generate species_list from list of species and/or genera"""

    ftp = ftp_login()
    complete_species_list = ftp_complete_species_list(genbank_mirror)
    species_list = []

    for i in species_and_or_genera:
        for species in complete_species_list:
            if species.startswith(i):
                species_list.append(species)

    return species_list

def cp_files(source, destination):

    """Copy new files to destination."""

    for root, dirs, files in os.walk(source):
        for f in files:
            try:
                source = os.path.join(root, f)
                copy_to = os.path.join(destination, f)
                shutil.copyfile(source, copy_to)
            except shutil.SameFileError: # Skip files that already exist in destination
                continue

def gunzip(target_dir):

    import gzip

    """Unzip fastas"""

    for root, dirs, files in os.walk(target_dir):
        for f in files:
            if f.endswith(".gz"):
                f = os.path.join(root, f)
                destination = os.path.splitext(f)[0]
                zipped = gzip.open(f)
                unzipped = open(destination, 'wb')
                decoded = zipped.read()
                unzipped.write(decoded)
                zipped.close()
                unzipped.close()
                os.remove(f)

def move_and_unzip(genbank_mirror, organism, assembly_summary_df):

    renamed_dir = "{}_renamed".format(genbank_mirror)
    species_dir = os.path.join(renamed_dir, organism)

    if os.path.isdir(species_dir):
        source = os.path.join(genbank_mirror, organism)
        cp_files(source, species_dir)
        gunzip(species_dir)

    else:
        os.mkdir(species_dir)
        source = os.path.join(genbank_mirror, organism)
        cp_files(source, species_dir)
        gunzip(species_dir)

    print("\nFiles renamed to:  ")
    rename_fastas.rename(species_dir, assembly_summary_df)

def main():

    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror.strip("/")
    ymdt = strftime("%y.%m.%d_%H%M")
    assembly_summary_location = os.path.join(genbank_mirror, "assembly_summary.txt")
    assembly_summary_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
    check_dirs(genbank_mirror)
    complete_species_list = ftp_complete_species_list()
    assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_location, assembly_summary_url)
    latest_assembly_versions = get_latest_assembly_versions(genbank_mirror, complete_species_list, ymdt)
 #  clean_up_files_and_dirs(genbank_mirror, assembly_summary_df)

 #  if args.rename:
 #      if args.rename_target:
 #          rename_target = os.path.join(renamed_dir, args.rename_target)
 #          rename_fastas.rename(rename_target, assembly_summary_df)
 #      else:
 #          complete_species_list = ftp_complete_species_list(genbank_mirror)
 #          get_accessions_in_latest_dirs(genbank_mirror, complete_species_list)
 #          mk_dir_structure(genbank_mirror, assembly_summary_df)
 #  else:
 #      ftp_paths = ftp_paths_from_assembly_summary(genbank_mirror, assembly_summary_df)
 #      complete_species_list = ftp_complete_species_list(genbank_mirror)
 #      dir_structure = get_dir_structure(genbank_mirror, complete_species_list)
 #      mk_dir_structure(genbank_mirror, assembly_summary_df)

main()
