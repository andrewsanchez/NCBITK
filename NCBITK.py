#!/usr/bin/env python3

import os, rename_fastas, shutil
from subprocess import call
from ftplib import FTP, error_temp
from time import strftime

def get_assembly_summary(wget, local_mirror, assembly_summary='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'):

    """Get current version of assembly_summary.txt"""

    assembly_summary_destination = os.path.join(local_mirror, "assembly_summary.txt")
    if os.path.isfile(assembly_summary_destination):
        print("Removing {}".format(assembly_summary_destination))
        os.remove(assembly_summary_destination)

    if wget:
        print("Fetching current version of assembly_summary.txt.\n")
        call(['wget', '-nv', '-O', assembly_summary_destination, assembly_summary])
    else:
        print("assembly_summary.txt will not be downloaded.\n")

def clean_up_files_and_dirs(local_mirror, assembly_summary_df):

    # Remove fastas no longer in assembly_summary.txt
    for root, dirs, files in os.walk(local_mirror):
        for name in files:
            if name.startswith("GCA"):
                accession_id = name.split("_")[0:2]
                accession_id = "_".join(accession_id)
                if accession_id not in assembly_summary_df.index:
                    print("{} no longer in assembly_summary.txt and will be removed.".format(name))
                    os.remove(os.path.join(root, name))

        # remove empty directories
        for name in dirs:
            name = os.path.join(root, name)
            if not os.listdir(name):
                os.rmdir(name)

    complete_ftp_paths = os.path.join(local_mirror, "fasta_list.txt")
    if os.path.isfile(complete_ftp_paths):
        os.remove(complete_ftp_paths)

    species_and_accession_ids = os.path.join(local_mirror, "species_and_accession_ids.csv")
    if os.path.isfile(species_and_accession_ids):
        print("Removing {}".format(species_and_accession_ids))
        os.remove(species_and_accession_ids)

    # Remove old filter_file
    filter_files_dir = os.path.join(local_mirror, "filter_files")
    for root, dirs, files in os.walk(filter_files_dir):
        for f in files:
            os.remove(os.path.join(root, f))

def assembly_summary_to_df(local_mirror):

    import pandas as pd

    path_to_summary = os.path.join(local_mirror, "assembly_summary.txt")
    assembly_summary_df = pd.read_csv(path_to_summary, sep="\t", index_col=0, skiprows=1)

    return assembly_summary_df

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    print("Logging into ftp.ncbi.nlm.nih.gov/{}/".format(directory))
    ftp.login()
    ftp.cwd(directory)

    return ftp

def ftp_complete_species_list(local_mirror):

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    print("Getting list of all bacteria directories.")
    print("Estimated wait time:  1 minute.\n")
    ftp = ftp_login()
    complete_species_list = ftp.nlst()

    ftp_stats = os.path.join(local_mirror, "ftp_stats.txt")
    time_stamp = strftime("%y/%m/%d_%H:%M")
    with open(ftp_stats, "a") as stats:
        stats.write("{} - {} dirs at genbank/bacteria/".format(time_stamp, len(complete_species_list)))

    return complete_species_list

def get_accessions_in_latest_dirs(local_mirror, complete_species_list, ftp):

    ftp_stats = os.path.join(local_mirror, "ftp_stats.txt")
    species_and_accession_ids = os.path.join(local_mirror, "species_and_accession_ids.csv")
    print("Writing species and their accession ids to {}".format(species_and_accession_ids))
    for species in complete_species_list:
        latest = os.path.join(species, "latest_assembly_versions")
        try:
            accessions_in_latest_dirs = [accession.split("/")[-1] for accession in ftp.nlst(latest)]
            with open(species_and_accession_ids, "a") as f:
                for accession in accessions_in_latest_dirs:
                    f.write("{},{}\n".format(species, accession))
                print("Got accession ids for {}.".format(species))
        except error_temp:
            print("No latest_assembly_versions dir for {}".format(species))
            with open(ftp_stats, "a") as stats:
                stats.write("{} - No latest_assembly_versions dir:  {}".format(strftime("%y/%m/%d"), species))

def get_dir_structure(local_mirror):

    complete_species_list = ftp_complete_species_list(local_mirror)
    species_and_accession_ids = os.path.join(local_mirror, "species_and_accession_ids.csv")
    ftp = ftp_login()
    try:
        get_accessions_in_latest_dirs(local_mirror, complete_species_list, ftp)
    except BrokenPipeError:
        ftp = ftp_login()
        get_accessions_in_latest_dirs(local_mirror, complete_species_list, ftp)

def mk_dir_structure(local_mirror, assembly_summary_df):

    import gzip

    renamed_dir = "{}_renamed".format(local_mirror)
    species_and_accession_ids = os.path.join(local_mirror, "species_and_accession_ids.csv")

    with open(species_and_accession_ids) as f:
        for line in f:
            species = line.split(",")[0]
            accession = line.split(",")[1].strip()
            assembly_id = "_".join(accession.split("_")[0:2])
            if assembly_id in assembly_summary_df.index:
                species_dir = os.path.join(renamed_dir, species)
                print("Updating {}\n".format(species_dir))
                source = os.path.join(local_mirror, accession, "{}_genomic.fna.gz".format(accession))
                destination = os.path.join(species_dir, "{}.fasta".format(accession))

                if not os.path.isdir(species_dir):
                    print("Creating {}.\n".format(species_dir))
                    os.mkdir(species_dir)

                if not os.path.isfile(destination):
                    zipped = gzip.open(source)
                    unzipped = open(destination, "wb")
                    decoded = zipped.read()
                    unzipped.write(decoded)
                    print("Successfully unzipped {}".format(destination))
                    zipped.close()
                    unzipped.close()

def ftp_paths_from_assembly_summary(local_mirror, assembly_summary_df):

    """Write full ftp paths for fastas to file 'local_mirror/fasta_list.txt'
    for use by rsync's --files-from=fasta_list"""

    ftp_paths = os.path.join(local_mirror, "ftp_paths_from_assembly_summary.txt")

    with open(ftp_paths, "a") as f:
        for ftp_path in assembly_summary_df.ftp_path:
            accession = ftp_path.split("/")[-1]
            fasta = "{}_genomic.fna.gz".format(accession)
            full_fasta_path = os.path.join(accession, fasta)
            f.write(full_fasta_path+"\n")

    return ftp_paths

def get_species_list_from_file(source):

    """Read organism list from file"""

    print("Getting list of directories to sync with from specified file.\n")

    species_list = [] 
    with open(source) as f:
        for line in f:
            line = line.strip()
            species_list.append(line)

    return species_list

def get_species_from_list(species_and_or_genera):

    """Generate species_list from list of species and/or genera"""

    ftp = ftp_login()
    complete_species_list = ftp_complete_species_list(local_mirror)
    species_list = []

    for i in species_and_or_genera:
        for species in complete_species_list:
            if species.startswith(i):
                species_list.append(species)

    return species_list

def check_dirs(local_mirror):

    """Make directories to store fastas if they don't already exist."""

    if not os.path.isdir(local_mirror):
        os.mkdir(local_mirror)

    renamed_dir =  "{}_renamed".format(local_mirror)
    if not os.path.isdir(renamed_dir):
        os.mkdir(renamed_dir)

    filter_files_dir = os.path.join(local_mirror, "filter_files")
    if not os.path.isdir(filter_files_dir):
        os.mkdir(filter_files_dir)

def rsync_latest_fastas_from_assembly_summary(local_mirror, fasta_list):

    """
    Get all latest fastas latest assembly_summary.txt
    """

    rsync_log = os.path.join(local_mirror, "rsync_log_{}.txt".format(strftime("%y-%m-%d_%H_%M")))

    call(['rsync',
        '--chmod=ugo=rwX', # Change permissions so files can be copied/renamed
        '--times',
        '--progress',
        '--itemize-changes',
        '--stats',
        '--files-from='+fasta_list,
        '--log-file='+rsync_log,
        'ftp.ncbi.nlm.nih.gov::genomes/all/',
        local_mirror])

def write_filter_list(local_mirror, organism, ftp):
    
    """Create accepted files list for each organism for rsync filter."""

    print("Creating accepted files list for {}.\n".format(organism))

    filter_files_dir = os.path.join(local_mirror, "filter_files")
    filter_file = os.path.join(filter_files_dir, organism+".txt")
    latest = os.path.join(organism, 'latest_assembly_versions')
    try:
        full_paths = ftp.nlst(latest)
    except error_temp:
        print("{} doesn't have a latest_assembly_versions/ directory and will be skipped".format(organism))
    with open(filter_file, "a") as f:
        for full_path in full_paths:
            accession = full_path.split("/")[-1]
            fasta = "_".join([accession, "genomic.fna.gz"])
            fasta = "/".join([accession, fasta])
            f.write(fasta+"\n")

    return filter_file

def rsync_latest_fastas_from_list(local_mirror, organism_list):
    
    renamed_dir = "{}_renamed".format(local_mirror)
    ftp = ftp_login()

    for organism in organism_list:
        try:
            filter_file = write_filter_list(local_mirror, organism, ftp)
        except BrokenPipeError:
            ftp = ftp_login()
            filter_file = write_filter_list(local_mirror, organism, ftp)

        rsync_log = os.path.join(local_mirror, organism, "rsync_log.txt")
        latest = os.path.join(organism, 'latest_assembly_versions')
        organism_dir = os.path.join(local_mirror, organism)
        call(['rsync',
            '--chmod=ugo=rwX', # Change permissions so files can be copied/renamed
            '--copy-links', # Follow symbolic links
            '--recursive',
            '--times',
            '--prune-empty-dirs',
            '--progress',
            '--itemize-changes',
            '--stats',
            '--files-from='+filter_file,
            '--log-file='+rsync_log,
            'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/'+latest,
            organism_dir])

        move_and_unzip(local_mirror, organism)
        organism_index = organism_list.index(organism)+1
        print("{} out of {} organisms updated.\n".format(organism_index, len(organism_list)))

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

def move_and_unzip(local_mirror, organism):

    renamed_dir = "{}_renamed".format(local_mirror)
    species_dir = os.path.join(renamed_dir, organism)

    if os.path.isdir(species_dir):
        source = os.path.join(local_mirror, organism)
        cp_files(source, species_dir)
        gunzip(species_dir)

    else:
        os.mkdir(species_dir)
        source = os.path.join(local_mirror, organism)
        cp_files(source, species_dir)
        gunzip(species_dir)

    print("\nFiles renamed to:  ")
    rename_fastas.rename(species_dir)

def Main():

    import argparse 

    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument('local_mirror', help = 'Directory to save fastas', type=str)
    parser.add_argument('-W', '--no-wget', help = "Skip downloading current assembly_summary.txt", action='store_false')
    parser.add_argument('-f', '--from_file', help = 'Full path to file containing directories to sync with.')
    parser.add_argument('-l', '--from_list', help = 'List of species and/or genera to download.  Must exactly match the name of \
            a directory at ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/ OR a genus, and be separated by spaces.  \
            e.g.: "Clostridium, Escherichia, Bacillus_anthracis"', nargs="+")
    parser.add_argument('-r', "--reorganize", action="store_true")
    args = parser.parse_args()

    local_mirror = args.local_mirror.strip("/")
    get_assembly_summary(args.no_wget, local_mirror)
    check_dirs(local_mirror)
    assembly_summary_df = assembly_summary_to_df(local_mirror)
    clean_up_files_and_dirs(local_mirror, assembly_summary_df)

    if args.from_file:
        organism_list = get_species_list_from_file(args.from_file)
        rsync_latest_fastas_from_list(local_mirror, organism_list)
    elif args.from_list:
        organism_list = get_species_from_list(args.from_list)
        rsync_latest_fastas_from_list(local_mirror, organism_list)
    else:
        fasta_list = ftp_paths_from_assembly_summary(local_mirror, assembly_summary_df)
        rsync_latest_fastas_from_assembly_summary(local_mirror, fasta_list)
        get_dir_structure(local_mirror)
        mk_dir_structure(local_mirror, assembly_summary_df)

Main()
