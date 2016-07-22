#!/usr/bin/env python3

import os, subprocess, rename_fastas, time, shutil
from ftplib import FTP, error_temp

def get_assembly_summary(wget, destination, summary_file='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'):

    """Get current version of assembly_summary.txt"""

    assembly_summary = summary_file
    assembly_summary_destination = os.path.join(destination, "assembly_summary.txt")
    if wget:
        print("Fetching current version of assembly_summary.txt.\n")
        try:
            os.remove('assembly_summary.txt')
            print("Removing old assembly_summary.txt")
            subprocess.call(['wget', '-nv', '-O', assembly_summary_destination, assembly_summary])
        except OSError:
            subprocess.call(['wget', '-nv', '-O', assembly_summary_destination, assembly_summary])
    else:
        print("assembly_summary.txt will not be downloaded.\n")

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd(directory)

    return ftp

def ftp_complete_species_list():

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    print("Getting list of all bacteria directories.")
    print("Estimated wait time:  1 minute.\n")
    ftp = ftp_login()
    complete_species_list = ftp.nlst()
    with open(ftp_stats, "a") as stats:
        stats.write("{} dirs at genbank/bacteria/.".format(len(complete_species_list)))

    return complete_species_list

def get_accessions_in_latest_dirs(complete_species_list, species_dict, ftp):

    for species in complete_species_list:
        latest = os.path.join(species, "latest_assembly_versions")
        try:
            accessions_in_latest_dirs = [accession.split("/")[-1] for acccession in ftp.nlst(latest)]
            species_dict[species] = accessions_in_latest_dirs
        except error_temp:
            with open(ftp_stats, "a") as stats:
                stats.write("No latest_assembly_versions dir:  {}".format(species))

def get_dir_structure():

    ftp = ftp_login()
    complete_species_list = ftp.nlst()
    species_and_ids = {}

    try:
        get_accessions_in_latest_dirs(complete_species_list, species_and_ids, ftp)
    except BrokenPipeError:
        ftp = ftp_login()
        get_accessions_in_latest_dirs(complete_species_list, species_and_ids, ftp)

    return species_and_ids

def mk_dir_structure(local_mirror):

    import pandas as pd

    renamed_dir = "{}_renamed".format(local_mirror.strip("/"))
    for root, dirs, files in os.walk(renamed_dir):
        for fasta in files:
            if fasta.startswith("GCA"):
                accession_id = fasta.split("_")[0:2]
                accession_id = "_".join(accession_id)
                if accession_id not in assembly_summary_df.index:
                    os.remove(os.path.join(root, fasta))

    ftp = ftp_login()
    complete_species_list = ftp.nlst()

    path_to_summary = os.path.join(local_mirror, "assembly_summary.txt")
    assembly_summary_df = pd.read_csv(path_to_summary, sep="\t", index_col=0, skiprows=1)

    species_and_accession_ids = get_dir_structure()

    for species, accessions in species_and_accession_ids:
        species_dir = os.path.join(renamed_dir, species)
        print("Updating {}\n".format(species_dir))

        if not os.path.isdir(species_dir):
            print("Creating {}.\n".format(species_dir))
            os.mkdir(species_dir)
        for accession in accessions:
            if accession in assembly_summary_df.index:
                source = os.path.join(local_mirror, accession, "{}_genomic.fna.gz".format(accession))
                destination = os.path.join(species_dir, "{}.fasta".format(accession))
                if not os.path.isfile(destination):
                    zipped = gzip.open(source)
                    unzipped = open(destination, "wb")
                    decoded = zipped.read()
                    unzipped.write(decoded)
                    print("Successfully unzipped {}".format(destination)
                    zipped.close()
                    unzipped.close()

  # for species in complete_species_list:
  #     species_index = complete_species_list.index(species)+1
  #     print("Updating {} out of {} directories.\n".format(species_index, len(complete_species_list)))

  #     species_dir = os.path.join(renamed_dir, species)
  #     latest = os.path.join(species, "latest_assembly_versions")

  #     try:
  #         accessions_in_latest_dirs = ftp.nlst(latest)
  #         # Included in try block to avoid creating empty dirs
  #         if not os.path.isdir(species_dir):
  #             os.mkdir(species_dir)
  #             print("Creating {}.\n".format(species_dir))

  #         for accession in accessions_in_latest_dirs:
  #             accession_id = accession.split("/")[-1]
  #             try:
  #                 source_fasta = os.path.join(local_mirror, accession_id, "{}_genomic.fna.gz".format(accession_id))
  #                 destination = os.path.join(species_dir, "{}_genomic.fna.gz".format(accession_id))
  #                 shutil.copyfile(source_fasta, destination)
  #                 print("Copied {} to {}\n.".format(source_fasta, destination))
  #             except shutil.SameFileError: # Skip files that already exist in destination
  #                 continue
  #             except FileNotFoundError:
  #                 missing_files = os.path.join(local_mirror, "missing_files.txt")
  #                 with open(missing_files, "a") as f:
  #                     f.write(accession)

  #     except BrokenPipeError: # just in case our ftp connection times out
  #         ftp = ftp_login()
  #         accessions_in_latest_dirs = ftp.nlst(latest)

  #         if not os.path.isdir(species_dir):
  #             os.mkdir(species_dir)
  #             print("Creating {}.\n".format(species_dir))

  #         for accession in accessions_in_latest_dirs:
  #             accession_id = accession.split("/")[-1]
  #             try:
  #                 source_fasta = os.path.join(local_mirror, accession_id, "{}_genomic.fna.gz".format(accession_id))
  #                 destination = os.path.join(species_dir, "{}_genomic.fna.gz".format(accession_id))
  #                 shutil.copyfile(source_fasta, destination)
  #                 print("Copied {} to {}\n.".format(source_fasta, destination))
  #             except shutil.SameFileError: # Skip files that already exist in destination
  #                 continue
  #             except FileNotFoundError:
  #                 missing_files = os.path.join(local_mirror, "missing_files.txt")
  #                 with open(missing_files, "a") as f:
  #                     f.write(accession)
  #     except error_temp: # skip organisms with no latest_assembly_dir
  #         continue

def ftp_paths_from_assembly_summary(local_mirror):

    """Write full ftp paths for fastas to file 'local_mirror/fasta_list.txt'
    for use by rsync's --files-from=fasta_list"""

    import pandas as pd

    path_to_summary = os.path.join(local_mirror, "assembly_summary.txt")
    assembly_summary_df = pd.read_csv(path_to_summary, sep="\t", index_col=0, skiprows=1)

    complete_ftp_paths = os.path.join(local_mirror, "fasta_list.txt")

    if os.path.isfile(complete_ftp_paths):
        os.remove(complete_ftp_paths)

    with open(complete_ftp_paths, "a") as f:
        for ftp_path in assembly_summary_df.ftp_path:
            accession = ftp_path.split("/")[-1]
            fasta = "_".join([accession, "genomic.fna.gz"])
            full_fasta_path = os.path.join(accession, fasta)
            f.write(full_fasta_path+"\n")

    return complete_ftp_paths

def get_species_list_from_file(source):

    """Read organism list from file"""

    print("Getting list of directories to sync with from specified file.\n")

    organism_list = [] 
    with open(source) as f:
        for line in f:
            line = line.strip()
            organism_list.append(line)

    return organism_list

def get_species_list_from_genus(genera):

    """Generate organism list from a gunus or list of genera"""

    ftp = ftp_login()
    complete_species_list = ftp_complete_species_list()
    species_list_genera = []

    for genus in genera:
        for organism in complete_species_list:
            if organism.startswith(genus):
                species_list_genera.append(organism)

    return species_list

def check_dirs(local_mirror):

    """Make directories to store fastas if they don't already exist."""

    if not os.path.isdir(local_mirror):
        os.mkdir(local_mirror)

    renamed_fastas_dir =  "{}_renamed".format(local_mirror.strip("/"))
    if not os.path.isdir(renamed_fastas_dir):
        os.mkdir(renamed_fastas_dir)

def write_filter_list(local_mirror, organism, ftp):
    
    """Create accepted files list for each organism for rsync filter."""

    filter_files_dir = os.path.join(local_mirror, "filter_files")

    if not os.path.isdir(filter_files_dir):
        os.mkdir(filter_files_dir)
        print("Created {} to store filter lists.\n".format(filter_files_dir))

    print("Creating accepted files list for {}.\n".format(organism))

    filter_file = os.path.join(filter_files_dir, organism+".txt")

    if os.path.isfile(filter_file):
        os.remove(filter_file)

    latest = os.path.join(organism, 'latest_assembly_versions')
    full_paths = ftp.nlst(latest)

    with open(filter_file, "a") as f:
        for full_path in full_paths:
            accession = full_path.split("/")[-1]
            fasta = "_".join([accession, "genomic.fna.gz"])
            fasta = "/".join([accession, fasta])
            f.write(fasta+"\n")

    return filter_file

def get_latest_fastas_from_assembly_summary(local_mirror, fasta_list):

    import pandas as pd

    """
    Get all latest fastas and remove deprecated ones based
    on ftp paths generated from latest assembly_summary.txt
    """

    rsync_log = os.path.join(local_mirror, "rsync_log.txt")

    path_to_summary = os.path.join(local_mirror, "assembly_summary.txt")
    assembly_summary_df = pd.read_csv(path_to_summary, sep="\t", index_col=0, skiprows=1)

    # Remove files no longer in assembly_summary.txt
    for root, dirs, files in os.walk(local_mirror):
        for fasta in files:
            if fasta.startswith("GCA"):
                accession_id = fasta.split("_")[0:2]
                accession_id = "_".join(accession_id)
                if accession_id not in assembly_summary_df.index:
                    print("{} is no longer assembly_summary.txt and will be removed.".format(fasta))
                    os.remove(os.path.join(root, fasta))

    subprocess.call(['rsync',
        '--chmod=ugo=rwX', # Change permissions so files can be copied/renamed
        '--times',
        '--progress',
        '--itemize-changes',
        '--stats',
        '--files-from='+fasta_list,
        '--log-file='+rsync_log,
        'ftp.ncbi.nlm.nih.gov::genomes/all/',
        local_mirror])

def get_latest_fastas_from_list(local_mirror, organism_list):
    
    renamed_fasta_dir = "%s_renamed" % (local_mirror.strip("/"))
    ftp = ftp_login()

    for organism in organism_list:
        try:
            filter_file = write_filter_list(local_mirror, organism, ftp)
        except BrokenPipeError:
            ftp = ftp_login()
            filter_file = write_filter_list(local_mirror, organism, ftp)
        except error_temp:
            print("{} doesn't have a latest_assembly_versions/ directory and will be skipped".format(organism))
            continue

        rsync_log = os.path.join(local_mirror, organism, "rsync_log.txt")
        latest = os.path.join(organism, 'latest_assembly_versions')
        renamed_organism_dir = os.path.join(renamed_fasta_dir, organism)
        organism_dir = os.path.join(local_mirror, organism)
        subprocess.call(['rsync',
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

    renamed_fasta_dir = "_".join([local_mirror.strip("/"),"renamed"])
    renamed_organism_dir = os.path.join(renamed_fasta_dir, organism)

    if os.path.isdir(renamed_organism_dir):
        source = os.path.join(local_mirror, organism)
        cp_files(source, renamed_organism_dir)
        gunzip(renamed_organism_dir)

    else:
        os.mkdir(renamed_organism_dir)
        source = os.path.join(local_mirror, organism)
        cp_files(source, renamed_organism_dir)
        gunzip(renamed_organism_dir)

    print("\nFiles renamed to:  ")
    rename_fastas.rename(renamed_organism_dir)

def new_and_removed_log(local_mirror):

    """Write a simple log file indicating new files and files that were removed."""

    time_stamp = time.strftime("%y/%m/%d_%H:%M")
    new = os.path.join(local_mirror, "{}_new.txt.format(time_stamp)")
    removed = os.path.join(local_mirror, "{}_removed.txt.format(time_stamp)")
    path_to_summary = os.path.join(local_mirror, "assembly_summary.txt")
    assembly_summary_df = pd.read_csv(path_to_summary, sep="\t", index_col=0, skiprows=1)

    if os.path.isdir(local_mirror):
        current_fasta_count = len(os.listdir(local_mirror))
        current_fasta_collection = os.listdir(local_mirror)
        new_fastas = []
        for fasta in current_fasta_collection:
            if fasta not in assembly_summary_df.index:
                new_fastas.append(fasta)

    with open(new, "a") as log:
        log.write("start time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")

    with open(removed, "a") as log:
        log.write("start time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")

def monitor_changes(local_mirror, organism_list):
    renamed_fastas_dir = local_mirror.strip("/") + "_renamed"
    dirs = get_species_list(organism_list)
    countdirs = str(len(dirs))
    local_dir_count = str(len(os.listdir(renamed_fastas_dir)))
    missingdirs = int(countdirs) - int(local_dir_count)
    changes_log = os.path.join(local_mirror, "changes_log.txt")
    with open(changes_log, "a") as log:
        log.write("Finish time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")
        log.write('Dirs at ftp.ncbi:  ' + countdirs)
        log.write("\n")
        log.write('Dirs in {}:  {}\n'.format(renamed_fastas_dir, local_dir_count))
        log.write('Missing dirs = ' + str(missingdirs))
        log.write("\n")
        for i in dirs:
            if i not in os.listdir(local_mirror):
                log.write(i+"\n")

def Main():

import argparse 

    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument('local_mirror', help = 'Directory to save fastas', type=str)
    parser.add_argument('-W', '--no-wget', help = "Skip downloading current assembly_summary.txt", action='store_false')
    parser.add_argument('-f', '--from_file', help = 'Full path to file containing directories to sync with.')
    parser.add_argument('-l', '--from_list', help = 'List of organisms to download.  Must exactly match the name of \
            a directory at ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/ and be separated by spaces.  \
            Use -g to get an enetire genus without having to list every species.', nargs="+")
    parser.add_argument('-g', '--genus', help = 'Retrieve all species of a genus or a list of genera eg. "Clostridium"')
    parser.add_argument('-r', "--reorganize", action="store_true")
    args = parser.parse_args()

    if args.from_file:
        check_dirs(args.local_mirror)
        get_assembly_summary(args.no_wget, args.local_mirror)
        organism_list = get_species_list_from_file(args.from_file)
        get_latest_fastas_from_list(args.local_mirror, organism_list)
    elif args.from_list:
        check_dirs(args.local_mirror)
        get_assembly_summary(args.no_wget, args.local_mirror)
        organism_list = args.from_list
        get_latest_fastas_from_list(args.local_mirror, organism_list)
    elif args.genus:
        check_dirs(args.local_mirror)
        get_assembly_summary(args.no_wget, args.local_mirror)
        organism_list = get_species_list_from_genus(args.genus)
        get_latest_fastas_from_list(args.local_mirror, organism_list)
    elif args.reorganize:
        check_dirs(args.local_mirror)
        get_assembly_summary(args.no_wget, args.local_mirror)
        mk_dir_structure(args.local_mirror)
    else:
        check_dirs(args.local_mirror)
        get_assembly_summary(args.no_wget, args.local_mirror)
        fasta_list = ftp_paths_from_assembly_summary(args.local_mirror)
        get_latest_fastas_from_assembly_summary(args.local_mirror, fasta_list)

    #monitor_changes(args.local_mirror, organism_list)

Main()
