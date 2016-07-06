#!/usr/bin/env python3

import os, subprocess, argparse, rename_fastas, gzip, time, shutil
from ftplib import FTP
from ftplib import error_temp

def get_assembly_summary(wget):

    """Get current version of assembly_summary.txt"""

    assembly_summary = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
    if wget:
        print("Fetching current version of assembly_summary.txt.")
        try:
            os.remove('assembly_summary.txt')
            subprocess.call(['wget', '-nv', assembly_summary])
        except OSError:
            subprocess.call(['wget', '-nv', assembly_summary])
    else:
        print("assembly_summary.txt will not be downloaded.")

def ftp_login(directory="genomes/genbank/bacteria"):

    """Login to ftp.ncbi.nlm.nih.gov"""

    ftp_site = 'ftp.ncbi.nlm.nih.gov'
    ftp = FTP(ftp_site)
    ftp.login()
    ftp.cwd(directory)

    return ftp

def get_complete_organism_list():

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    print("Getting list of all bacteria directories to sync with.")
    print("Estimated wait time:  1 minute")
    ftp = ftp_login()
    organism_list = ftp.nlst()

    return organism_list

def get_organism_list_from_file(source):

    """Read organism list from file"""

    print("Getting list of directories to sync with.")
    organism_list = [] 
    with open(source, 'r') as f:
        for line in f:
            line = line.strip()
            organism_list.append(line)
    return organism_list

def check_dirs(local_mirror):

    """Make directories to store fastas if they don't already exist."""

    renamed_fastas_dir =  "_".join([local_mirror, "renamed"])
    if not os.path.isdir(local_mirror):
        os.mkdir(local_mirror)

    if not os.path.isdir(renamed_fastas_dir):
        os.mkdir(renamed_fastas_dir)

def write_filter_list(local_mirror, organism, ftp):
    
    """Create accepted files list for each organism for rsync filter."""

    filter_files_dir = os.path.join(local_mirror, "filter_files")

    if not os.path.isdir(filter_files_dir):
        os.mkdir(filter_files_dir)
        print("Created {} to store filter lists.".format(filter_files_dir))

    print("Creating accepted files list for organism.")

    latest = os.path.join(organism, 'latest_assembly_versions')
    filter_file = os.path.join(filter_files_dir, organism+".txt")

    if os.path.isfile(filter_file):
        os.remove(filter_file)

    with open(filter_file, "a") as f:
        try:
            for full_path in ftp.nlst(latest):
                accession = full_path.split("/")[-1]
                fasta = "_".join([accession, "genomic.fna.gz"])
                fasta = "/".join([accession, fasta])
                f.write(fasta+"\n")
        except error_temp: 
            pass

    return filter_file

def changes_log(local_mirror):
    changes_log = os.path.join(local_mirror, "changes_log.txt")
    with open(changes_log, "a") as log:
        log.write("start time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")

def get_fastas(local_mirror, organism_list):

    renamed_fasta_dir = "_".join([local_mirror,"renamed"])
    ftp = ftp_login()

    for organism in organism_list:
        filter_file = write_filter_list(local_mirror, organism, ftp)
        print(filter_file)
        rsync_log = os.path.join(organism, "rsync_log.txt")
        latest = os.path.join(organism, 'latest_assembly_versions')
        renamed_organism_dir = os.path.join(renamed_fasta_dir, organism)
        try:
            organism_dir = os.path.join(local_mirror, organism)
            subprocess.call(['rsync',
                            '--chmod=ugo=rwX', # Change permissions so files can be copied/renamed
                            '--info=progress2',
                            '--copy-links', # Follow symbolic links
                            '--recursive',
                            '--times',
                            '--itemize-changes',
                            '--prune-empty-dirs',
                            '--files-from='+filter_file,
                            #'--log-file='+rsync_log,
                            'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/'+latest,
                            organism_dir])
        except error_temp: # Skips organisms that don't have latest_assembly_versions directory
            continue

        move_and_unzip(local_mirror, organism)
        organism_index = organism_list.index(organism)+1
        print("\n{} out of {} organisms downloaded.\n".format(organism_index, len(organism_list)))

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

    renamed_fasta_dir = "_".join([local_mirror,"renamed"])
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

def monitor_changes(local_mirror, organism_list):
    renamed_fastas_dir = local_mirror + "_renamed"
    dirs = get_organism_list(organism_list)
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
    parser = argparse.ArgumentParser(description = "Sync with NCBI's database and organize them in a sane way.")
    parser.add_argument('local_mirror', help = 'Directory to save fastas', type=str)
    parser.add_argument('-W', '--no-wget', help = "Skip downloading current assembly_summary.txt", action='store_false')
    parser.add_argument('-f', '--from_file', help = 'Full path to file containing directories to sync with.')
    parser.add_argument('-l', '--from_list', help = 'List of organisms to download', nargs="+")
    args = parser.parse_args()

    if args.from_file:
        organism_list = get_organism_list_from_file(args.from_file)
    elif args.from_list:
        organism_list = args.from_list
    else:
        organism_list = get_complete_organism_list()

    check_dirs(args.local_mirror)
    get_assembly_summary(args.no_wget)
    get_fastas(args.local_mirror, organism_list)
    #monitor_changes(args.local_mirror, organism_list)

Main()
