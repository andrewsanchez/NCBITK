#!/usr/bin/env python3

import os, subprocess, argparse, rename_fastas, gzip, time, shutil
from ftplib import FTP
from ftplib import error_temp

# Get current version of assembly_summary.txt
def get_assembly_summary(wget):
    assembly_summary = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
    if wget:
        try:
            os.remove('assembly_summary.txt')
            subprocess.call(['wget', assembly_summary])
        except OSError:
            subprocess.call(['wget', assembly_summary])
    else:
        pass

def get_organism_list(input_file):
    if input_file:
        # Read from file
        if 'txt' in input_file:
            dirs = [] 
            with open(input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    dirs.append(line)
            return dirs
        # Read from list given as an arg to -l
        else:
            dirs = input_file.split(',')
            return dirs

    # Connect with ftp.ncbi and get list of directories
    else:
        ftp_site = 'ftp.ncbi.nlm.nih.gov'
        ftp = FTP(ftp_site)
        ftp.login()
        ftp.cwd('genomes/genbank/bacteria')
        dirs = ftp.nlst()
        return dirs

# Make directories if they don't already exist
def check_dirs(local_mirror):
    all_fastas =  local_mirror+'_renamed/'
    if not os.path.isdir(local_mirror):
        os.mkdir(local_mirror)

    if not os.path.isdir(all_fastas):
        os.mkdir(all_fastas)

def cp_files(source, destination):
    # Copy new files to destination
    for root, dirs, files in os.walk(source):
        for f in files:
            try:
                source = os.path.join(root, f)
                copy_to = os.path.join(destination, f)
                shutil.copyfile(source, copy_to)
            # Skip files that already exist in destination
            except shutil.SameFileError:
                continue

# Unzip fastas
def gunzip(target_dir):
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

def get_fastas(local_mirror, organism_list):
    dirs = get_organism_list(organism_list)
    rsync_log = os.path.join(local_mirror, "rsync_log.txt")
    changes_log = os.path.join(local_mirror, "changes_log.txt")
    renamed = local_mirror+'_renamed'

    with open(changes_log, "a") as log:
        log.write("start time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")

    for organism in dirs:
        ftp_site = 'ftp.ncbi.nlm.nih.gov'
        ftp = FTP(ftp_site)
        ftp.login() # Reestablish ftp connection for each organism
        ftp.cwd('genomes/genbank/bacteria') # Connection may timeout otherwise
        latest = os.path.join(organism, 'latest_assembly_versions')
        path_to_renamed = os.path.join(renamed, organism)
        try:
            for parent in ftp.nlst(latest):
                accession = parent.split("/")[-1]
                fasta = accession+"_genomic.fna.gz" # The file to be downloaded
                organism_dir = os.path.join(local_mirror, organism)
                subprocess.call(['rsync',
                                '--chmod=ugo=rwX', # Change permissions so files can be copied/renamed
                                '--info=progress2',
                                '--copy-links', # Follow symbolic links
                                '--recursive',
                                '--times',
                                '--itemize-changes',
                                '--prune-empty-dirs',
                                '-f=+ '+accession,
                                '-f=+ '+fasta,
                                '--exclude=*',
                                'ftp.ncbi.nlm.nih.gov::genomes/genbank/bacteria/'+parent,
                                organism_dir])
        # Skips organisms that don't have the /latest_assembly_versions/ directory
        except error_temp:
            continue

        # Copy files to different directory for renaming/decompression
        if os.path.isdir(path_to_renamed):
            source = os.path.join(local_mirror, organism)
            cp_files(source, path_to_renamed)
            gunzip(path_to_renamed)

        else:
            os.mkdir(path_to_renamed)
            source = os.path.join(local_mirror, organism)
            cp_files(source, path_to_renamed)
            gunzip(path_to_renamed)

        rename_fastas.rename(path_to_renamed)

def monitor_changes(local_mirror, organism_list):
    all_fastas = local_mirror + "_renamed"
    dirs = get_organism_list(organism_list)
    countdirs = str(len(dirs))
    local_dir_count = str(len(os.listdir(all_fastas)))
    missingdirs = int(countdirs) - int(local_dir_count)
    changes_log = os.path.join(local_mirror, "changes_log.txt")
    with open(changes_log, "a") as log:
        log.write("Finish time:  " + time.strftime("%m/%d/%y - %H:%M"))
        log.write("\n")
        log.write('Dirs at ftp.ncbi:  ' + countdirs)
        log.write("\n")
        log.write('Dirs in {}:  {}\n'.format(all_fastas, local_dir_count))
        log.write('Missing dirs = ' + str(missingdirs))
        log.write("\n")
        for i in dirs:
            if i not in os.listdir(local_mirror):
                log.write(i+"\n")

def Main():
    parser = argparse.ArgumentParser(description = "Sync with NCBI's database, give the files useful names,"\
            "and organize them in a sane way.")
    parser.add_argument('local_mirror', help = 'directory to save fastas to', type=str)
    parser.add_argument('-w', '--wget', help = "Fetch assembly_summary.txt", action='store_true', default=False)
    parser.add_argument('-i', '--from_file', help = 'Input file containing directories to sync with.')
    parser.add_argument('-l', '--from_list', help = 'Comma separated list of directories to be downloaded')
    args = parser.parse_args()

    if args.from_list:
        get_assembly_summary(args.wget)
        check_dirs(args.local_mirror)
        get_fastas(args.local_mirror, args.from_list)
        monitor_changes(args.local_mirror, args.from_list)
            
    elif args.from_file:
        get_assembly_summary(args.wget)
        check_dirs(args.local_mirror)
        get_fastas(args.local_mirror, args.from_file)
        monitor_changes(args.local_mirror, args.from_file)

    else:
        get_assembly_summary(args.wget)
        check_dirs(args.local_mirror)
        get_fastas(args.local_mirror, args.from_file)
        monitor_changes(args.local_mirror, args.from_file)

Main()
