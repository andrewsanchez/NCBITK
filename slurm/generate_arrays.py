#!/usr/bin/env python

import os
import argparse
import pandas as pd
from re import sub
from time import strftime
from ftplib import FTP, error_temp

ymd = strftime("%y.%m.%d")

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

def instantiate_df(path, cols):
    df = pd.read_csv(path, index_col=0, header=None)
    df.columns = cols
    return df

def clean_up(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")
    for d in [genbank_mirror, info_dir, slurm, out]:
        if not os.path.isdir(d):
            os.mkdir(d)

    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions.csv")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    for f in [latest_assembly_versions, latest_assembly_versions_array, slurm_script]:
        if os.path.isfile(f):
            os.remove(f)

def check_species_dirs(genbank_mirror, species):

    species_dir = os.path.join(genbank_mirror, species)
    if not os.path.isdir(species_dir):
        os.mkdir(species_dir)

    return species_dir

def check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids):

    for genome_id in local_genome_ids:
        if genome_id not in latest_genome_ids:
            fasta = glob("{}*".format(genome_id))
            os.remove(os.path.join(genbank_mirror, species, fasta[0]))

def gen_latest_assembly_versions_array(genbank_mirror, complete_species_list):

    print("Generating slurm array.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    groups = [complete_species_list[n:n+44] for n in range(0, len(complete_species_list), 44)]
    with open(latest_assembly_versions_array, "a") as f:
        for group in groups:
            group = ' '.join(group).strip()
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/get_latest_assembly_versions.py {} {}\n".format(genbank_mirror, group))

    return latest_assembly_versions_array

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    print("Generating slurm script.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", "get_latest_%a.out")
    get_latest_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    array_len = len(list(open(latest_assembly_versions_array)))
    with open(get_latest_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=00:30\n")
        f.write("#SBATCH --job-name=get_latest\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --array=1-{}\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(latest_assembly_versions_array))
        f.write("srun $cmd")

    return get_latest_script 

def gen_sync_array_new(genbank_mirror):

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    if os.path.isfile(sync_array):
        os.remove(sync_array)
    info = [i.strip() for i in list(open(latest_assembly_versions))]
    groups = [info[n:n+30] for n in range(0, len(info), 30)]
    print(len(groups))
    with open(sync_array, "a") as f:
        for group in groups:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n".format(genbank_mirror, ' '.join(group)))

    return sync_array

def gen_sync_array(genbank_mirror):

    latest_assembly_versions = instantiate_df(os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv"), ["id", "dir"])
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    if os.path.isfile(sync_array):
        os.remove(sync_array)
    with open(sync_array, "a") as f:
        species_directories = list(set(latest_assembly_versions.index))
        info = []
        for species in species_directories:
            species_dir = check_species_dirs(genbank_mirror, species)
            print("Assessing {}".format(species_dir))
            local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
            latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
            latest_genome_paths = [sub("[\[\]]", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
            check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids)
            for id in latest_genome_ids:
                info.append((species, latest_genome_ids, latest_genome_paths))

            #ids_and_paths = list(zip(latest_genome_ids, latest_genome_paths))
        groups = [info[n:n+30] for n in range(0, len(info), 30)]
        for group in groups:
            print(group)
            grab_genome_args = []
            for tup in group:
                print(tup)
                species = tup[0]
                genome_id = tup[1]
                genome_path = tup[2]
                print(species, genome_id, genome_path)
                if genome_id not in local_genome_ids:
                    grab_genome_args.append(','.join([species,genome_id,genome_path]))

            print(grab_genome_args)
            for arg in grab_genome_args:
                print("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n"\
                        .format(genbank_mirror, ' '.join(arg)))
                f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n"\
                        .format(genbank_mirror, ' '.join(arg)))

           #for info in ids_and_paths:
           #    genome_id = info[0]
           #    genome_path = info[1]
           #    if genome_id not in local_genome_ids:
           #        f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {} {} {}\n"\
           #                .format(genbank_mirror, species, genome_id, genome_path))

    #return sync_array

def gen_sync_script(genbank_mirror, sync_array):

    print("Generating slurm script.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", "grab_genome_%a.out")
    sync_script = os.path.join(slurm, "sync_genbank.sh")
    if os.path.isfile(sync_script):
        os.remove(sync_script)
    array_len = len(list(open(sync_array)))
    with open(sync_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=00:30\n")
        f.write("#SBATCH --job-name=grab_genome\n")
        f.write("#SBATCH --output={}".format(out))
        f.write("#SBATCH --array=1-{}%5\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(sync_array))
        f.write("srun $cmd")

    return sync_script 

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    args = parser.parse_args()

    gen_sync_array_new(args.genbank_mirror)

if __name__ == "__main__":
    main()
