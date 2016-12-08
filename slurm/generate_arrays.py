#!/usr/bin/env python

import os
import argparse
import subprocess
import pandas as pd
from re import sub
from time import strftime, sleep

ymd = strftime("%y.%m.%d")

def submit_sbatch(slurm_script):

    job_id, errs = subprocess.Popen("sbatch {}".format(slurm_script), # submit array
            shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    print(job_id)
    job_id = job_id.split(" ")[-1].strip()

    return job_id


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

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    print('Generating {}'.format(latest_assembly_versions_array))
    groups = [complete_species_list[n:n+10] for n in range(0, len(complete_species_list), 10)]
    with open(latest_assembly_versions_array, "a") as f:
        for group in groups:
            group = ' '.join(group)
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/get_latest_assembly_versions.py {} {}\n".format(genbank_mirror, group))

    return latest_assembly_versions_array

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", "get_latest_%a.out")
    latest_assembly_versions_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    print('Generating {}'.format(latest_assembly_versions_script))
    array_len = len(list(open(latest_assembly_versions_array)))
    with open(latest_assembly_versions_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=get_latest\n")
        f.write("#SBATCH --output={}\n".format(out)) # can I remove this line to avoid getting out files?
        f.write("#SBATCH --array=1-{}%2\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(latest_assembly_versions_array))
        f.write("srun $cmd")

    return latest_assembly_versions_script 

def gen_sync_array_script(genbank_mirror, get_latest_job_id):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", 'gen_sync_array.out')
    sync_array_script = os.path.join(slurm, 'sync_array_script.sh')
    print('Generating {}'.format(sync_array_script))

    if os.path.isfile(sync_array_script):
        os.remove(sync_array_script)

    with open(sync_array_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=gen_sync_array\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(get_latest_job_id))
        f.write('cmd="python /common/contrib/tools/NCBITK/slurm/generate_arrays.py {}"\n'.format(genbank_mirror))
        f.write("srun $cmd")
    
    return sync_array_script 

def gen_grab_genomes_script(genbank_mirror, sync_array_job_id):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", 'grab_genomes.out')
    sync_array = os.path.join(slurm, "sync_array.txt")
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sh')
    print('Generating {}'.format(grab_genomes_script))

    if os.path.isfile(grab_genomes_script):
        os.remove(grab_genomes_script)

    while not os.path.isfile(sync_array):
        sleep(1)
    else:
        sync_array_len = len(list(open(sync_array)))

    with open(grab_genomes_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=grab_genomes\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(sync_array_job_id))
        f.write("#SBATCH --array=1-{}%2\n".format(sync_array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(sync_array))
        f.write("srun $cmd")
    
    return grab_genomes_script

def write_grab_genome_commands(genbank_mirror):

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")

    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    print('Generating {}'.format(sync_array))
    if os.path.isfile(sync_array):
        os.remove(sync_array)
    info = [i.strip() for i in list(open(latest_assembly_versions))]
    groups = [info[n:n+25] for n in range(0, len(info), 25)] 
    with open(sync_array, "a") as f:
        for group in groups:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n".format(genbank_mirror, ' '.join(group)))

    return sync_array


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    args = parser.parse_args()

    write_grab_genome_commands(args.genbank_mirror)

if __name__ == "__main__":
    main()
