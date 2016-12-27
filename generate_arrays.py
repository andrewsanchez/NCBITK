#!/usr/bin/env python

import os
import argparse
import subprocess
import pandas as pd
import curate
from re import sub
from time import strftime, sleep

ymd = strftime("%y.%m.%d")

def get_local_and_latest(genbank_mirror, latest_assembly_versions):

    species = list(set(latest_assembly_versions['species']))
    species_dirs = []
    for name in species:
        species_dirs.append(os.path.join(genbank_mirror, name))

    print(species_dirs)
    local_genome_ids = []
    for species_dir in species_dirs:
        for local_genome in os.listdir(species_dir):
            local_genome_id = "_".join(local_genome.split("_")[:2])
            local_genome_ids.append(local_genome_id)

    #  local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
    latest_genome_ids = latest_assembly_versions.index[latest_assembly_versions['species'] == species].tolist()
    latest_genome_paths = latest_assembly_versions.dir[latest_assembly_versions['species'] == species].tolist()
    ids_and_paths = zip(latest_genome_ids, latest_genome_paths)
    print(local_genome_ids)
    print(latest_genome_ids)
    print(latest_genome_paths)

    return local_genome_ids, latest_genome_ids, ids_and_paths

def instantiate_path_vars(genbank_mirror):
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")
    return info_dir, slurm, out

def instantiate_df(path, cols):
    df = pd.read_csv(path, index_col=0, header=None)
    df.columns = cols
    return df

def gen_latest_assembly_versions_array(genbank_mirror, complete_species_list):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    print('Generating {}'.format(latest_assembly_versions_array))
    groups = [complete_species_list[n:n+1000] for n in range(0, len(complete_species_list), 1000)]
    with open(latest_assembly_versions_array, "a") as f:
        for group in groups:
            group = ' '.join(group)
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/get_latest_assembly_versions.py {} {}\n".format(genbank_mirror, group))
    array_len = len(list(open(latest_assembly_versions_array)))

    return latest_assembly_versions_array, array_len

def gen_sbatch_script(genbank_mirror, array, job_name, time):
    None

def gen_sbatch_array_script(genbank_mirror, array, job_name, time):
    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    out = os.path.join(out, "{}_%a.out".format(job_name))
    sbatch_script = os.path.join(slurm, "{}.sbatch".format(job_name))
    print('Generating {}'.format(sbatch_script))
    array_len = len(list(open(sbatch_script)))
    with open(sbatch_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time={}\n".format(time))
        f.write("#SBATCH --job-name={}\n".format(job_name))
        f.write("#SBATCH --output={}\n".format(out)) # can I remove this line to avoid getting out files?
        f.write("#SBATCH --array=1-{}%2\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(array))
        f.write("srun $cmd")

    return sbatch_script 

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    out = os.path.join(out, "get_latest_%a.out")
    latest_assembly_versions_script = os.path.join(slurm, "get_latest_assembly_versions.sbatch")
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
    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    sync_array_script = os.path.join(slurm, 'sync_array_script.sbatch')
    print('Generating {}'.format(sync_array_script))

    # consolidate into function
    with open(sync_array_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=01:00\n")
        f.write("#SBATCH --job-name=gen_sync_array\n")
        f.write("#SBATCH --output={}/gen_sync_array.out\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(get_latest_job_id))
        f.write('cmd="python /common/contrib/tools/NCBITK/slurm/generate_arrays.py {}"\n'.format(genbank_mirror))
        f.write("srun $cmd")
    
    return sync_array_script 

def gen_grab_genomes_script(genbank_mirror, sync_array_job_id):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    sync_array = os.path.join(slurm, "sync_array.txt")
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')
    out = os.path.join(out, 'grab_genomes%a.out')
    print('Generating {}'.format(grab_genomes_script))

    while not os.path.isfile(sync_array):
        sleep(10)
    else:
        sync_array_len = len(list(open(sync_array)))

    with open(grab_genomes_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=04:30\n")
        f.write("#SBATCH --job-name=grab_genomes\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(sync_array_job_id))
        f.write("#SBATCH --array=1-{}%2\n".format(sync_array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(sync_array))
        f.write("srun $cmd")
    
    return grab_genomes_script, sync_array_len

def write_grab_genomes_array(genbank_mirror):

    latest_assembly_versions = curate.read_latest_assembly_versions(genbank_mirror)
    grab_genomes_array = os.path.join(genbank_mirror, ".info", "slurm", "grab_genomes_array.txt")
    print('Generating {}'.format(grab_genomes_array))
    new_genomes = curate.get_new_genomes(genbank_mirror, latest_assembly_versions)
    args = []
    for name in new_genomes:
        species = latest_assembly_versions.loc[name, 'species']
        path = latest_assembly_versions.loc[name, 'dir']
        args.append(','.join([species, name, path]))

    groups = [args[n:n+2000] for n in range(0, len(args), 2000)] 
    with open(grab_genomes_array, "a") as f:
        for group in groups:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n".format(genbank_mirror, ' '.join(group)))

    return grab_genomes_array

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("-g", '--grab_genomes', action='store_true')
    args = parser.parse_args()

    if args.grab_genomes:
        write_grab_genomes_array(args.genbank_mirror)

if __name__ == "__main__":
    main()
