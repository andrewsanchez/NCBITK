#!/usr/bin/env python

import os
from time import strftime
from ftp_functions.ftp_functions import ftp_complete_species_list

ymd = strftime("%y.%m.%d")

def gen_latest_assembly_versions_array(genbank_mirror, complete_species_list):

    print("Generating slurm array.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array_{}.txt".format(ymd))
    with open(latest_assembly_versions_array, "a") as f:
        for species in complete_species_list:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/get_latest_assembly_versions.py {} {}\n".format(genbank_mirror, species))

    return latest_assembly_versions_array

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    print("Generating slurm script.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    array_len = len(list(open(latest_assembly_versions_array)))
    with open(slurm_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=00:30\n")
        f.write("#SBATCH --job-name=get_latest\n")
        f.write("#SBATCH --array=1-{}%3\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(latest_assembly_versions_array))
        f.write("srun $cmd")

    return slurm_script 

def clean_up(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    for d in [genbank_mirror, info_dir, slurm]:
        if not os.path.isdir(d):
            os.mkdir(d)

    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions_{}.txt".format(ymd))
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array_{}.txt".format(ymd))
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    for f in [latest_assembly_versions, latest_assembly_versions_array, slurm_script]:
        if os.path.isfile(f):
            os.remove(f)

def generate_sync_array(latest_assembly_versions):

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions_{}.txt".format(ymd))
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None).head()
    latest_assembly_versions.columns = ["dir", "id"]
    species_directories = list(set(latest_assembly_versions.index))
    for species in species_directories:
        None
