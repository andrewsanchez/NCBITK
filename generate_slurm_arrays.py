#!/usr/bin/env python

import argparse, os, sys
from ftp_funs.ftp_funs import *
from subprocess import Popen

def gen_latest_assembly_versions_array(genbank_mirror, complete_species_list):

    print("Generating slurm array.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    with open(latest_assembly_versions_array, "a") as f:
        for species in complete_species_list:
            f.write("python /common/contrib/tools/NCBITK/get_latest_assembly_versions.py {} {}\n".format(genbank_mirror, species))

    return latest_assembly_versions_array

def gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array):

    print("Generating slurm script.")
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    array_len = len(list(open(latest_assembly_versions_array)))
    with open(slurm_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=03:00\n")
        f.write("#SBATCH --job-name=get_latest\n")
        f.write("#SBATCH --workdir={}\n".format(genbank_mirror))
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

    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sh")
    for f in [latest_assembly_versions_array, slurm_script]:
        if os.path.isfile(f):
            os.remove(f)

def main():

    parser = argparse.ArgumentParser(description = "Generate slurm array commands")
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-l", "--latest", help = "Generate get latest assembly versions array", action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    clean_up(genbank_mirror)
    complete_species_list = ftp_complete_species_list()
    if args.latest:
        latest_assembly_versions_array = gen_latest_assembly_versions_array(genbank_mirror, complete_species_list)
        slurm_script = gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array)
        #get_latest_assembly_versions_cmd = "bash {}".format(get_latest_assembly_versions)
        #print("Submitting slurm array.")
        #Popen(slurm_script, shell="True")

main()
