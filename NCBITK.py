#!/usr/bin/env python

import os
import argparse
from time import strftime
import subprocess
from re import sub
from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list 
from slurm.generate_arrays import *
from time import sleep

def get_latest(genbank_mirror):

    complete_species_list = ftp_complete_species_list()[:100]
    latest_assembly_versions_array = gen_latest_assembly_versions_array(genbank_mirror, complete_species_list)
    latest_assembly_versions_script = gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array)

    job_id, errs = subprocess.Popen("sbatch {}".format(latest_assembly_versions_script), # submit array
            shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    print(job_id)
    job_id = job_id.split(" ")[-1].strip()

    return job_id

def update_genomes(genbank_mirror, job_id):

    # array_len used to identify the task number of the last job in the array for --dependency
    array_len = len(list(open(os.path.join(genbank_mirror, '.info', 'slurm', 'latest_assembly_versions_array.txt'))))

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out", 'gen_sync_array.out')
    sync_array_script = os.path.join(slurm, 'sync_array_script.sh')
    if os.path.isfile(sync_array_script):
        os.remove(sync_array_script)

    def gen_sync_array_script(job_id, array_len):
        with open(sync_array_script, 'a') as f:
            f.write("#!/bin/sh\n")
            f.write("#SBATCH --time=05:00\n")
            f.write("#SBATCH --job-name=gen_sync_array\n")
            f.write("#SBATCH --output={}\n".format(out))
            f.write("#SBATCH --dependency={}_{}\n".format(job_id, array_len))
            f.write('cmd="python /common/contrib/tools/NCBITK/slurm/generate_arrays.py {}"\n'.format(genbank_mirror))
            f.write("srun $cmd")

    gen_sync_array_script(job_id, array_len)
    sleep(40) # make sure the following job doesn't get submitted too soon.
    job_id, errs = subprocess.Popen("sbatch {}".format(sync_array_script), # submit array
            shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    job_id = job_id.split(" ")[-1].strip()
    print(job_id)
    # generate the grab genomes array
#   subprocess.Popen("salloc --dependency=after:{}_{} python /common/contrib/tools/NCBITK/slurm/generate_arrays.py -t 05:00 {} {}".format(job_id, array_len, genbank_mirror), shell="True")

    # run this from a script and submit via sbatch so that that job_id is accessible.
  # job_id, errs = subprocess.Popen("salloc --dependency=after:{}_{} -t 05:00 python /common/contrib/tools/NCBITK/slurm/generate_arrays.py \
  #         {}".format(job_id, array_len, genbank_mirror), shell="True", stdout=subprocess.PIPE, universal_newlines=True).communicate()

  # job_id = job_id.split(" ")[-1].strip()
    #sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    #sync_array = gen_sync_array(genbank_mirror, job_id) # (needs to depend on completion of latest_assembly_versions_script)
   #sync_script = gen_sync_script(genbank_mirror, sync_array)
   #subprocess.Popen(cmd, shell="True")

def parallel(genbank_mirror):
    job_id = get_latest(genbank_mirror)
    update_genomes(genbank_mirror, job_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--parallel", action="store_true")
    parser.add_argument("-s", "--sync", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror
    clean_up(genbank_mirror)
    if args.parallel:
        parallel(genbank_mirror)
    if args.sync:
        update_genomes(genbank_mirror) # submit this via sbatch in the cron job - make it depend on latest_assembly_versions_script

main()
