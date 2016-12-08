#!/usr/bin/env python

import os
import argparse
from re import sub
from time import sleep
from time import strftime
from slurm.generate_arrays import *
from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list 
from sync.sync_genbank import unzip_genbank_mirror, get_assembly_summary
from sync.rename import *

def get_latest(genbank_mirror):

    complete_species_list = ftp_complete_species_list()[:50]
    latest_assembly_versions_array = gen_latest_assembly_versions_array(genbank_mirror, complete_species_list)
    latest_assembly_versions_script = gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array)
    get_latest_job_id = submit_sbatch(latest_assembly_versions_script)
    last_task_number = len(list(open(latest_assembly_versions_array )))
    get_latest_job_id = '{}_{}'.format(get_latest_job_id, last_task_number)

    return get_latest_job_id

def update_genomes(genbank_mirror, get_latest_job_id):

    sync_array_script = gen_sync_array_script(genbank_mirror, get_latest_job_id)
    sleep(30) # make sure the following job doesn't get submitted too soon.
    sync_array_job_id = submit_sbatch(sync_array_script)
    grab_genomes_script = gen_grab_genomes_script(genbank_mirror, sync_array_job_id)
    job_id = submit_sbatch(grab_genomes_script)

def parallel(genbank_mirror):
    get_latest_job_id = get_latest(genbank_mirror)
    update_genomes(genbank_mirror, get_latest_job_id)

def dir_vars(genbank_mirror):
    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")
    dirs = info_dir, slurm, out
    return dirs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--parallel", action="store_true")
    parser.add_argument("-s", "--sync", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror
    dirs = dir_vars(genbank_mirror)

    clean_up(genbank_mirror)
    assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
    if args.parallel:
        parallel(genbank_mirror)
    if args.sync:
        update_genomes(genbank_mirror) # submit this via sbatch in the cron job - make it depend on latest_assembly_versions_script
    # might have to submit these via sbatch as well...
    # or just submit them a few hours after the above jobs
    unzip_genbank_mirror(genbank_mirror)
    rename(genbank_mirror, assembly_summary)

main()
