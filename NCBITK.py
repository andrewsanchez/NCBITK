#!/usr/bin/env python

import os
import argparse
from time import strftime
from subprocess import Popen
from re import sub
from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list 
from slurm.generate_arrays import *

def get_latest(genbank_mirror):
    latest_assembly_versions_array = gen_latest_assembly_versions_array(genbank_mirror, complete_species_list)
    get_latest_script = gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array)
    job_id, errs = Popen("sbatch {}".format(get_latest_script), shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    return job_id

def update_genomes(genbank_mirror, job_id):
    sync_array = gen_sync_array(genbank_mirror) # (needs to depend on completion of latest_assembly_versions_script)
    sync_script = gen_sync_script(genbank_mirror, sync_array)
    cmd = "sbatch -d=after:{} {}".format(job_id, sync_script)
    print(cmd)
    Popen(cmd, shell="True")

def parallel(genbank_mirror):
    clean_up(genbank_mirror)
    complete_species_list = ftp_complete_species_list()
    job_id = get_latest(genbank_mirror)
    update_genomes(genbank_mirror)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--parallel", action="store_true")
    parser.add_argument("-s", "--sync", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror
    if args.parallel:
        get_latest(genbank_mirror)
        update_genomes(genbank_mirror)
    if args.sync:
        update_genomes(genbank_mirror) # submit this via sbatch in the cron job - make it depend on latest_assembly_versions_script

main()
