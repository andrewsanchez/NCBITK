#!/usr/bin/env python

import os
import argparse
import prun
import curate
from re import sub
from time import sleep, strftime
from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list 

def parallel(genbank_mirror):
    paths = curate.instantiate_path_vars(genbank_mirror)
    get_latest_job_id = prun.get_latest(genbank_mirror, paths)
    prun.update_genomes(genbank_mirror, get_latest_job_id)
    #dependency_id = update_genomes(genbank_mirror, get_latest_job_id)
#   cmd = 'python /common/contrib/tools/NCBITK/sync/rename.py {}'.format(genbank_mirror)
#   salloc(cmd, dependency_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--parallel", help = 'Submit jobs in parallel via SLURM arrays.', action="store_true")
    parser.add_argument("-s", "--sync", action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    curate.clean_up(genbank_mirror)
    path_vars = curate.instantiate_path_vars(genbank_mirror)

    if args.parallel:
        parallel(genbank_mirror)
    elif args.sync:
        update_genomes(genbank_mirror) # submit this via sbatch in the cron job - make it depend on latest_assembly_versions_script

main()
