#!/usr/bin/env python

import os
import argparse
from NCBITK import config
from NCBITK import sync
from NCBITK import curate
from NCBITK import get_resources
from re import sub
from time import sleep, strftime
# from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list

def parallel(genbank_mirror):
    path_vars = config.instantiate_path_vars(genbank_mirror)
    get_latest_job_id = prun.get_latest(genbank_mirror, path_vars)
    prun.update_genomes(genbank_mirror, path_vars, get_latest_job_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--slurm", help = 'Submit jobs in parallel via SLURM arrays.', action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    path_vars = config.instantiate_path_vars(genbank_mirror)
    curate.clean_up(genbank_mirror, path_vars)
    assembly_summary = get_resources.get_resources(genbank_mirror)
    local_genomes, new_genomes = assess_genbank_mirror(genbank_mirror, assembly_summary)
    # wrap the following into curate.update_genbank_mirror(genbank_mirror, assembly_summary, local_genomes, new_genomes)
    curate.create_species_dirs(genbank_mirror, assembly_summary)
    curate.remove_old_genomes(genbank_mirror, assembly_summary, local_genomes)
    sync.sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes)

main()
