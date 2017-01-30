#!/usr/bin/env python

import os
import argparse
import prun
import curate
import sync.sync_genbank as sync
from re import sub
from time import sleep, strftime
from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list

def get_resources(genbank_mirror):

    """
    Get assembly_summary.txt for bacteria and taxonomy dump file.
    Parse and load into Pandas DataFrames.
    """

    assembly_summary = sync.get_assembly_summary(genbank_mirror)
    names = sync.get_species_names(genbank_mirror)

    return assembly_summary, names

def parallel(genbank_mirror):
    path_vars = curate.instantiate_path_vars(genbank_mirror)
    get_latest_job_id = prun.get_latest(genbank_mirror, path_vars)
    prun.update_genomes(genbank_mirror, path_vars, get_latest_job_id)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--slurm", help = 'Submit jobs in parallel via SLURM arrays.', action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    curate.clean_up(genbank_mirror)
    path_vars = curate.instantiate_path_vars(genbank_mirror)

    if args.slurm:
        parallel(genbank_mirror)
    else:
        assembly_summary, names = get_resources(genbank_mirror)
        species_taxids = curate.get_species_taxids(assembly_summary)
        species_list = curate.species_list_from_taxdmp(species_taxids, names)
        curate.check_species_dirs(genbank_mirror, species_list)
        sync.sync_latest_genomes(genbank_mirror, assembly_summary, names)

main()
