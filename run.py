#!/usr/bin/env python

import os
import argparse
import logging
from NCBITK import config
from NCBITK import sync
from NCBITK import curate
from NCBITK import get_resources
from re import sub
# from ftp_functions.ftp_functions import ftp_login, ftp_complete_species_list

def parallel(genbank_mirror):
    path_vars = config.instantiate_path_vars(genbank_mirror)
    get_latest_job_id = prun.get_latest(genbank_mirror, path_vars)
    prun.update_genomes(genbank_mirror, path_vars, get_latest_job_id)

def update_genbank_mirror(genbank_mirror, species_list="all", fetch_new=True):

    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror, logger, fetch_new)
    curate.create_species_dirs(genbank_mirror, assembly_summary, logger, species_list)
    local_genomes, new_genomes, sketch_files, missing_sketch_files = curate.assess_genbank_mirror(genbank_mirror, assembly_summary, species_list)
    curate.remove_old_genomes(genbank_mirror, assembly_summary, local_genomes, logger)
    sync.sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes, logger)

    logger.info('{} genomes in assembly_summary.txt'.format(len(assembly_summary)))
    logger.info("{} genomes present in local collection.".format(len(local_genomes)))
    logger.info("{} genomes not in local collection.".format(len(new_genomes)))
    logger.info('{} sketch files present in local collection.'.format(len(sketch_files)))
    logger.info('{} sketch files not in local collection.'.format(len(missing_sketch_files)))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-s", "--species", help = 'List of species', nargs='+', default='all')
    parser.add_argument("-p", "--slurm", help = 'Submit jobs in parallel via SLURM arrays.', action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    update_genbank_mirror(genbank_mirror, args.species)

if __name__ == "__main__":
    main()
