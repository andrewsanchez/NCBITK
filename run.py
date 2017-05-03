#!/usr/bin/env python

import os
import argparse
import logging

from NCBITK import config
from NCBITK import sync
from NCBITK import curate
from NCBITK import get_resources

def setup(genbank_mirror, species_list, update_assembly_summary):

    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror, logger, update_assembly_summary)
    species = curate.get_species_list(assembly_summary, species_list)

    return path_vars, assembly_summary, species

def assess_genbank(genbank_mirror, assembly_summary, species_list, logger):

    genbank_status = curate.assess_genbank_mirror(genbank_mirror, assembly_summary, species_list, logger)
    local_genomes, local_genome_paths, new_genomes, old_genomes = genbank_status

    return genbank_status

def update(genbank_mirror, genbank_status, path_vars, assembly_summary, species_list):

    info_dir, slurm, out, logger = path_vars
    curate.create_species_dirs(genbank_mirror, assembly_summary, logger, species_list)
    local_genomes, new_genomes, old_genomes = genbank_status

    curate.remove_old_genomes(genbank_mirror, assembly_summary, old_genomes, logger)
    sync.sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes, logger)
    curate.unzip_genbank_mirror(genbank_mirror)
    rename.rename(genbank_mirror, assembly_summary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-s", "--species", help = 'List of species', nargs='+', default='all')
    parser.add_argument("-p", "--parallel", help = 'Submit jobs in parallel via SLURM pipeline tool.', action="store_true")
    parser.add_argument("-r", "--rsync", action="store_true")
    parser.add_argument("--use_local_assembly", help = 'Use local copy new assembly_summary.txt and names.dmp',\
                        action="store_true")
    args = parser.parse_args()

    update_assembly_summary=True
    if args.use_local_assembly:
        update_assembly_summary=False

    genbank_mirror = args.genbank_mirror
    path_vars, assembly_summary, species = setup(genbank_mirror, args.species, update_assembly_summary)
    info_dir, slurm, out, logger = path_vars
    genbank_status = assess_genbank_mirror(genbank_mirror, assembly_summary, species, logger)
    local_genomes, new_genomes, old_genomes = genbank_status

    curate.create_species_dirs(genbank_mirror, logger, species)
    curate.remove_old_genomes(genbank_mirror, assembly_summary, local_genomes, old_genomes, logger)
    sync.rsync_latest_genomes(genbank_mirror, assembly_summary, new_genomes)
    post_rsync_cleanup(genbank_mirror, assembly_summary, logger)
    curate.unzip_genbank_mirror(genbank_mirror)
    rename.rename(genbank_mirror, assembly_summary)

if __name__ == "__main__":
    main()
