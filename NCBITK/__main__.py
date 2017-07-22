import argparse
import logging
import os

from NCBITK import config, curate, get_resources, sync


def setup(genbank_mirror, species_list, update_assembly_summary):

    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror,
                                                   update_assembly_summary)
    species = curate.get_species_list(assembly_summary, species_list)
    genbank_status = curate.assess_genbank_mirror(
        genbank_mirror, assembly_summary, species, logger)

    return path_vars, assembly_summary, species, genbank_status


def update(genbank_mirror, genbank_status, path_vars, assembly_summary,
           species_list):

    info_dir, slurm, out, logger = path_vars
    curate.create_species_dirs(genbank_mirror, assembly_summary, logger,
                               species_list)
    local_genomes, new_genomes, old_genomes = genbank_status

    curate.remove_old_genomes(genbank_mirror, assembly_summary, old_genomes,
                              logger)
    sync.sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes,
                             logger)
    curate.unzip_genbank_mirror(genbank_mirror)
    rename.rename(genbank_mirror, assembly_summary)


def show_genbank_status(genbank_status):

    local_genomes, new_genomes, old_genomes = genbank_status
    print('{} local genome(s)'.format(len(local_genomes)))
    print('{} new genome(s)'.format(len(new_genomes)))
    print('{} old genome(s)'.format(len(old_genomes)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'genbank_mirror', help='Directory to save genomes to', type=str)
    parser.add_argument(
        '-s',
        '--species',
        help='List of species to download genomes for',
        nargs='+',
        default='all')
    parser.add_argument('-u', '--update', action='store_true')
    parser.add_argument(
        '--use_local_assembly',
        help='Use local assembly_summary and taxonomy dump',
        action='store_true')
    parser.add_argument(
        '--status',
        help='Show the current status of your genome collection',
        action='store_true')
    args = parser.parse_args()

    if args.use_local_assembly:
        update_assembly_summary = False
    else:
        update_assembly_summary = True

    genbank_mirror = args.genbank_mirror
    path_vars, assembly_summary, species, genbank_status = setup(
        genbank_mirror, args.species, update_assembly_summary)
    info_dir, slurm, out, logger = path_vars
    local_genomes, new_genomes, old_genomes = genbank_status

    if args.status:
        show_genbank_status(genbank_status)

    if args.update:
        curate.create_species_dirs(genbank_mirror, logger, species)
        curate.remove_old_genomes(genbank_mirror, assembly_summary,
                                  local_genomes, old_genomes, logger)
        sync.rsync_latest_genomes(genbank_mirror, assembly_summary,
                                  new_genomes)
        curate.post_rsync_cleanup(genbank_mirror, assembly_summary, logger)
        curate.unzip_genbank_mirror(genbank_mirror)
        curate.rename(genbank_mirror, assembly_summary)


if __name__ == "__main__":
    main()
