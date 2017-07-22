import logging
import os
import click

from NCBITK import config, curate, get_resources, sync


def setup(genbank_mirror, species, update_assembly_summary):
    path_vars = config.instantiate_path_vars(genbank_mirror)
    info_dir, slurm, out, logger = path_vars
    assembly_summary = get_resources.get_resources(genbank_mirror,
                                                   update_assembly_summary)
    species = curate.get_species(assembly_summary, species)
    genbank_status = curate.assess_genbank_mirror(
        genbank_mirror, assembly_summary, species, logger)

    return path_vars, assembly_summary, species, genbank_status


def update(genbank_mirror, genbank_status, path_vars, assembly_summary,
           species):

    info_dir, slurm, out, logger = path_vars
    curate.create_species_dirs(genbank_mirror, assembly_summary, logger,
                               species)
    local_genomes, new_genomes, old_genomes = genbank_status

    curate.remove_old_genomes(genbank_mirror, assembly_summary, old_genomes,
                              logger)
    sync.sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes,
                             logger)
    curate.unzip_genbank(genbank_mirror)
    rename.rename(genbank_mirror, assembly_summary)


def show_genbank_status(genbank_status):

    local_genomes, new_genomes, old_genomes = genbank_status
    print('{} local genome(s)'.format(len(local_genomes)))
    print('{} new genome(s)'.format(len(new_genomes)))
    print('{} old genome(s)'.format(len(old_genomes)))


@click.command()
@click.option('--update/--no-update',
              help='Sync your collection with '
              'the latest assembly versions',
              default=True)
@click.option('--update-assembly/--local-assembly',
              help='Download the latest assembly summary and taxonomy dump'
              'Or use your local copies.',
              default=True)
@click.option('--from-file', type=click.File('r'))
@click.option('--status',
              help='Show the current status of your genome collection',
              is_flag=True,
              default=False)
@click.argument('genbank')
@click.argument('species', nargs=-1, required=False)
def main(update, update_assembly, from_file, status, genbank, species):
    if from_file:
        species = from_file
    path_vars, assembly_summary, species, genbank_status = setup(
        genbank, species, update_assembly)
    info_dir, slurm, out, logger = path_vars
    local_genomes, new_genomes, old_genomes = genbank_status
    if status:
        show_genbank_status(genbank_status)
    if update:
        curate.create_species_dirs(genbank, logger, species)
        curate.remove_old_genomes(genbank, assembly_summary,
                                  local_genomes, old_genomes, logger)
        sync.rsync_latest_genomes(genbank, assembly_summary,
                                  new_genomes)
        curate.post_rsync_cleanup(genbank, assembly_summary, logger)
        curate.unzip_genbank(genbank)
        curate.rename_genbank(genbank, assembly_summary)

if __name__ == '__main__':
    main()
