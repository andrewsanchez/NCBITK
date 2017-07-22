import os
import gzip
import re
import logging
import shutil
import pandas as pd
from io import TextIOWrapper


def get_species(assembly_summary, species):

    if not species:
        species = assembly_summary.scientific_name[
            assembly_summary.scientific_name.notnull()]
        species = set(species.tolist())
        return species

    elif type(species) is tuple:
        return species

    elif type(species) is TextIOWrapper:
        species = [name.strip() for name in species]
        return species


def create_species_dirs(genbank_mirror, logger, species_list):

    for species in species_list:
        try:
            species_dir = os.path.join(genbank_mirror, species)
        except TypeError:
            continue
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)


def parse_genome_id(genome):

    genome_id = re.match('GCA_\d+\.\d', genome)

    return genome_id


def diff(a, b):

    diff = set(a) - set(b)

    return list(diff)


def get_local_genomes(genbank_mirror):

    local_genomes = {}

    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if re.match('GCA.*fasta', f):
                genome_id = parse_genome_id(f).group(0)
                genome_path = os.path.join(root, f)
                local_genomes[genome_id] = genome_path

    return local_genomes


def get_latest_assembly_versions(assembly_summary, species_list):

    latest_assembly_versions = assembly_summary.index[
        assembly_summary.scientific_name.isin(species_list)]

    return latest_assembly_versions.tolist()


def get_new_genomes(latest_assembly_versions, local_genomes):

    new_genomes = diff(latest_assembly_versions, list(local_genomes.keys()))

    return new_genomes


def get_old_genomes(local_genomes, latest_assembly_versions):

    old_genomes = diff(list(local_genomes.keys()), latest_assembly_versions)

    return old_genomes


def assess_genbank_mirror(genbank_mirror, assembly_summary, species_list,
                          logger):

    local_genomes = get_local_genomes(genbank_mirror)
    latest_assembly_versions = get_latest_assembly_versions(
        assembly_summary, species_list)
    new_genomes = get_new_genomes(latest_assembly_versions, local_genomes)
    old_genomes = get_old_genomes(local_genomes, latest_assembly_versions)

    logger.info(
        "{} genomes present in local collection.".format(len(local_genomes)))
    logger.info(
        "{} genomes missing from local collection.".format(len(new_genomes)))
    logger.info("{} old genomes to be removed.".format(len(old_genomes)))
    if not new_genomes:
        logger.info(
            "Local collection is up to date with assembly_summary.txt.")

    return local_genomes, new_genomes, old_genomes


def remove_old_genomes(genbank_mirror, assembly_summary, local_genomes,
                       old_genomes, logger):

    for genome_id in old_genomes:
        genome_path = local_genomes[genome_id]
        os.remove(genome_path)
        logger.info("Removed {}".format(genome_id))


def unzip_genome(root, f, genome_id):
    """
    Decompress genome and remove the compressed genome.
    """

    zipped_src = os.path.join(root, f)
    zipped = gzip.open(zipped_src)
    decoded = zipped.read()
    unzipped = "{}.fasta".format(genome_id)
    unzipped = os.path.join(root, unzipped)
    unzipped = open(unzipped, "wb")
    zipped.close()
    os.remove(zipped_src)
    unzipped.write(decoded)
    unzipped.close()


def unzip_genbank(genbank_mirror):

    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("gz"):
                genome_id = "_".join(f.split("_")[:2])
                try:
                    unzip_genome(root, f, genome_id)
                except OSError:
                    continue


def post_rsync_cleanup(genbank_mirror, assembly_summary, logger):

    incoming = os.path.join(genbank_mirror, 'incoming')
    for root, dirs, files in os.walk(incoming):
        for f in files:
            accession = '_'.join(f.split('_')[:2])
            try:
                species = assembly_summary.scientific_name.loc[accession]
            except KeyError:
                logger.info('KeyError for {}'.format(accession))
                continue

            src = os.path.join(root, f)
            dst = os.path.join(genbank_mirror, species, f)
            shutil.move(src, dst)

    shutil.rmtree(incoming)


def rm_duplicates(seq):
    """
    remove duplicate strings during renaming
    """

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def clean_up_name(name):
    rm_words = re.compile(
        r'((?<=_)(sp|sub|substr|subsp|str|strain)(?=_))')
    name = rm_words.sub('_', name)
    name = re.sub('_+', '_', name)
    name = rm_duplicates(name.split('_'))
    name = '_'.join(name)
    return name


def rename_genome(genome, assembly_summary):
    """
    Rename FASTA's.
    """

    genome_id = parse_genome_id(genome).group(0)
    if genome_id in assembly_summary.index:
        scientific_name = assembly_summary.get_value(
            genome_id, 'scientific_name')
        infraspecific_name = assembly_summary.get_value(
            genome_id, 'infraspecific_name')
        organism_name = assembly_summary.get_value(
            genome_id, 'organism_name')
        if type(infraspecific_name) == float:
            infraspecific_name = ''
        isolate = assembly_summary.get_value(
            genome_id, 'isolate')
        if type(isolate) == float:
            isolate = ''
        assembly_level = assembly_summary.get_value(
            genome_id, 'assembly_level')
        name = '{}_{}_{}_{}_{}_{}.fasta'.format(
            genome_id, organism_name, scientific_name, infraspecific_name, isolate, assembly_level)
        name = clean_up_name(name)
        return name

def rename_genbank(target_dir, assembly_summary):

    for root, dirs, files in os.walk(target_dir):
        genomes = [f for f in files if re.match('GCA.*fasta', f)]
        for genome in genomes:
            name = rename_genome(genome, assembly_summary)
            if name:
                old = os.path.join(root, genome)
                new = os.path.join(root, name)
                os.rename(old, new)
