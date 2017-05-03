import os
import glob
import gzip
import re
import logging
import shutil
import pandas as pd

def get_species_list(assembly_summary, species_list):

    if species_list == "all":

        species_list = assembly_summary.scientific_name[assembly_summary.scientific_name.notnull()]
        species_list = set(species_list.tolist())

        return species_list

    elif type(species_list) is list:

        return species_list

    elif type(species_list) is str:

        return [species_list]

def create_species_dirs(genbank_mirror, logger, species_list):

    for species in species_list:
        try:
            species_dir = os.path.join(genbank_mirror, species)
        except TypeError:
            continue
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
            logger.info("Directory created: {}".format(species))

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

    latest_assembly_versions = assembly_summary.index[assembly_summary.scientific_name.isin(species_list)]

    return latest_assembly_versions.tolist()

def get_new_genome_list(latest_assembly_versions, local_genomes):

    new_genomes = diff(latest_assembly_versions, list(local_genomes.keys()))

    return new_genomes

def get_old_genomes(local_genomes, latest_assembly_versions):

    old_genomes = diff(list(local_genomes.keys()), latest_assembly_versions)

    return old_genomes

def assess_genbank_mirror(genbank_mirror, assembly_summary, species_list, logger):

    local_genomes = get_local_genomes(genbank_mirror)
    latest_assembly_versions = get_latest_assembly_versions(assembly_summary, species_list)
    new_genomes = get_new_genome_list(latest_assembly_versions, local_genomes)
    old_genomes = get_old_genomes(local_genomes, latest_assembly_versions)

    logger.info("{} genomes present in local collection.".format(len(local_genomes)))
    logger.info("{} genomes missing from local collection.".format(len(new_genomes)))
    logger.info("{} old genomes to be removed.".format(len(old_genomes)))
    if not new_genomes:
        logger.info("Local collection is up to date with assembly_summary.txt.")

    return local_genomes, new_genomes, old_genomes

def remove_old_genomes(genbank_mirror, assembly_summary, local_genomes, old_genomes, logger):

    for genome_id in old_genomes:
        genome_path = local_genomes[genome_id]
        os.remove(genome_path)
        logger.info("Removed {}".format(f))

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

def unzip_genbank_mirror(genbank_mirror):

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
            accession = '_'.join(f.split('_') [:2])
            try:
                species = assembly_summary.scientific_name.loc[accession]
            except KeyError:
                logger.info('KeyError for {}'.format(accession))
                continue

            src = os.path.join(root, f)
            dst = os.path.join(genbank_mirror, species, f)
            shutil.move(src, dst)

    shutil.rmtree(incoming)


def rename(target_dir, assembly_summary):

    """
    Clean up assembly_summary.txt and renamed FASTA's.
    """

    def rm_duplicates(seq):

        """
        remove duplicate strings during renaming
        """

        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    # If infraspecific_name and isolate columns are empty, fill infraspecific_name with "NA"
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].isnull())].fillna('NA'))

    # If infraspecific_name column is empty and isolate column is not empty, fill infraspecific_name with the value of isolate.
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].notnull())].fillna(assembly_summary['isolate']))

    assembly_summary.assembly_level.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({'[\W]': '_'}, regex=True, inplace=True)
    assembly_summary.infraspecific_name.replace({'[\W]': '_'}, regex=True, inplace=True)

    for root, dirs, files in os.walk(target_dir):
        for f in files:
            if f.startswith("GCA"):
                accession_id = '_'.join(f.split('_') [:2])
                if accession_id in assembly_summary.index:
                    org_name = assembly_summary.get_value(accession_id, 'scientific_name')
                    strain = assembly_summary.get_value(accession_id, 'infraspecific_name')
                    assembly_level  = assembly_summary.get_value(accession_id, 'assembly_level')
                    new_name = '{}_{}_{}_{}.fasta'.format(accession_id, org_name, strain, assembly_level)
                    rm_words = re.compile( r'((?<=_)(sp|sub|substr|subsp|str|strain)(?=_))' )
                    new_name = rm_words.sub('_', new_name)
                    new_name = re.sub(r'_+', '_', new_name )
                    new_name = new_name.split('_')
                    new_name = rm_duplicates(new_name)
                    new_name = '_'.join(new_name)
                    old = os.path.join(root, f)
                    new = os.path.join(root, new_name)
                    os.rename(old, new)
