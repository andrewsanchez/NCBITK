import argparse
import os
import glob
import gzip
import re
import logging
import shutil
import pandas as pd


def clean_up(genbank_mirror, path_vars):

    info_dir, slurm, out, logger = path_vars

    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions.csv")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sbatch")
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    sync_array_script = os.path.join(slurm, 'sync_array_script.sbatch')
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')

    for f in [latest_assembly_versions, latest_assembly_versions_array, slurm_script, sync_array_script, grab_genomes_script, sync_array]:
        if os.path.isfile(f):
            os.remove(f)

def get_species_list(assembly_summary, species_list):

    if species_list == "all":

        species_list = assembly_summary.scientific_name[assembly_summary.scientific_name.notnull()]
        species_list = set(species_list.tolist())

    return species_list

def create_species_dirs(genbank_mirror, assembly_summary, logger, species_list):

    for species in species_list:
        try:
            species_dir = os.path.join(genbank_mirror, species)
        except TypeError:
            continue
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
            logger.info("Directory created: {}".format(species))

def get_local_genomes(genbank_mirror):

    local_genomes = []

    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.startswith('GCA'):
                genome_id = '_'.join(f.split('_')[:2])
                local_genomes.append(genome_id)

    return local_genomes

def get_new_genome_list(genbank_mirror, assembly_summary, local_genomes, species_list):

    # TODO: Faster way to do this in pandas?

    new_genomes = []

    for species in species_list:
        latest_assembly_versions = assembly_summary.index[assembly_summary.scientific_name == species].tolist()
        for genome in latest_assembly_versions:
            if genome not in local_genomes:
                new_genomes.append(genome)

    return new_genomes

def remove_old_genomes(genbank_mirror, assembly_summary, old_genomes, logger):

    # TODO: there might be a faster way to do this with pandas
    for genome_id in old_genomes:
        # Would have to keep the old assembly summary file in order to avoid globbing the species dir
        associated_files = glob.glob("{}/*/{}*".format(genbank_mirror, genome_id)) # globs sketch files as well
        for f in associated_files:
            os.remove(f)
            logger.info("Removed {}".format(f))

def get_sketch_files(genbank_mirror):

    sketch_files = []
    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.endswith('msh'):
                genome_id = re.sub(r'.msh', '', f)
                sketch_files.append(genome_id)

    return sketch_files

def get_missing_sketch_files(local_genomes, new_genomes, sketch_files):

    missing_sketch_files = []

    for genome_id in local_genomes:
        if genome_id not in sketch_files:
            missing_sketch_files.append(genome_id)

    for genome_id in new_genomes:
        missing_sketch_files.append(genome_id)

    return missing_sketch_files

def get_old_genomes(genbank_mirror, assembly_summary, local_genomes):

    old_genomes = []

    latest_assembly_versions = assembly_summary.index.tolist()
    for genome_id in local_genomes:
        if genome_id not in latest_assembly_versions:
            old_genomes.append(genome_id)

    return old_genomes

def assess_genbank_mirror(genbank_mirror, assembly_summary, species_list):

    local_genomes = get_local_genomes(genbank_mirror)
    new_genomes = get_new_genome_list(genbank_mirror, assembly_summary, local_genomes, species_list)
    old_genomes = get_old_genomes(genbank_mirror, assembly_summary, local_genomes)
    sketch_files = get_sketch_files(genbank_mirror)
    missing_sketch_files = get_missing_sketch_files(local_genomes, new_genomes, sketch_files)

    return local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files

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

    # shutil.rmtree(incoming)


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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    path_vars = config.instantiate_path_vars(genbank_mirror)
    clean_up(genbank_mirror, path_vars)

if __name__ == "__main__":
    main()
