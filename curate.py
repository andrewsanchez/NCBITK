import os
import pandas as pd
from re import sub

def instantiate_path_vars(genbank_mirror):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    out = os.path.join(slurm, "out")

    return info_dir, slurm, out

def instantiate_log(info_dir):

    log = os.path.join(info_dir, 'genbank.log')

    try:
        log = pd.read_csv(log, index_col=0)
    except OSError:
        log = pd.DataFrame(columns=["Status", "Date"])

    return record

def instantiate_df(path, cols):

    df = pd.read_csv(path, index_col=0, header=None)
    df.columns = cols

    return df


def clean_up(genbank_mirror):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    for d in [genbank_mirror, info_dir, slurm, out]:
        if not os.path.isdir(d):
            os.mkdir(d)

    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions.csv")
    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt")
    slurm_script = os.path.join(slurm, "get_latest_assembly_versions.sbatch")
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    sync_array_script = os.path.join(slurm, 'sync_array_script.sbatch')
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')
    for f in [latest_assembly_versions, latest_assembly_versions_array, slurm_script, sync_array_script, grab_genomes_script, sync_array]:
        if os.path.isfile(f):
            os.remove(f)

def check_species_dirs(genbank_mirror, complete_species_list):

    print('Checking directories for each species in complete_species_list')

    for species in complete_species_list:
        species_dir = os.path.join(genbank_mirror, species)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)

def remove_old_genomes(genbank_mirror):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    log = instantiate_log(info_dir)
    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0)
    latest_assembly_versions.columns = ["id", "dir"]
    species_directories = list(set(latest_assembly_versions.index))
    for species in species_directories:
        species_dir = os.path.join(genbank_mirror, species)
        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
        latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        latest_genome_paths = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)
        check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids)

def read_latest_assembly_versions(genbank_mirror, ix_col=1):
    
    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=ix_col, header=None)
    if ix_col == 1:
        latest_assembly_versions.columns = ['species', 'dir']
    elif ix_col == 0:
        latest_assembly_versions.columns = ['id', 'dir']

    return latest_assembly_versions

def get_ids_and_paths(latest_assembly_versions):

    species_directories = list(set(latest_assembly_versions.index))

    for species in species_directories:
        species_dir = os.path.join(genbank_mirror, species)
        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
        latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        latest_genome_paths = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)

        return ids_and_paths

def get_new_genome_list(latest_assembly_versions):
    
    ids_and_paths = get_ids_and_paths(latest_assembly_versions)
    new_genomes = []
    for info in ids_and_paths:
        genome_id = info[0]
        genome_path = info[1]
        if genome_id not in local_genome_ids:
            new_genomes.append(genome_id)

    return new_genomes

def check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids):

    for genome_id in local_genome_ids:
        if genome_id not in latest_genome_ids:
            fasta = glob("{}*".format(genome_id))
            os.remove(os.path.join(genbank_mirror, species, fasta[0]))
