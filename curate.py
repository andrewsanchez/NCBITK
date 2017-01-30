import os
import tarfile
import pandas as pd
from re import sub
from urllib.request import urlretrieve
from urllib.error import URLError

def get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")
    urlretrieve(assembly_summary_url, assembly_summary_dst)
    assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

def get_species_names(genbank_mirror, taxdump_url="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"):

    """
    Get names.dmp from the taxonomy dump
    """

    info_dir = os.path.join(genbank_mirror, ".info")
    taxdump = urlretrieve(taxdump_url)
    taxdump_tar = tarfile.open(taxdump[0])
    taxdump_tar.extract('names.dmp', info_dir)
    names_dmp = os.path.join(genbank_mirror, ".info", 'names.dmp')
    names = pd.read_csv(names_dmp, sep='\t|', skiprows=3, index_col=0, header=None, usecols=[0,2,6], engine='python', iterator=True)
    names = pd.concat(chunk[chunk[6] == 'scientific name'] for chunk in names)
    names.drop(6, axis=1, inplace=True)
    names.columns = ['species']
    names.index.name = 'species_taxid'
    names.species.replace({' ': '_'}, regex=True, inplace=True)
    names.species.replace({'/': '_'}, regex=True, inplace=True)
    names.to_csv(names_dmp)

    return names

def get_resources(genbank_mirror):

    """
    Get assembly_summary.txt for bacteria and taxonomy dump file.
    Parse and load into Pandas DataFrames.
    """

    assembly_summary = get_assembly_summary(genbank_mirror)
    names = get_species_names(genbank_mirror)

    return assembly_summary.head(), names

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
            print("removing old {}".format(f))

def species_list_from_taxdmp(species_taxids, names):

    species_list = []
    for species_taxid in species_taxids:
        species = names.species.loc[species_taxid]
        species_list.append(species)

    return species_list

def check_species_dirs(genbank_mirror, species_list):

    print('Checking directories for each species in complete_species_list')

    for species in species_list:
        species_dir = os.path.join(genbank_mirror, species)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)

def get_species_taxids(assembly_summary):

    species_taxids = set(assembly_summary.species_taxid.tolist())

    return species_taxids

def remove_old_genomes(genbank_mirror):

    info_dir, slurm, out = instantiate_path_vars(genbank_mirror)
    log = instantiate_log(info_dir)
    latest_assembly_versions = read_latest_assembly_versions(genbank_mirror)
    species_directories = list(set(latest_assembly_versions.index))
    for species in species_directories:
        species_dir = os.path.join(genbank_mirror, species)

        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
        latest_genome_ids = latest_assembly_versions.index[latest_assembly_versions['species'] == species].tolist()
        latest_genome_paths = latest_assembly_versions.dir[latest_assembly_versions['species'] == species].tolist()
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)
        #  latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        #  latest_genome_paths = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]


        for genome_id in local_genome_ids:
            if genome_id not in latest_genome_ids:
                fasta = glob("{}*".format(genome_id))
                os.remove(os.path.join(genbank_mirror, species, fasta[0]))

def get_local_genome_ids(species_list):

        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]

        return local_genome_ids

def check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids):

    for genome_id in local_genome_ids:
        if genome_id not in latest_genome_ids:
            fasta = glob("{}*".format(genome_id))
            os.remove(os.path.join(genbank_mirror, species, fasta[0]))

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("-r", '--remove_old', action='store_true')
    args = parser.parse_args()

    if args.remove_old:
        remove_old_genomes(args.genbank_mirror)

if __name__ == "__main__":
    main()
