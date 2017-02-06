import os
import glob
import tarfile
import subprocess
import pandas as pd
from re import sub
from urllib.request import urlretrieve
from urllib.error import URLError

def clean_up(genbank_mirror, path_vars):

    info_dir, slurm, out = path_vars
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

def get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")
    urlretrieve(assembly_summary_url, assembly_summary_dst)
    assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

def update_assembly_summary(genbank_mirror, assembly_summary, names):

    for taxid in names.index:
        scientific_name = names.scientific_name.loc[taxid]
        # get the list of indices that share the same species_taxid in assembly_summary
        ixs = assembly_summary.index[assembly_summary.species_taxid == taxid].tolist()
        assembly_summary.loc[ixs, 'scientific_name'] = scientific_name

    updated_assembly_summary = os.path.join(genbank_mirror, '.info', 'assembly_summary_updated.csv')
    assembly_summary.to_csv(updated_assembly_summary)

    return assembly_summary

def get_scientific_names(genbank_mirror, assembly_summary, taxdump_url="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"):

    """
    Get names.dmp from the taxonomy dump
    """
    info_dir = os.path.join(genbank_mirror, ".info")
    taxdump = urlretrieve(taxdump_url)
    taxdump_tar = tarfile.open(taxdump[0])
    taxdump_tar.extract('names.dmp', info_dir)
    names_dmp = os.path.join(genbank_mirror, ".info", 'names.dmp')
    sed_cmd = "sed -i '/scientific name/!d' {}".format(names_dmp) # we only want rows with the scientific name
    subprocess.Popen(sed_cmd, shell='True').wait()
    names = pd.read_csv(names_dmp, sep='\t', index_col=0, header=None, usecols=[0,2])
    names = names.loc[set(assembly_summary.species_taxid.tolist())]
    names.index.name = 'species_taxid'
    names.columns = ['scientific_name']
    names.scientific_name.replace({' ': '_'}, regex=True, inplace=True)
    names.scientific_name.replace({'/': '_'}, regex=True, inplace=True)
    names.to_csv(names_dmp)

    return names

def get_resources(genbank_mirror):

    """
    Get assembly_summary.txt for bacteria and taxonomy dump file.
    Parse and load into Pandas DataFrames.
    """

    assembly_summary = get_assembly_summary(genbank_mirror)
    names = get_scientific_names(genbank_mirror, assembly_summary)
    assembly_summary = update_assembly_summary(genbank_mirror, assembly_summary, names)

    return assembly_summary

def check_species_dirs(genbank_mirror, assembly_summary):

    print('Checking directories for each species in complete_species_list')

    for species in set(assembly_summary.scientific_name.tolist()):
        species_dir = os.path.join(genbank_mirror, species)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)
        # TODO: also remove directories without any FASTA files
        # the directories may not actually be empty because of MASH files

def get_local_genomes(genbank_mirror):

    local_genomes = []
    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.startswith('GCA'):
                genome_id = '_'.join(f.split('_')[:2])
                local_genomes.append(genome_id)

    return local_genomes

def remove_old_genomes(genbank_mirror, assembly_summary, local_genomes):

    for genome_id in local_genomes:
        if genome_id not in assembly_summary.index.tolist():
            species = assembly_summary.scientific_name[genome_id]
            fasta = glob.glob("{}/*/{}*".format(genbank_mirror, genome_id))
            os.remove(fasta[0])
            print("Removed {}".format(fasta[0]))

def get_new_genome_list(genbank_mirror, assembly_summary, local_genomes):

    new_genomes = []
    for genome in assembly_summary.index.tolist():
        if genome not in local_genomes:
            new_genomes.append(genome)

    print("{} new genomes.".format(len(new_genomes)))

    return new_genomes

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("-r", '--remove_old', action='store_true')
    args = parser.parse_args()

    if args.remove_old:
        remove_old_genomes(args.genbank_mirror)

if __name__ == "__main__":
    main()
