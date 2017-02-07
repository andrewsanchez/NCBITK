import os
import glob
import gzip

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

def create_species_dirs(genbank_mirror, assembly_summary):

    print('Checking directories for each species in complete_species_list')

    for species in set(assembly_summary.scientific_name.tolist()):
        species_dir = os.path.join(genbank_mirror, species)
        if not os.path.isdir(species_dir):
            os.mkdir(species_dir)

def get_local_genomes(genbank_mirror):

    local_genomes = []
    for root, dirs, files in os.walk(genbank_mirror):
        for f in files:
            if f.startswith('GCA'):
                genome_id = '_'.join(f.split('_')[:2])
                local_genomes.append(genome_id)

    return local_genomes

def remove_old_genomes(genbank_mirror, assembly_summary, local_genomes):

    # TODO: there might be a faster way to do this with pandas
    for genome_id in local_genomes:
        if genome_id not in assembly_summary.index.tolist():
            associated_files = glob.glob("{}/*/{}*".format(genbank_mirror, genome_id)) # globs sketch files as well
            for f in associated_files:
                os.remove(associated_files)
                print("Removed {}".format(associated_files))

def get_new_genome_list(genbank_mirror, assembly_summary, local_genomes):

    new_genomes = []
    for genome in assembly_summary.index.tolist():
        if genome not in local_genomes:
            new_genomes.append(genome)

    print("{} new genomes.".format(len(new_genomes)))

    return new_genomes

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
    unzipped.write(decoded)
    zipped.close()
    unzipped.close()
    os.remove(zipped_src)

def unzip_genbank_mirror(genbank_mirror):

    for root, files, dirs, in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("gz"):
                genome_id = "_".join(f.split("_")[:2])
                unzip_genome(root, zipped_src, genome_id)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    parser.add_argument("-r", '--remove_old', action='store_true')
    args = parser.parse_args()

    if args.remove_old:
        remove_old_genomes(args.genbank_mirror)

if __name__ == "__main__":
    main()
