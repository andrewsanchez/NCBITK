def read_latest_assembly_versions(genbank_mirror, ix_col=1):

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=ix_col, header=None)
    if ix_col == 1:
        latest_assembly_versions.columns = ['species', 'dir']
    elif ix_col == 0:
        latest_assembly_versions.columns = ['id', 'dir']

    return latest_assembly_versiondef get_new_genomes(genbank_mirror, latest_assembly_versions):

    new_genomes = []
    species_directories = list(set(latest_assembly_versions['species']))
    for species in species_directories:
        species_dir = os.path.join(genbank_mirror, species)


#  def get_new_genome_list(latest_assembly_versions):
    
    #  ids_and_paths = get_ids_and_paths(latest_assembly_versions)
    #  new_genomes = []
    #  for info in ids_and_paths:
        #  genome_id = info[0]
        #  genome_path = info[1]
        #  if genome_id not in local_genome_ids:
            #  new_genomes.append(genome_id)

    #  return new_genomes

def ftp_complete_species_list():

    """Connect to NCBI's ftp site and retrieve complete list of bacteria."""

    ftp = ftp_login()

    print("Getting list of all bacteria directories.")
    print("Estimated wait time:  1 minute.")
    try:
        complete_species_list = ftp.nlst()
    except error_temp:
        sleep(30)
        complete_species_list = ftp.nlst()

    print("All bacteria directories succesfully read from ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/")

    return complete_species_list

def write_latest_assembly_versions(genbank_mirror, species, ftp):

    latest_dir = os.path.join(species, "latest_assembly_versions")
    latest_dirs = [accession.split("/")[-1] for accession in ftp.nlst(latest_dir)]
    accession_ids = ["_".join(id.split("_")[:2]) for id in latest_dirs]
    latest_assembly_versions_list = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    dirs_and_ids = zip(latest_dirs, accession_ids)

    with open(latest_assembly_versions_list, "a") as f:
        for item in dirs_and_ids:
            genome_path = item[1]
            accession_id = item[0]
            f.write("{},{},{}\n".format(species, genome_path, accession_id))

def get_latest_assembly_versions(genbank_mirror, complete_species_list, genbank_stats, ymdt):

    """
    Create DataFrame to represent genbank's directory structure.
    """

    ftp = ftp_login()

    print("Getting latest assembly versions for {} species.".format(len(complete_species_list)))
    print("Estimated wait time:  2 hours.")

    for species in complete_species_list:
        try:
            write_latest_assembly_versions(genbank_mirror, species, ftp)
        except error_temp:
            log_error("latest", species, genbank_stats, ymdt)
        except BrokenPipeError:
            try:
                ftp = ftp_login()
                write_latest_assembly_versions(genbank_mirror, species, ftp)
            except error_temp:
                log_error("latest", species, genbank_stats, ymdt)

    latest_assembly_versions = os.path.join(genbank_mirror, ".info", "latest_assembly_versions.csv")
    latest_assembly_versions = pd.read_csv(latest_assembly_versions, index_col=0, header=None)
    latest_assembly_versions.columns = ["id", "dir"]

    return latest_assembly_versions

def check_species_dirs(genbank_mirror, species):

    species_dir = os.path.join(genbank_mirror, species)
    if not os.path.isdir(species_dir):
        os.mkdir(species_dir)

    return species_dir

def grab_and_organize_genomes(genbank_mirror, genbank_stats, latest_assembly_versions):

    species_directories = os.listdir(genbank_mirror)

    for species in species_directories:
        local_genome_ids = get_local_genome_ids(species_dir)
        latest_genome_ids = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["id"]].values.tolist()]
        latest_genome_paths = [sub("[\[\]']", "", str(i)) for i in latest_assembly_versions.loc[species, ["dir"]].values.tolist()]
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)
        print("species", species)
        print("local_genome_ids", local_genome_ids)
        print("latest_genome_ids", latest_genome_ids)
        print("latest_genome_paths", latest_genome_paths)

        check_local_genomes(genbank_mirror, species, local_genome_ids, latest_genome_ids, genbank_stats)
        sync_latest_genomes(genbank_mirror, species, local_genome_ids, ids_and_paths, genbank_stats)

def log_error(msg, species, genbank_stats, ymdt):
    if msg == "latest":
        msg = "no latest_assembly_versions dir for"
        with open(genbank_stats, "a") as stats:
            stats.write("{} - {} {}\n".format(species, msg, ymdt))
