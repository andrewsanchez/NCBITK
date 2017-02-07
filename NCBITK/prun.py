import os
import argparse
import subprocess
import pandas as pd
from re import sub
from time import strftime, sleep
import curate
from ftp_functions import ftp_functions

def gen_array(array, cmd, in_list, groups):

    with open(array, 'a') as f:
        for group in groups:
            group = ' '.join(group)
            f.write("{} {}\n".format(cmd, group))

    array_len = len(list(open(array)))

    return array_len

def get_groups(in_list, group_size):

    groups = [in_list[n:n+group_size] for n in range(0, len(in_list), group_size)]

    return groups

def gen_sbatch_dependent(path_vars, dp_id, job_name, cmd, time='01:00'):

    info_dir, slurm, out = path_vars
    sbatch_dependent_script = os.path.join(slurm, '{}.sbatch'.format(job_name))
    out = os.path.join(out, '{}_%a.out'.format(job_name))
    print('Generating {}'.format(sbatch_dependent_script))

    with open(sbatch_dependent_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time={}\n".format(time))
        f.write("#SBATCH --job-name={}\n".format(job_name))
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(dp_id))
        f.write('cmd="{}"\n'.format(cmd))
        f.write("srun $cmd")
    
    return sbatch_dependent_script

def gen_sbatch_array_script(path_vars, array, array_len, job_name, time='01:00'):

    info_dir, slurm, out = path_vars
    out = os.path.join(out, "{}_%a.out".format(job_name))
    sbatch_script = os.path.join(slurm, "{}.sbatch".format(job_name))
    print('Generating {}'.format(sbatch_script))
    with open(sbatch_script, "a") as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time={}\n".format(time))
        f.write("#SBATCH --job-name={}\n".format(job_name))
        f.write("#SBATCH --output={}\n".format(out)) # can I remove this line to avoid getting out files?
        f.write("#SBATCH --array=1-{}%2\n".format(array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(array))
        f.write("srun $cmd")

    return sbatch_script 

def submit_sbatch(slurm_script):

    job_id, errs = subprocess.Popen("sbatch {}".format(slurm_script), # submit array
            shell="True",  stdout=subprocess.PIPE, universal_newlines=True).communicate()
    print(job_id, end='')
    job_id = job_id.split(" ")[-1].strip()

    return job_id

def submit_salloc(cmd, dependency_id, time='20:00'):
    subprocess.Popen('salloc -t {} --dependency={} {}'.format(time, dependency_id, cmd), shell=True)

def get_ids_and_paths(genbank_mirror, latest_assembly_versions):

    species_directories = list(set(latest_assembly_versions['species']))

    local = []
    latest_ids = []
    latest_paths = []
    for species in species_directories:
        species_dir = os.path.join(genbank_mirror, species)
        for genome_id in os.listdir(species_dir):
            genome_id = "_".join(genome_id.split("_")[:2])
            local.append(genome_id)
        local_genome_ids = ["_".join(genome_id.split("_")[:2]) for genome_id in os.listdir(species_dir)]
        latest_genome_ids = latest_assembly_versions.index[latest_assembly_versions['species'] == species].tolist()
        latest_genome_paths = latest_assembly_versions.dir[latest_assembly_versions['species'] == species].tolist()
        ids_and_paths = zip(latest_genome_ids, latest_genome_paths)

        return local_genome_ids, ids_and_paths

def get_new_genome_list(genbank_mirror, latest_assembly_versions):
    
    local_genome_ids, ids_and_paths = get_ids_and_paths(genbank_mirror, latest_assembly_versions)
    new_genomes = []
    for info in ids_and_paths:
        genome_id = info[0]
        genome_path = info[1]
        if genome_id not in local_genome_ids:
            new_genomes.append(genome_id)

    return new_genomes

ymd = strftime("%y.%m.%d")

def gen_grab_genomes_script(genbank_mirror, grab_genomes_array_job_id):

    info_dir, slurm, out = curate.instantiate_path_vars(genbank_mirror)
    grab_genomes_array = os.path.join(slurm, "grab_genomes_array.txt")
    grab_genomes_script = os.path.join(slurm, 'grab_genomes_script.sbatch')
    out = os.path.join(out, 'grab_genomes%a.out')
    print('Generating {}'.format(grab_genomes_script))

    # Make sure generation of array is complete
    while not os.path.isfile(grab_genomes_array):
        sleep(10)
    else:
        grab_genomes_array_len = len(list(open(grab_genomes_array)))

    with open(grab_genomes_script, 'a') as f:
        f.write("#!/bin/sh\n")
        f.write("#SBATCH --time=180:00\n")
        f.write("#SBATCH --job-name=grab_genomes\n")
        f.write("#SBATCH --output={}\n".format(out))
        f.write("#SBATCH --dependency={}\n".format(grab_genomes_array_job_id))
        f.write("#SBATCH --array=1-{}%2\n".format(grab_genomes_array_len))
        f.write('cmd=$(sed -n "$SLURM_ARRAY_TASK_ID"p "{}")\n'.format(grab_genomes_array))
        f.write("srun $cmd")
    
    return grab_genomes_script, grab_genomes_array_len

def write_grab_genomes_array(genbank_mirror):

    latest_assembly_versions = curate.read_latest_assembly_versions(genbank_mirror)
    sync_array = os.path.join(genbank_mirror, ".info", "slurm", "sync_array.txt")
    print('Generating {}'.format(sync_array))
    new_genomes = get_new_genomes(genbank_mirror, latest_assembly_versions)
    args = []
    for id in new_genomes:
        species = latest_assembly_versions.loc[id, 'species']
        path = latest_assembly_versions.loc[id, 'dir']
        args.append(','.join([species, id, path]))

    groups = [args[n:n+2000] for n in range(0, len(args), 2000)] 
    with open(sync_array, "a") as f:
        for group in groups:
            f.write("python /common/contrib/tools/NCBITK/ftp_functions/ftp_functions.py -g {} {}\n".format(genbank_mirror, ' '.join(group)))

    return sync_array

def get_latest(genbank_mirror, path_vars):

    '''
    Get the contents of each species' latest_assembly_versions directory in order
    to map the accession IDs of every genome to their respective species directory.
    '''

    info_dir, slurm, out = path_vars
    complete_species_list = ftp_functions.ftp_complete_species_list()#[:30]
    curate.create_species_dirs(genbank_mirror, complete_species_list)

    latest_assembly_versions_array = os.path.join(slurm, "latest_assembly_versions_array.txt") # should I unpack these paths as in instantiate_path_vars?

    array_len = gen_array(latest_assembly_versions_array,
            'python /common/contrib/tools/NCBITK/ftp_functions/get_latest_assembly_versions.py {}'.format(genbank_mirror),
            complete_species_list, get_groups(complete_species_list, 1000))

    latest_assembly_versions_script = gen_sbatch_array_script(path_vars, latest_assembly_versions_array, array_len, 'get_latest', time='90:00')
    get_latest_job_id = submit_sbatch(latest_assembly_versions_script)
    get_latest_job_id = '{}_{}'.format(get_latest_job_id, array_len)

    return get_latest_job_id

def update_genomes(genbank_mirror, path_vars, get_latest_job_id):
    
    info_dir, slurm, out = path_vars
    gen_grab_genomes_sbatch = gen_sbatch_dependent(path_vars, get_latest_job_id,
            'gen_grab_genomes_array',
            'python /common/contrib/tools/NCBITK/generate_arrays.py -g {}'.format(genbank_mirror),
            time = '15:00')

    # submit the sbatch script that creates the array containing commands for grabbing genomes
    gen_grab_genomes_array_job_id = submit_sbatch(gen_grab_genomes_sbatch)
    grab_genomes_script, grab_genomes_array_len = gen_grab_genomes_script(genbank_mirror, gen_grab_genomes_array_job_id)
    job_id = submit_sbatch(grab_genomes_script)
    #  remove_old_genomes(genbank_mirror)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror")
    args = parser.parse_args()

    write_grab_genomes_array(args.genbank_mirror)

if __name__ == "__main__":
    main()
