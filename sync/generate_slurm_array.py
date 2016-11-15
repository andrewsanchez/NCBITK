#!/usr/bin/env python

import os
from re import sub
from global_functions import *

def generate_dir_structure_array(all_species):

    ftp = ftp_login()
    with open("get_dir_structure_array.txt", "a") as f:
        for species in all_species:
            f.write('python /common/contrib/tools/NCBITK/organize_genomes.py -t -s {}\n'.format(species))

def get_dir_tree(species):

    ftp = ftp_login()
    latest = "{}/latest_assembly_versions".format(species)
    with open("{}_dir_structure.csv".format(species), "a") as f:
        try:
            ids = ftp.nlst(latest)
            for name in ids:
                id = name.split("/")[-1]
                f.write("{},{}\n".format(species, id))
        except error_temp:
            f.write("{} doesn't have a latest_assembly_versions/ directory and will be skipped".format(species))

def main():

    import argparse 

    parser = argparse.ArgumentParser(description = "Get the latest assembly versions for each species.")
    parser.add_argument("-s", "--species", type=str)
    parser.add_argument("-a", "--array", help="Generate slurm array commands.", action="store_true")
    parser.add_argument("-t", "--tree", help="Get the directory tree structure of */latest_assembly_versions", action="store_true")
    args = parser.parse_args()

    all_species = ftp_complete_species_list()
    if args.array:
        generate_dir_structure_array(all_species)
    elif args.tree:
        get_dir_tree(args.species)

main()
