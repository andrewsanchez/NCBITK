#!/usr/bin/env python

import os
import argparse
from time import strftime
from subprocess import Popen
from re import sub
from ftp_functions.ftp_functions import *
from slurm.generate_arrays import *
#from sync.sync import *

def generate_and_submit_arrays(genbank_mirror):
    clean_up(genbank_mirror)
    complete_species_list = ftp_complete_species_list()
    latest_assembly_versions_array = gen_latest_assembly_versions_array(genbank_mirror, complete_species_list)
    slurm_script = gen_latest_assembly_versions_script(genbank_mirror, latest_assembly_versions_array)
    Popen("sbatch {}".format(slurm_script), shell="True").wait()

def parallel(genbank_mirror):
    generate_and_submit_arrays(genbank_mirror)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    parser.add_argument("-p", "--parallel", action="store_true")
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    if args.parallel:
        parallel(genbank_mirror)

main()
