#!/usr/bin/env python

import os
import argparse
from time import strftime
from re import sub
from ftp_functions import *

ymd = strftime("%y.%m.%d")

def get_latest_assembly_versions(genbank_mirror, species):

    info_dir = os.path.join(genbank_mirror, ".info")
    slurm = os.path.join(info_dir, "slurm")
    latest_assembly_versions = os.path.join(info_dir, "latest_assembly_versions_{}.txt".format(ymd))
    ftp = ftp_login()
    latest_dir = os.path.join(species, "latest_assembly_versions")
    with open(latest_assembly_versions, "a") as f:
        complete_ids = [complete_id.split("/")[-1] for complete_id in ftp.nlst(latest_dir)]
        short_ids = ["_".join(accession_id.split("_")[:2]) for accession_id in complete_ids]
        complete_and_short = zip(complete_ids, short_ids)
        for item in complete_and_short:
            complete_id = item[0]
            short_id = item[1]
            f.write("{},{},{}\n".format(species, short_id, complete_id))
def main():

        parser = argparse.ArgumentParser(description = "Get the latest assembly versions for each species.")
        parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
        parser.add_argument("species", type=str)
        args = parser.parse_args()
        get_latest_assembly_versions(args.genbank_mirror, args.species)

main()

