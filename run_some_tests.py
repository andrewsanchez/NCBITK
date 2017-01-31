#!/usr/bin/env python

# from NCBITK import instantiate_path_vars
import NCBITK.config as config

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_mirror", help = "Directory to save fastas", type=str)
    args = parser.parse_args()

    genbank_mirror = args.genbank_mirror
    # path_vars = instantiate_path_vars(genbank_mirror)
    path_vars = config.instantiate_path_vars(genbank_mirror)
    print(path_vars)

main()
