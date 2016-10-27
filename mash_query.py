#!/usr/bin/env python

import os, argparse
import pandas as pd
from subprocess import Popen

def query_genome(sketch_files, reference_genome, threshold,  mash_exe):
    output = reference_genome.strip(".fasta.msh")
    output = "{}_distances.csv".format(output)
    output = os.path.join(sketch_files, output)
    distance_cmd = "{} dist -t {} *msh > {}".format(mash_exe, reference_genome, output)
    Popen(distance_command, shell="True").wait()
    distances = pd.read_csv(output, delimiter="\t", index_col=0, header=0)
    passed = distances[distances <=  threshold]
    for genome in passed.index:


def main():
    parser = argparse.ArgumentParser(description = "Identify genomes within a provided distance from a single reference genome")
    parser.add_argument("sketch_files", help = "directory containing MASH sketch files")
    parser.add_argument("reference_genome", help = "directory containing MASH sketch files")
    parser.add_argument("-t", "--threshold", help = "Find genomes that are within this distance", default=.02, type=float)
    parser.add_argument("-x", "--mash_exe", help = "Path to MASH exectuable", default="/common/contrib/bin/mash-Linux64-v1.1.1/mash")
    args = parser.parse_args()

    query_genome(args.sketch_files, args.referenc_genome, args.threshold, args.mash_exe)

main()
