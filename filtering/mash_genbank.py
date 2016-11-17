#!/usr/bin/env python

import os, argparse
from glob import glob

def get_fastas(genbank_mirror):

    sketch_ids=[]
    genome_ids=[]
    genome_paths=[]

    for root, dirs, files, in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("msh"):
                sketch_id = "_".join(f.split(".")[1:3])
                sketch_ids.append(sketch_id)
            elif f.endswith("fasta"):
                genome_path = os.path.join(root, f)
                genome_paths.append(genome_path)
                genomes.append((genome_id, genome_path))

    sketch_commands = os.path.join(genbank_mirror, ".info", "sketch_commands.sh")
    if os.path.isfile(sketch_commands):
        os.remove(sketch_commands)

    with open(sketch_commands) as cmds:
        for genome in genome_paths:
            genome_id = genome.split("/")[-1]
            genome_id = "_".join(genome.split("_")[:2])
            genome_ids.append(genome_id)
            species_dir = "/".join(genome.split("/")[:-1])
            if genome_id not in sketch_ids:
                sketch_dst = os.path.join(species_dir, ".{}.msh".format(genome_id))
                cmd = "/common/contrib/bin/mash-Linux64-v1.1.1/mash sketch {} -o {}\n".format(genome, sketch_dst)
                cmds.write(cmd)

    for sketch_id in sketch_ids:
        if sketch_id not in genome_ids:
            sketch_path = glob(os.path.join(genbank_mirror, "*", ".{}.msh".format(sketch_id)))
            os.remove(sketch_path[0])

    return sketch_commands # get the line count of this file

def run_mash_array(sketch_commands):

def main():
    parser = argparse.ArgumentParser(description = "Run MASH on entire genbank collection.")
    parser.add_argument("genbank_mirror", help = "directory containing your FASTA files")
    args = parser.parse_args()
    get_fastas(args.genbank_mirror)

main()
