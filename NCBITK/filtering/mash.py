#!/usr/bin/env python

import os, argparse
import glob

def write_sketch_commands(genbank_mirror, assembly_summary, new_genomes):

    sketch_commands = os.path.join(genbank_mirror, ".info", "sketch_commands.sh")
    if os.path.isfile(sketch_commands):
        os.remove(sketch_commands)

    with open(sketch_commands, 'a') as cmds:
        for genome in new_genomes:
            species_dir = assembly_summary.scientific_name.loc[genome]
            fasta = os.path.join(genbank_mirror, species_dir, "{}*fasta".format(genome))
            fasta = glob.glob(fasta)

            try:
                fasta = fasta[0]
            except IndexError:
                continue

            sketch_dst = os.path.join(genbank_mirror, species_dir, "{}.msh".format(genome))
            cmd = "/common/contrib/bin/mash-Linux64-v1.1.1/mash sketch {} -o {}\n".format(fasta, sketch_dst)
            cmds.write(cmd)

    return sketch_commands # get the line count of this file

def main():
    parser = argparse.ArgumentParser(description = "Run MASH on entire genbank collection.")
    parser.add_argument("genbank_mirror", help = "directory containing your FASTA files")
    args = parser.parse_args()
    get_fastas(args.genbank_mirror)

if __name__ == "__main__":
    main()
