#!/usr/bin/env python

import os, argparse

def get_fastas(genbank):
    with open("/scratch/aas229/jobs/arrays/mash_genbank_commands.sh", "a") as cmds:
        for root, dirs, files, in os.walk(genbank):
            for f in files:
                name = f.split("/")[-1]
                if f.endswith("fasta"):
                    cmds.write("/common/contrib/bin/mash-Linux64-v1.1.1/mash sketch {} -o sketch_files/{}\n".format(os.path.join(root,f), name))

def main():
    parser = argparse.ArgumentParser(description = "Run MASH on entire genbank collection.")
    parser.add_argument("genbank", help = "directory containing your FASTA files")
    args = parser.parse_args()
    get_fastas(args.genbank)

main()
