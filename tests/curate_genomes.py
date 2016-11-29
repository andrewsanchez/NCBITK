#!/usr/bin/env python

import os, argparse
from random import shuffle

def mixer(test_species, near_neighbors, target):
    dirty_collection = test_species.strip("/").split("/")[-1]
    dirty_collection = "{}_dirty".format(dirty_collection)
    print(dirty_collection)
    target_dir = os.path.join(target, dirty_collection)
    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)
    if os.path.isfile(os.path.join(target_dir, "near_neighbors.txt")):
        os.remove(os.path.join(target_dir, "near_neighbors.txt"))

    genomes = [f for f in os.listdir(test_species)]
    near_neighbor_genomes = [f for f in os.listdir(near_neighbors)]
    shuffle(near_neighbor_genomes)
    fraction_of_test_species = int(len(genomes) * .015)

    for name in genomes:
        src = os.path.join(test_species, name)
        dst = os.path.join(target_dir, name)
        os.link(src, dst)

    with open(os.path.join(target_dir, "near_neighbors.txt"), "a") as f:
        for name in near_neighbor_genomes[:fraction_of_test_species]:
            f.write(name + "\n")
            src = os.path.join(near_neighbors, name)
            dst = os.path.join(target_dir, name)
            os.link(src, dst)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("test_species", help = "near_neighbors genbank")
    parser.add_argument("near_neighbors", help = "near_neighbors genbank")
    parser.add_argument("target", help = "where to link mixed genomes")
    args = parser.parse_args()

    mixer(args.test_species, args.near_neighbors, args.target)

main()
