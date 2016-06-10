#!/usr/bin/env python

import pandas as pd
from os import path, getcwd, mkdir
from subprocess import call
import argparse

# make dirs to hold reorganize fastas into passed/failed groups
def mkdirs():
    cwd = getcwd()
    pass_dir = path.join(cwd, "pass")
    fail_dir = path.join(cwd, "fail")

    if not path.isdir(pass_dir): 
        mkdir(pass_dir)
    if not path.isdir(fail_dir):
        mkdir(fail_dir)

def read_csv(organism):
    df = pd.read_csv(path.join(organism, "stats.csv"))
    return df

def describe(df):
    print("Summary of unfiltered fastas:")
    description = df.describe()
    print(description)

def N_count(df, N_count=200):
    N_count = df[df["N Count"] < N_count]
    return N_count

def contigs(df, contigs=200):
    contigs = df[df["Contigs"] < contigs]
    return contigs

def df_filters(df, fsize=5000000, contigs=500, N_count=500):
    df = df[ (df["File Size"] >= fsize) & (df["Contigs"] <= contigs) & (df["N Count"] <= N_count) ]
    return df

def dialogue():
    print("There are {} files matching the filter".format(len(df)))
    print("Summary of the unfiltered data frame:")
    print(description)
    print("Please enter the max values for N, contigs, and file size.")
    print("Press enter at each prompt to use default values:\n\
    N_count=500\n\
    contigs=500\n\
    file_size=5000000\n")
    N_count = int(input("  > Max N count:  "))
    contigs = int(input("  > Max N count:  "))
    File_size = int(input("  > Max N count:  "))

#def mash(reference):
    # subprocess(path to mash executable)
    # make sketch, etc.

#def sort_files(df, organism):
#   for f in df["Accession"]:
#       source = join(organism, f)
#       dest = join(organism, "pass", f)
#       call.([ "ln", "-s", source, dest ])


def Main():
    parser = argparse.ArgumentParser(description = "Functions to explore qc stats")
    parser.add_argument("organism", help="Path to the fasta dir you want to know about")
    parser.add_argument("-d", "--describe", help="Summary stats.csv", action='store_true')
    parser.add_argument("-n", "--N_count", default=500, type=int)
    parser.add_argument("-c", "--contigs", default=500, type=int)
    parser.add_argument("-f", "--file_size", default=5000000, type=int)
    args = parser.parse_args()

    organism = args.organism
    df = read_csv(organism)
    describe(df)
    description = describe(df)

    std_dict = {i:description[i]["std"] for i in description.columns}
    mean_dict = {i:description[i]["mean"] for i in description.columns}
    min_dict = {i:description[i]["min"] for i in description.columns}
    max_dict = {i:description[i]["max"] for i in description.columns}

    """

    

    """

    df = df_filters(df, fsize=args.file_size, contigs=args.contigs, N_count=args.N_count)
    print("Files matching this criteria:  {}".format(len(df)))

Main()
