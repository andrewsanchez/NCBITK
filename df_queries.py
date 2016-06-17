#!/usr/bin/env python

import pandas as pd
from os import path, getcwd, mkdir
from subprocess import call
import argparse

# make dirs to hold reorganize fastas into passed/failed groups
def mkdirs(organism):
    pass_dir = path.join(organism, "pass")
    fail_dir = path.join(organism, "fail")

    if not path.isdir(pass_dir): 
        mkdir(pass_dir)
    if not path.isdir(fail_dir):
        mkdir(fail_dir)

def read_csv(organism):
    df = pd.read_csv(path.join(organism, "stats.csv"))
    print("Number of fastas in this directory:  {}".format(len(df)))
    return df

def describe(df):
    print("Summary of fastas prior to filtering:")
    description = df.describe()
    print(description)
    return description

def N_count(df, N_count=200):
    N_count = df[df["N Count"] < N_count]
    return N_count

def contigs(df, contigs=200):
    contigs = df[df["Contigs"] < contigs]
    return contigs

def df_filters(df, fsize=5000000, contigs=500, N_count=500):
    description = df.describe()
    std_dict = {i:description[i]["std"] for i in description.columns}
    mean_dict = {i:description[i]["mean"] for i in description.columns}
    min_dict = {i:description[i]["min"] for i in description.columns}
    max_dict = {i:description[i]["max"] for i in description.columns}

    fsize=(int(mean_dict["File Size"]))-(int(std_dict["File Size"]*2))
    contigs=200
    N_count=200

    print("Enter the max values for number of N's, contigs, and file size.")
    print("Press enter to use the default values file_size={}, contigs={}, N_count={}".format(fsize, contigs, N_count))
    fsize = int(input("  > Max File size:  ") or fsize)
    contigs = int(input("  > Max contigs:  ") or contigs)
    N_count = int(input("  > Max N count:  ") or N_count)
    df = df[ (df["File Size"] >= fsize) & (df["Contigs"] <= contigs) & (df["N Count"] <= N_count) ]
    print("There are {} fastas passing the filter".format(len(df)))
    again = input("Try again with different values?\n  > y/n:  ")
    if again == "y":
        dialogue(df)
    else:
        pass
    return df

def dialogue(df):
    print("Enter the max values for number of N's, contigs, and file size.")
    print("Press enter to use the default values file_size=5000000, contigs=500, N_count=500")
    fsize = int(input("  > Max File size:  ") or "5000000")
    contigs = int(input("  > Max contigs:  ") or "500")
    N_count = int(input("  > Max N count:  ") or "500")
    df = df_filters(df, fsize=fsize, contigs=contigs, N_count=N_count)
    print("There are {} fastas passing the filter".format(len(df)))
    again = input("Try again with different values?\n  > y/n:  ")
    if again == "y":
        dialogue(df)
    else:
        pass
    return df

#def mash(reference):
    # subprocess(path to mash executable)
    # make sketch, etc.

def sort_files(df, organism):
    organize = input("Organize files based on these parameters?\n  > y/n:  ")
    if organize == "y":
        for f in df["Accession"]:
            source = path.join(organism, f)
            dest = path.join(organism, "pass", f)
            call([ "ln", "-s", source, dest ])
        print("Fastas that passed the quality control filters are now linked in the folder 'pass'")
        print("Fastas that failed are linked in the folder 'fail'")
    elif organize == "n":
        print("Fastas were not reorganized.  Ciao!")
    else:
        print("Please answer y or n")


def Main():
    parser = argparse.ArgumentParser(description = "Functions to explore qc stats")
    parser.add_argument("organism", help="Path to the fasta dir you want to know about")
    parser.add_argument("-d", "--describe", help="Summary stats.csv", action='store_true')
    parser.add_argument("-n", "--N_count", default=500, type=int)
    parser.add_argument("-c", "--contigs", default=500, type=int)
    parser.add_argument("-f", "--file_size", default=5000000, type=int)
    args = parser.parse_args()

    organism = args.organism
    mkdirs(organism)
    df = read_csv(organism)
    description = describe(df)
    std_dict = {i:description[i]["std"] for i in description.columns}
    mean_dict = {i:description[i]["mean"] for i in description.columns}
    min_dict = {i:description[i]["min"] for i in description.columns}
    max_dict = {i:description[i]["max"] for i in description.columns}
    df = df_filters(df)
    sort_files(df, organism)

Main()
