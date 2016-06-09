#!/usr/bin/env python

import pandas as pd
from os import path
import argparse

def read_csv(organism):
    df = pd.read_csv(path.join(organism, "stats.csv"))
    return df

def describe(df):
    return df.describe()

def N_count(df, N_count=200):
    N_count = df[df["N Count"] < N_count]
    return N_count

def contigs(df, contigs=200):
    contigs = df[df["Contigs"] < contigs]
    return contigs

def df_filters(df, fsize=5000000, contigs=500, N_count=500):
    df = df[ (df["File Size"] >= fsize) & (df["Contigs"] <= contigs) & (df["N_count"] <= N_count) ]
    return df

df mash(reference):
    # subprocess(path to mash executable)
    # make sketch, etc.

def Main():
    parser = argparse.ArgumentParser(description = "Functions to explore qc stats")
    parser.add_argument("organism", help="Path to the fasta dir you want to know about")
    parser.add_argument("-d", "--describe", help="Summary stats.csv", action='store_true')
    args = parser.parse_args()

    df = read_csv(args.organism)
    description = describe(df)
    std_dict = {i:description[i]["std"] for i in description.columns}
    mean_dict = {i:description[i]["mean"] for i in description.columns}
    min_dict = {i:description[i]["min"] for i in description.columns}
    max_dict = {i:description[i]["max"] for i in description.columns}

    if args.describe:
        describe(df)

    elif args.N_count:
        N_count(df)

    else:
        print("wtf")

Main()
