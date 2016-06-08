#!/usr/bin/env python

import pandas as pd
from os import path
import argparse

def read_csv(organism):
    df = pd.read_csv(path.join(organism, "stats.csv"))
    return df


def describe(df):
    print(df.describe())

def Main():
    parser = argparse.ArgumentParser(description = "Functions to explore qc stats")
    parser.add_argument("organism", help="Path to the fasta dir you want to know about")
    parser.add_argument("-d", "--describe", help="Summary stats.csv")
    args = parser.parse_args()

    read_csv(args.organism)

    if args.describe:
        describe(df)
    else:
        print("wtf")

Main()
