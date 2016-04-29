#!/usr/bin/env python3

import argparse
import pandas as pd

def Main():
    parser = argparse.ArgumentParser(description = 'Get information for a given accession number.  '\
                                    'Columns 6 and 7 are displayed.')
    parser.add_argument('accession_number', help='The accession number you want to know about.')
    parser.add_argument('-w', '--wget_assembly_summary', help='Download assembly_summary.txt')
    args = parser.parse_args()

    if args.wget_assembly_summary:
        import os, subprocess
        wget = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'

        try:
            os.remove('assembly_summary.txt')
            subprocess.call(['wget', assembly_summary])
        except OSError:
            subprocess.call(['wget', assembly_summary])

        df = pd.read_csv('assembly_summary.txt', delimiter = '\t', index_col=0)
        print(df.loc[args.accession_number][6:7][0])
        print(df.loc[args.accession_number][7:8][0])

    else:
        df = pd.read_csv('assembly_summary.txt', delimiter = '\t', index_col=0)
        print(df.loc[args.accession_number][6:7][0])
        print(df.loc[args.accession_number][7:8][0])

if __name__ == "__main__":
    Main()
