#!/usr/bin/env python3

import pandas as pd

def findaccession(accession, wget='nowget'):
    if wget == 'wget':
        import os, subprocess

        if os.path.isfile('assembly_summary.txt'):
           os.remove('assembly_summary.txt')
           subprocess.call(['wget',
                           'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])
        else:
           subprocess.call(['wget',
                           'ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'])

        df = pd.read_csv('assembly_summary.txt', delimiter = '\t', index_col=0)
        print(df.loc[accession][6:7][0])
        print(df.loc[accession][7:8][0])

    else:
        df = pd.read_csv('assembly_summary.txt', delimiter = '\t', index_col=0)
        print(df.loc[accession][6:7][0])
        print(df.loc[accession][7:8][0])

# if __name__ == "__main__":
#   import sys
#   findaccession(sys.argv[1], sys.argv[2])
