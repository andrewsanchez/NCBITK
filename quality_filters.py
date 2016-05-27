#!/usr/bin/env python

import os, argparse
import pandas as pd
from Bio import SeqIO
#import matplotlib.pyplot as plt
#import numpy as np

def quality_control(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):

        file_names = []
        file_sizes = []
        contig_totals = []
        lengths = []
        N_percent_totals = []

        for f in files:

            if "fasta" in f:
                fasta = (os.path.join(root, f))
                file_names.append(f)
                file_sizes.append(os.path.getsize(fasta))

                # Read all contigs for current fasta into list
                # Append the total number of contigs to contig_totals
                contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                # Append the sum of all contig lengths to lengths
                contig_lengths = [ len(str(seq)) for seq in contigs ]
                total_length  = sum(contig_lengths)
                lengths.append(total_length)

                # Read the N count for each contig into al ist
                # Append the total percent N to N_percent_totals
                N_count = [ float(str(seq.upper()).count("N")) for seq in contigs ]
                N_percent = (sum(N_count)/float(total_length))*100
                N_percent_totals.append(N_percent)

            SeqDataSet = list(zip(file_names, file_sizes, contig_totals, lengths, N_percent_totals))
            seq_df = pd.DataFrame(data = SeqDataSet, columns=["Accession", "File Size", "Contigs", "Total Length", "% N"])
            seq_df.to_csv(os.path.join(root, "stats.csv"), index=False)
            with open(os.path.join(root, "stats_summary.txt"), "w") as stats:
                stats.write(str(seq_df.describe()))

def plot_data():
    fig, (ax0, ax1) = plt.subplots(ncols=2)
    ax0.hist(seq_totals)
    ax0.set_title("Contigs")
    ax1.hist(file_sizes)
    ax1.set_title("file_sizes")
    plt.savefig("stats.png")

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("--percent_N", help = "percent of N's")
    args = parser.parse_args()
    quality_control(args.fasta_dir)
    #plot_data()

Main()
