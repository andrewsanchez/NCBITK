#!/usr/bin/env python

import os, argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
#import matplotlib.pyplot as plt

def quality_control(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):

        file_names = []
        file_sizes = []
        contig_totals = []
        lengths = []
        N_percent_totals = []

        for f in files:
            fasta = (os.path.join(root, f))
            file_names.append(f)
            file_sizes.append(os.path.getsize(fasta))

            # Read all contigs for this fasta into a list
            contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]

            # Append the total number of contigs to contig_totals
            contig_totals.append(len(contigs))

            # Create a list to hold the length of each contig
            # Append the total to lengths
            contig_lengths = [ len(str(seq)) for seq in contigs ]
            lengths.append(sum(contig_lengths)

            # Create a list to hold the N count for each contig
            # Append the total percent N to N_percent_totals
            N_count = [ float(str(seq)).count("N") for seq in contigs ]
            N_percent = sum(N_count)/float(total_length)
            N_percent_totals.append(N_percent)

           #contigs = [seq.seq for seq in SeqIO.parse(fasta, "fasta")]
           #total_contigs.append(len(contigs))
           #contig_lengths = [len(str(seq.seq)) for seq in SeqIO.parse(fasta, "fasta")]
           #contig_lengths_totals.append(sum(contig_lengths))
           #total_Ns = [float(str(seq.seq).count("N")) for seq in SeqIO.parse(fasta, "fasta")]
           #percent_N_totals.append("%.1f" % sum(total_Ns)/sum(contig_lengths))

        SeqDataSet = list(zip(file_names, file_sizes, contig_totals, lengths, N_percent_totals))
        seq_df = pd.DataFrame(data = SeqDataSet, columns=["Accession", "File Size", "Contigs", "Total Length", "% N"])
        seq_df.to_csv(os.path.join(root, "stats.csv"), index=False)

#def save_stats():

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
