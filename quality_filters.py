#!/usr/bin/env python

import os, argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
#import matplotlib.pyplot as plt

file_names = []
file_sizes = []
total_contigs = []
contig_lengths_totals = []
percent_N_totals = []

def quality_control(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):
        for f in files:

            """ information for whole files
            number of sequences, file size """

            fasta = (os.path.join(root, f))
            accession = f.split("_")[0:2]
            accession = "_".join(accession)
            file_names.append(accession)
            file_sizes.append(os.path.getsize(fasta))
            contigs = [seq.seq for seq in SeqIO.parse(fasta, "fasta")]
            total_contigs.append(len(contigs))
            contig_lengths = [len(str(seq.seq)) for seq in SeqIO.parse(fasta, "fasta")]
            contig_lengths_totals.append(sum(contig_lengths))
            percent_Ns = [float(str(seq.seq).count("N"))/float(len(str(seq.seq).upper()))*100 for seq in SeqIO.parse(fasta, "fasta")]
            percent_N_totals.append("%.2f" % sum(percent_Ns))

            """ length of individual seqs; percentage Ns for individual seqs """

        SeqDataSet = list(zip(file_names, file_sizes, total_contigs, contig_lengths_totals, percent_N_totals))
        seq_df = pd.DataFrame(data = SeqDataSet, columns=["Accession", "File Size", "Contigs", "Length", "% N"])
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
