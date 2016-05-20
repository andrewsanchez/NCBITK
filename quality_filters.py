#!/usr/bin/env python

import os, argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

file_names = []
file_sizes = []
seq_totals = []
seq_ids = []
seq_lengths = []
percentage_N = []

def quality_control(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):
        for f in files:

            """ information for whole files
            number of sequences, file size """

            fasta = (os.path.join(root, f))
            accession = f.split("_")[0:2]
            accession = "_".join(accession)
            seqs = [seq_record.seq for seq_record in SeqIO.parse(fasta, "fasta")]

            """ length of individual seqs; percentage Ns for individual seqs """

            for seq_record in SeqIO.parse(fasta, "fasta"):
                seq = str(seq_record.seq).upper()
                percent_N = float(seq.count("N"))/float(len(seq))*100

                file_names.append(accession)
                file_sizes.append(os.path.getsize(fasta))
                seq_totals.append(len(seqs))
                seq_ids.append(seq_record.id)
                seq_lengths.append(len(seq_record.seq))
                percentage_N.append(percent_N)

        SeqDataSet = list(zip(file_names, file_sizes, seq_totals, seq_ids, seq_lengths, percentage_N))
        seq_df = pd.DataFrame(data = SeqDataSet, columns=["File Name", "Size", "Total Seqs", "ID", "Length", "N's"])
        seq_df.to_csv(os.path.join(root, "stats.csv"), index=False)

def plot_data():

    fig, (ax0, ax1) = plt.subplots(ncols=2)
    ax0.hist(seq_totals)
    ax0.set_title("Contigs")
    ax1.hist(file_sizes)
    ax1.set_title("file_sizes")
    plt.show()

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("--percent_N", help = "percent of N's")
    args = parser.parse_args()
    quality_control(args.fasta_dir)
    plot_data()

Main()
