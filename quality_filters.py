#!/usr/bin/env python

import os, argparse
import pandas as pd
from Bio import SeqIO

def quality_control(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):

        file_names = []
        file_sizes = []
        contig_totals = []
        lengths = []
        N_counts = []

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
                lengths.append(sum(contig_lengths))

                # Read the N count for each contig into a list
                # Append the total N count to N_counts
                N_count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
                N_counts.append(sum(N_count))

            SeqDataSet = list(zip(file_names, file_sizes, contig_totals, lengths, N_counts))
            seq_df = pd.DataFrame(data = SeqDataSet, columns=["Accession", "File Size", "Contigs", "Total Length", "N Count"])
            seq_df.to_csv(os.path.join(root, "stats.csv"), index=False)

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("--percent_N", help = "percent of N's")
    args = parser.parse_args()
    quality_control(args.fasta_dir)

Main()
