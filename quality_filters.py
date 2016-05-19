#!/usr/bin/env python

import os
import pandas as pd
from Bio import SeqIO

for root, dirs, files, in os.walk("Escherichia_coli"):
    file_names = []
    file_sizes = []
    seq_totals = []
    seq_ids = []
    seq_lengths = []
    percentage_N = []

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
print(seq_df.head(20))

