#!/usr/bin/env python

import os, argparse
from subprocess import Popen
import pandas as pd
from Bio import SeqIO

def mash(fasta_dir):
    mash="/home/asanchez/local/bin/mash-Linux64-v1.1/mash"
    sketch_file = os.path.join(fasta_dir, "all.msh")
    all_fastas = os.path.join(fasta_dir, "*.fasta")
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.csv")
    sketch_command = "{} sketch -o {} {}".format(mash, sketch_file, all_fastas)
    distance_command = "{} dist -p 4 -t {} {} > {}".format(mash, sketch_file, all_fastas, distance_matrix)
    Popen(sketch_command, shell="True").wait()
    Popen(distance_command, shell="True").wait()
    clean_up_matrix(distance_matrix)

def clean_up_matrix(distance_matrix):
    distance_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    new_index = []
    for i in distance_df.index:
        name = i.split("/")[-1]
        name = name.split("_")[0:2]
        new_index.append("_".join(name))

    distance_df.index = new_index
    distance_df.columns = new_index
    distance_df.to_csv(distance_matrix, sep="\t")

def generate_fasta_stats(fasta_dir, distance_matrix):
    for root, dirs, files, in os.walk(fasta_dir):
        accessions = []
        file_sizes = []
        contig_totals = []
        lengths = []
        N_Counts = []

        for f in files:
            if f.endswith(".fasta"):
                accession = f.split("_")[0:2]
                accessions.append("_".join(accession))
                fasta = (os.path.join(root, f))
                file_sizes.append(os.path.getsize(fasta))

                # Read all contigs for current fasta into list
                # Append the total number of contigs to contig_totals
                contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                # Append the sum of all contig lengths to lengths
                contig_lengths = [ len(str(seq)) for seq in contigs ]
                lengths.append(sum(contig_lengths))

                # Read the N_Count for each contig into a list
                # Append the total N_Count to N_Counts
                N_Count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
                N_Counts.append(sum(N_Count))

        SeqDataSet = list(zip(file_sizes, contig_totals, lengths, N_Counts))
        stats_df = pd.DataFrame(data = SeqDataSet, index=accessions, columns=["File_Size", "Contigs", "Total_Length", "N_Count"])
        distance_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
        average_distances = distance_df.mean()
        stats_df["Avg_Distances"] = distance_df.mean()
        stats_df.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")
        print(stats_df)
        print("\n\n")

        return stats_df

def assess_stats_df(fasta_dir, stats_df, max_n_count, max_contigs):
    quantiles = stats_df.quantile([.1, .9])
    quantiles.index = [10, 90]
    lower_percentile = quantiles.loc[10]
    upper_percentile = quantiles.loc[90]
    print("Percentiles:  ", end="\n\n")
    print(quantiles)

    passed_df = stats_df[stats_df["N_Count"] <= max_n_count ]
    passed_df = passed_df[passed_df["Contigs"] <= max_contigs ]

    passed_df = passed_df[(passed_df["Total_Length"] >= lower_percentile["Total_Length"]) &
            (passed_df["Total_Length"] <= upper_percentile["Total_Length"])]

    passed_df = passed_df[(passed_df["Avg_Distances"] >= lower_percentile["Avg_Distances"]) &
            (passed_df["Avg_Distances"] <= upper_percentile["Avg_Distances"])]

    passed_files = os.path.join(fasta_dir, "passed.txt")
    if os.path.isfile(passed_files):
        os.remove(passed_files)
    for accession in passed_df.index:
        with open(passed_files, "a") as file:
            file.write(accession+"\n")

def bool_df_for_failed(fasta_dir, stats_df, max_n_count, max_contigs):
    quantiles = stats_df.quantile([.1, .9])
    quantiles.index = [10, 90]
    lower_percentile = quantiles.loc[10]
    upper_percentile = quantiles.loc[90]

    n_count_bool = (stats_df.iloc[:]["N_Count"]) <= max_n_count
    contigs_bool = (stats_df.iloc[:]["Contigs"]) <= max_contigs

    distances_bool = (stats_df.iloc[:]["Avg_Distances"] >= lower_percentile["Avg_Distances"]) \
                   & (stats_df.iloc[:]["Avg_Distances"] <= upper_percentile["Avg_Distances"])

    lengths_bool = (stats_df.iloc[:]["Total_Length"] >= lower_percentile["Total_Length"]) \
                   & (stats_df.iloc[:]["Total_Length"] <= upper_percentile["Total_Length"])

    for i in stats_df.index:
        N_Count = stats_df.loc[i, "N_Count"]
        Contigs = stats_df.loc[i, "Contigs"]
        Avg_Distances = stats_df.loc[i, "Avg_Distances"]
        print("avg_distances")
        print(Avg_Distances)
        Total_Length = stats_df.loc[i, "Total_Length"]

        if N_Count > max_n_count:
            stats_df.loc[i, "N_Count"] = "*{}*".format(str(stats_df.loc[i, "N_Count"]))
        if Contigs > max_contigs:
            stats_df.loc[i, "Contigs"] = "*{}*".format(str(stats_df.loc[i, "Contigs"]))
        if Avg_Distances <= lower_percentile["Avg_Distances"] or Avg_Distances >= upper_percentile["Avg_Distances"]:
            stats_df.loc[i, "Avg_Distances"] = "*{}*".format(str(stats_df.loc[i, "Avg_Distances"]))
        if Total_Length <= lower_percentile["Total_Length"] or Total_Length >= upper_percentile["Total_Length"]:
            stats_df.loc[i, "Total_Length"] = "*{}*".format(str(stats_df.loc[i, "Total_Length"]))

    bool_df = pd.concat([lengths_bool, n_count_bool, contigs_bool, distances_bool], axis=1)
    bool_df.to_csv(os.path.join(fasta_dir, "failed_tf.csv"))
    stats_df.to_csv(os.path.join(fasta_dir, "failed.csv"))
    print("bool_df", end="\n\n")
    print(bool_df)
    print("stats_df", end="\n\n")
    print(stats_df)

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=2000)
    parser.add_argument("-c", "--max_contigs", help = "Maximum number of contigs acceptable", type=int, default=500)
    parser.add_argument("-m", "--mash", help = "Create a sketch file and distance matrix", action="store_true")
    args = parser.parse_args()

    fasta_dir = args.fasta_dir
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.csv")

    if args.mash:
        mash(fasta_dir)

    stats_df = generate_fasta_stats(fasta_dir, distance_matrix)
    failed_df = assess_stats_df(fasta_dir, stats_df, args.max_n_count, args.max_contigs)
    print(failed_df)
    bool_df_for_failed(fasta_dir, stats_df, args.max_n_count, args.max_contigs)

Main()
