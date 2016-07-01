#!/usr/bin/env python

import os, argparse
import pandas as pd
from Bio import SeqIO

def mash(fasta_dir):
    mash="/home/asanchez/local/bin/mash-Linux64-v1.1/mash"
    sketch_file = os.path.join(fasta_dir, "all.msh")
    all_fastas = os.path.join(fasta_dir, "*.fasta")
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.tsv")
    sketch_command = "{} sketch -o {} {}".format(mash, sketch_file, all_fastas)
    distance_command = "{} dist -p 4 -t {} {} > {}".format(mash, sketch_file, all_fastas, distance_matrix)
    os.system(sketch_command)
    os.system(distance_command)

def clean_up_matrix(distance_matrix):
    distance_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    new_index = []
    for i in distance_df.index:
        name = i.split("/")[-1]
        name = name.split("_")[0:2]
        new_index.append("".join(name))

    distance_df.index = new_index
    distance_df.columns = new_index
    distance_df.to_csv(distance_matrix, sep="\t")

def generate_fasta_stats(fasta_dir):
    for root, dirs, files, in os.walk(fasta_dir):
        file_names = []
        file_sizes = []
        contig_totals = []
        lengths = []
        N_counts = []

        for f in files:
            if f.endswith(".fasta"):
                file_names.append(f)
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

                # Read the N count for each contig into a list
                # Append the total N count to N_counts
                N_count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
                N_counts.append(sum(N_count))

        SeqDataSet = list(zip(file_sizes, contig_totals, lengths, N_counts))
        seq_df = pd.DataFrame(data = SeqDataSet, index=file_names, columns=["File Size", "Contigs", "Total Length", "N Count"])

        new_index = []
        for i in seq_df.index:
            name = i.split("_")[0:2]
            new_index.append("".join(name))

        seq_df.index = new_index
        seq_df.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")

def concat_distances_to_stats(fasta_dir):
    """Incorporate the average distances from mash results into stats.csv"""
    stats_file = os.path.join(fasta_dir, "stats.csv")
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.tsv")
    stats_df = pd.read_csv(stats_file, index_col=0)
    distance_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    average_distances = distance_df.mean()
    stats_df["Average Distances"] = average_distances
    stats_df.to_csv(os.path.join(fasta_dir, "stats.csv"))

def quality_control_filters(fasta_dir, max_n_count, max_contigs):
    stats_file = os.path.join(fasta_dir, "stats.csv")
    stats_df = pd.read_csv(stats_file, index_col=0)
    quantiles_stats = stats_df.quantile([.1, .9])
    quantiles_stats.index = [10, 90]
    lower_percentile_stats = quantiles_stats.loc[10]
    upper_percentile_stats = quantiles_stats.loc[90]

    print("Fasta stats:  ", end="\n\n")
    print(stats_df, end="\n\n")

    print("Percentiles:  ", end="\n\n")
    print(quantiles_stats, end="\n\n")

    filtered_df = stats_df[stats_df["N Count"] <= max_n_count ]
    filtered_df = stats_df[stats_df["Contigs"] <= max_contigs ]
    filtered_df = stats_df[(stats_df["Total Length"] >= lower_percentile_stats["Total Length"]) &
            (stats_df["Total Length"] <= upper_percentile_stats["Total Length"])]
    filtered_df = stats_df[(stats_df["Average Distances"] >= lower_percentile_stats["Average Distances"]) &
            (stats_df["Average Distances"] <= upper_percentile_stats["Average Distances"])]
    print("Filtered stats:  ", end="\n\n")
    print(filtered_df, end="\n\n")
    return filtered_df

def record_pass_or_fail(fasta_dir, filtered_df, max_n_count, max_contigs):
    stats_file = os.path.join(fasta_dir, "stats.csv")
    stats_df = pd.read_csv(stats_file, index_col=0)
    quantiles_stats = stats_df.quantile([.1, .9])
    quantiles_stats.index = [10, 90]
    lower_percentile_stats = quantiles_stats.loc[10]
    upper_percentile_stats = quantiles_stats.loc[90]
    n_count_bool = (stats_df.iloc[:]["N Count"]) < max_n_count
    contigs_bool = (stats_df.iloc[:]["Contigs"]) < max_contigs
    distances_bool = (stats_df.iloc[:]["Average Distances"] >= lower_percentile_stats["Average Distances"]) \
                   & (stats_df.iloc[:]["Average Distances"] <= upper_percentile_stats["Average Distances"])
    lengths_bool = (stats_df.iloc[:]["Total Length"] >= lower_percentile_stats["Total Length"]) \
                   & (stats_df.iloc[:]["Total Length"] <= upper_percentile_stats["Total Length"])
    bool_df = pd.concat([lengths_bool, n_count_bool, contigs_bool, distances_bool], axis=1)
    print(bool_df)

    passed_files = os.path.join(fasta_dir, "passed.txt")

    try:
        os.remove(passed_files)
    except OSError:
        pass

    for accession in filtered_df.index:
        with open(passed_files, "a") as file:
            fasta = os.path.join(fasta_dir, accession)
            file.write(fasta)
            file.write("\n")

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=2000)
    parser.add_argument("-c", "--max_contigs", help = "Maximum number of contigs acceptable", type=int, default=500)
    parser.add_argument("-m", "--mash", help = "Create a sketch file and distance matrix", action="store_true")
    args = parser.parse_args()

    fasta_dir = args.fasta_dir

    if args.mash:
        mash(fasta_dir)

    distance_matrix = os.path.join(fasta_dir, "distance_matrix.tsv")
    clean_up_matrix(distance_matrix)
    generate_fasta_stats(fasta_dir)
    concat_distances_to_stats(fasta_dir)
    filtered_df = quality_control_filters(fasta_dir, args.max_n_count, args.max_contigs)
    record_pass_or_fail(fasta_dir, filtered_df, args.max_n_count, args.max_contigs)

Main()
