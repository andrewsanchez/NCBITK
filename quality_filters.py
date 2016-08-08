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
    clean_up_matrix_and_get_quantiles(distance_matrix)

def clean_up_matrix_and_get_quantiles(distance_matrix):
    distances_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
  # distances_pass_fail = os.path.join(fasta_dir, "distance_matrix_pass_fail.csv")
    new_index = []
    for i in distances_df.index:
        name = i.split("/")[-1]
        name = name.split("_")[0:2]
        new_index.append("_".join(name))

    distances_df.index = new_index
    distances_df.columns = new_index
    distances_df.to_csv(distance_matrix, sep="\t")

  # distances_upper_percentile = distances_df.quantile(upper_percentile)
  # pass_fail_distances = pd.DataFrame(index=distances_df.index)

  # for accession in distances_df.columns:
  #     for distance in distances_df.iloc[:][accession]:
  #         if distance > distances_upper_percentile.loc[accession]:
  #             pass_fail_distances[accession] = "F"
  #         else:
  #             pass_fail_distances[accession] = "P"

  # distances_df.to_csv(distances_pass_fail, sep="\t")

def generate_fasta_stats(fasta_dir, distance_matrix):

    distances_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    avg_distances = ["{:05.4f}".format(i) for i in distances_df.mean()]

    for root, dirs, files, in os.walk(fasta_dir):
        accessions = []
        contig_totals = []
        lengths = []
        N_Counts = []

        for f in files:
            if f.endswith(".fasta"):
                print("Getting stats for {}".format(f))
                accession = f.split("_")[0:2]
                accessions.append("_".join(accession))
                fasta = (os.path.join(root, f))

                # Read all contigs for current fasta into list
                # Append the total number of contigs to contig_totals
                try:
                    contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                except UnicodeDecodeError:
                    print("{} threw UnicodeDecodeError".format(f))
                    with open(os.path.join(fasta_dir, "quality_filters_log.txt"), "a") as log:
                        log.write("{} threw UnicodeDecodeError".format(f))

                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                # Append the sum of all contig lengths to lengths
                contig_lengths = [ len(str(seq)) for seq in contigs ]
                lengths.append(sum(contig_lengths))

                # Read the N_Count for each contig into a list
                # Append the total N_Count to N_Counts
                N_Count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
                N_Counts.append(sum(N_Count))

        SeqDataSet = list(zip(contig_totals, lengths, N_Counts, avg_distances))
        stats_df = pd.DataFrame(data = SeqDataSet, index=accessions, columns=["Contigs", "Total_Length", "N_Count", "Avg_Distances"], dtype="float64")
        stats_df.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")

        return stats_df

def assess_stats_df(fasta_dir, stats_df, max_n_count, max_contigs, lower_percentile, upper_percentile):
    print(stats_df.quantile([lower_percentile, upper_percentile]))
    quantiles_stats = stats_df.quantile([lower_percentile, upper_percentile])
    print("Quantiles:", end="\n\n")
    print(quantiles_stats)
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]
    print("Percentiles:  ", end="\n\n")
    print(quantiles_stats)

    passed_df = stats_df[stats_df["N_Count"] <= max_n_count ]
    print("passed_df\n\n")
    print(passed_df)
    passed_df = passed_df[passed_df["Contigs"] <= max_contigs ]
    passed_df = passed_df[(passed_df["Total_Length"] >= lower_percentiles["Total_Length"]) &
            (passed_df["Total_Length"] <= upper_percentiles["Total_Length"])]
    passed_df = passed_df[(passed_df["Avg_Distances"] >= lower_percentiles["Avg_Distances"]) &
            (passed_df["Avg_Distances"] <= upper_percentiles["Avg_Distances"])]

    passed_files = os.path.join(fasta_dir, "passed.txt")
    if os.path.isfile(passed_files):
        os.remove(passed_files)
    for accession in passed_df.index:
        with open(passed_files, "a") as f:
            f.write(accession+"\n")

def bool_df_for_failed(fasta_dir, stats_df, max_n_count, max_contigs, lower_percentile, upper_percentile):
    quantiles_stats = stats_df.quantile([lower_percentile, upper_percentile])
    quantiles_stats.index = [lower_percentile, upper_percentile]
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]

    n_count_bool = (stats_df.iloc[:]["N_Count"]) <= max_n_count
    contigs_bool = (stats_df.iloc[:]["Contigs"]) <= max_contigs

    distances_bool = (stats_df.iloc[:]["Avg_Distances"] >= lower_percentiles["Avg_Distances"]) \
                   & (stats_df.iloc[:]["Avg_Distances"] <= upper_percentiles["Avg_Distances"])

    lengths_bool = (stats_df.iloc[:]["Total_Length"] >= lower_percentiles["Total_Length"]) \
                   & (stats_df.iloc[:]["Total_Length"] <= upper_percentiles["Total_Length"])

    for i in stats_df.index:
        N_Count = stats_df.loc[i, "N_Count"]
        Contigs = stats_df.loc[i, "Contigs"]
        Avg_Distances = stats_df.loc[i, "Avg_Distances"]
        Total_Length = stats_df.loc[i, "Total_Length"]

        if N_Count > max_n_count:
            stats_df.loc[i, "N_Count"] = "*{}*".format(stats_df.loc[i, "N_Count"])
        if Contigs > max_contigs:
            stats_df.loc[i, "Contigs"] = "*{}*".format(stats_df.loc[i, "Contigs"])
        if Avg_Distances <= lower_percentiles["Avg_Distances"] or Avg_Distances >= upper_percentiles["Avg_Distances"]:
            stats_df.loc[i, "Avg_Distances"] = "*{}*".format(stats_df.loc[i, "Avg_Distances"])
        if Total_Length <= lower_percentiles["Total_Length"] or Total_Length >= upper_percentiles["Total_Length"]:
            stats_df.loc[i, "Total_Length"] = "*{}*".format(stats_df.loc[i, "Total_Length"])

    bool_df = pd.concat([lengths_bool, n_count_bool, contigs_bool, distances_bool], axis=1)
    bool_df.to_csv(os.path.join(fasta_dir, "failed_tf.csv"))
    stats_df.to_csv(os.path.join(fasta_dir, "failed.csv"))

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=2000)
    parser.add_argument("-c", "--max_contigs", help = "Maximum number of contigs acceptable", type=int, default=500)
    parser.add_argument("-l", "--lower", help = "Enter the lower percentile range to filter fastas.", type=float, default=.05)
    parser.add_argument("-u", "--upper", help = "Enter the upper percentile range to filter fastas.", type=float, default=.9)
    parser.add_argument("-m", "--mash", help = "Create a sketch file and distance matrix", action="store_true")
    args = parser.parse_args()

    fasta_dir = args.fasta_dir
    distance_matrix = os.path.join(fasta_dir, "distance_matrix.csv")

    if args.mash:
        mash(fasta_dir)

#   clean_up_matrix_and_get_quantiles(distance_matrix)
    stats_df = generate_fasta_stats(fasta_dir, distance_matrix)
    assess_stats_df(fasta_dir, stats_df, args.max_n_count, args.max_contigs, args.lower, args.upper)
    bool_df_for_failed(fasta_dir, stats_df, args.max_n_count, args.max_contigs, args.lower, args.upper)

Main()
