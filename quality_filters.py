#!/usr/bin/env python

import os, argparse
from subprocess import Popen
from shutil import rmtree
import pandas as pd
from Bio import SeqIO
from re import findall

def clean_up(FASTA_dir):
    sketch_file = os.path.join(FASTA_dir, "all.msh")
    distance_matrix = os.path.join(FASTA_dir, "distance_matrix.csv")
    filter_log = os.path.join(FASTA_dir, "filter_log.txt")

    files = [sketch_file, distance_matrix, filter_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

def mash(FASTA_dir):
    mash="/home/asanchez/local/bin/mash-Linux64-v1.1/mash"
    sketch_file = os.path.join(FASTA_dir, "all.msh")
    all_fastas = os.path.join(FASTA_dir, "*.fasta")
    distance_matrix = os.path.join(FASTA_dir, "distance_matrix.csv")
    sketch_command = "{} sketch -o {} {}".format(mash, sketch_file, all_fastas)
    distance_command = "{} dist -p 4 -t {} {} > {}".format(mash, sketch_file, all_fastas, distance_matrix)
    Popen(sketch_command, shell="True").wait()
    Popen(distance_command, shell="True").wait()
    distance_matrix = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    clean_up_matrix(FASTA_dir, distance_matrix)

    return distance_matrix

def clean_up_matrix(FASTA_dir, distance_matrix):

    """
    Set indices and headers to the accession ID's
    """

    new_index = []
    for i in distance_matrix.index:
        name = i.split("/")[-1].strip(".fasta")
        new_index.append(name)

    distance_matrix.index = new_index
    distance_matrix.columns = new_index
    distance_matrix.to_csv(os.path.join(FASTA_dir, "distance_matrix.csv"), sep="\t")

def generate_fasta_stats(FASTA_dir, distance_matrix):

    med_ad_distances = {}
    for i in distance_matrix:
        med_ad_distances[i] = abs(distance_matrix[i] - distance_matrix[i].median()).mean()

    med_ad_distances = pd.Series(med_ad_distances)

    for root, dirs, files, in os.walk(FASTA_dir):
        file_names = []
        contig_totals = []
        assembly_sizes = []
        n_counts = []

        for f in files:
            if f.endswith(".fasta"):
                print("Getting stats for {}".format(f))
                fasta = (os.path.join(root, f))
                name = fasta.split("/")[-1].strip(".fasta")
                file_names.append(name)

                # Read all contigs for current fasta into list
                try:
                    contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                except UnicodeDecodeError:
                    print("{} threw UnicodeDecodeError".format(f))
                    with open(os.path.join(FASTA_dir, "filter_log.txt"), "a") as log:
                        log.write("{} threw UnicodeDecodeError\n".format(f))

                # Append the total number of contigs to contig_totals
                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                assembly_size = [ len(str(seq)) for seq in contigs ]
                # Append the sum of all contig lengths to lengths
                assembly_sizes.append(sum(assembly_size))

                # Read the N_Count for each contig into a list
              # N_Count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
              # N_Count = [sum(1 for N in str(seq.upper()))for seq in contigs]
                N_Count = [len(findall("[^ATCG]", str(seq))) for seq in contigs]
                # Append the total N_Count to n_counts
                n_counts.append(sum(N_Count))

        SeqDataSet = list(zip(assembly_sizes, contig_totals, n_counts))
        stats = pd.DataFrame(data = SeqDataSet, index=file_names, columns=["Assembly_Size", "Contigs", "N_Count"], dtype="float64")
        stats["MASH"] = med_ad_distances
        stats.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")

        return stats

def filter_med_ad(FASTA_dir, stats, multiplier, max_n_count):

    # Create new log files and remove old ones
    passed_log = "passed_{}.txt".format(multiplier)
    passed_log = os.path.join(FASTA_dir, passed_log)
    failed_log = "failed_{}.txt".format(multiplier)
    failed_log = os.path.join(FASTA_dir, failed_log)
    files = [passed_log, failed_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

    # Filter based on N's first
    passed_I = stats[stats["N_Count"] <= max_n_count ]
    failed_df = pd.DataFrame(index=stats.index, columns=stats.columns)
    for i in stats.index:
        if i not in passed_I.index:
            failed_df["N_Count"][i] = stats["N_Count"][i]

    # Deal with contigs separately
    passed_II = filter_contigs(stats, passed_I, multiplier, failed_df, failed_log)

    axes = ["Assembly_Size", "MASH"]
    for axis in axes:
        med_ad = abs(passed_II[axis] - passed_II[axis].median()).mean()# Median absolute deviation
        deviation_reference = med_ad * multiplier
        passed_III = passed_II[abs(passed_II[axis] - passed_II[axis].median()) <= deviation_reference]
        for i in passed_II.index:
            if i not in passed_III.index:
                failed_df[axis][i] = stats[axis][i]

    failed_df.dropna(axis=0, how="all", inplace=True)
    failed_df.fillna(value="N/A", inplace=True)
    failed_df.to_csv(failed_log, header=True)
    passed_III.to_csv(passed_log, header=True)

def filter_contigs(stats, passed_I, multiplier, failed_df, failed_log):

    contigs = passed_I["Contigs"]
    contigs = contigs[contigs > 10] # Genomes with < 10 contigs automatically pass.
    contigs_lower = contigs[contigs < 10] # Save genomes with < 10 contigs to add them back in later.
    contigs_med_ad = abs(contigs - contigs.median()).mean() # Median absolute deviation
    contigs_deviation_reference = contigs_med_ad * multiplier
    # should this be <= contigs_deviation_reference] OR < contigs.median() ?
    # in order to avoid filtering genomes that have too few contigs?
    contigs = contigs[abs(contigs - contigs.median()) <= contigs_deviation_reference]
    contigs = pd.concat([contigs, contigs_lower])
    failed_contigs = [i for i in passed_I.index if i not in contigs.index]

    for i in passed_I.index:
        if len(contigs) == len(passed_I):
            passed_II = passed_I
        elif i not in contigs.index:
            passed_II = passed_I.drop(i) # Update the passed_I DataFrame
            failed_df["Contigs"][i] = stats["Contigs"][i]

    return passed_II

def generate_links_to_passed(passed_df, passed_dir, FASTA_dir):

    filter_log = os.path.join(FASTA_dir, "filter_log.txt")
    with open(filter_log, "a") as log:
        log.write("{} FASTA's passed".format(len(passed_df)))

    for source in passed_df.index:
        print("{} passed.".format(source.split("/")[-1]))
        dst = os.path.join(passed_dir, source.split("/")[-1])
        print("Linking to {}".format(dst))
        os.link(source, dst)

def make_passed_dir(FASTA_dir):
    
    passed_dir = os.path.join(FASTA_dir, "passed")
    if os.path.isdir(passed_dir):
        rmtree(passed_dir)
    if not os.path.isdir(passed_dir):
        os.mkdir(passed_dir)

    return passed_dir
    
def bool_df_for_failed(FASTA_dir, stats, max_n_count, max_contigs, lower_percentile, upper_percentile):

    quantiles_stats = stats.quantile([lower_percentile, upper_percentile])
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]

    n_count_bool = (stats.iloc[:]["N_Count"]) <= max_n_count
    contigs_bool = (stats.iloc[:]["Contigs"]) <= max_contigs

    distances_bool = (stats.iloc[:]["Avg_Distances"] >= lower_percentiles["Avg_Distances"]) \
                   & (stats.iloc[:]["Avg_Distances"] <= upper_percentiles["Avg_Distances"])

    lengths_bool = (stats.iloc[:]["Assembly_Size"] >= lower_percentiles["Total_Length"]) \
                   & (stats.iloc[:]["Assembly_Size"] <= upper_percentiles["Total_Length"])

    passed = os.listdir(os.path.join(FASTA_dir, passed))
    total = os.listdir(FASTA_dir)
    filter_log = os.path.join(FASTA_dir, "filter_log.txt")
    with open(filter_log, "a") as log:
        log.write("{} passed out of {}".format(len(passed), len(total)))


  # for i in stats.index:
  #     N_Count = stats.loc[i, "N_Count"]
  #     Contigs = stats.loc[i, "Contigs"]
  #     Avg_Distances = stats.loc[i, "Avg_Distances"]
  #     Assembly_Size = stats.loc[i, "Total_Length"]

  #     if N_Count > max_n_count:
  #         stats.loc[i, "N_Count"] = "*{}*".format(stats.loc[i, "N_Count"])
  #     if Contigs > max_contigs:
  #         stats.loc[i, "Contigs"] = "*{}*".format(stats.loc[i, "Contigs"])
  #     if Avg_Distances <= lower_percentiles["Avg_Distances"] or Avg_Distances >= upper_percentiles["Avg_Distances"]:
  #         stats.loc[i, "Avg_Distances"] = "*{}*".format(stats.loc[i, "Avg_Distances"])
  #     if Assembly_Size <= lower_percentiles["Total_Length"] or Total_Length >= upper_percentiles["Total_Length"]:
  #         stats.loc[i, "Assembly_Size"] = "*{}*".format(stats.loc[i, "Total_Length"])

    bool_df = pd.concat([lengths_bool, n_count_bool, contigs_bool, distances_bool], axis=1)
    bool_df.to_csv(os.path.join(FASTA_dir, "failed_tf.csv"))
    stats.to_csv(os.path.join(FASTA_dir, "failed.csv"))

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("FASTA_dir", help = "directory containing your FASTA files", nargs="+")
    parser.add_argument("-d", "--directories", help = "The complete path to one or more directories containing\
            FASTA's you want to run the filters on.", action="store_true")
    parser.add_argument("-p", "--parent_dir", help = "The parent directory containing subdirectories with FASTA's for\
            each of the collections you want to run the filters on.")
    parser.add_argument("--from_list", help = "Specify a list of one more more directories to fun filters on.", nargs="+")
    parser.add_argument("-f", "--file", help = "Specify a file containing one species directory per line.", action="store_true")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=100)
    parser.add_argument("-c", "--max_contigs", help = "Maximum number of contigs acceptable", type=int, default=500)
    parser.add_argument("-l", "--lower", help = "Enter the lower percentile range to filter fastas.", type=float, default=.05)
    parser.add_argument("-u", "--upper", help = "Enter the upper percentile range to filter fastas.", type=float, default=.9)
    parser.add_argument("-s", "--multiplier", help = "The number of standard deviations used to define the acceptable distance \
            between the mean and the median, min, and max value, when judging the homogeneity of the dataset", type=float, default=1.5)
    parser.add_argument("-m", "--mash", help = "Create a sketch file and distance matrix", action="store_true")
    parser.add_argument("--deviation_type", help = "Which measurement to be used for the deviation reference point.\
            options are `stds`, `med_ads`, and `mads`", type=str, default="med_ads")
    args = parser.parse_args()

    max_ns = args.max_n_count

    if args.directories:
        for name in args.FASTA_dir:
            FASTA_dir = name
            if len(os.listdir(FASTA_dir)) <= 5: # pass if there are <= 5 FASTA's
                continue
            clean_up(FASTA_dir)
            distance_matrix = mash(FASTA_dir)
            stats = generate_fasta_stats(FASTA_dir, distance_matrix)
            for num in [2, 2.5, 3, 3.5]:
                filter_med_ad(FASTA_dir, stats, multiplier=num, max_n_count=100)
    else:
        root = args.FASTA_dir[0]
        for name in os.listdir(root):
            FASTA_dir = os.path.join(root, name)
            if os.path.isdir(FASTA_dir):
                if len(os.listdir(FASTA_dir)) <= 5:
                    continue
                else:
                    clean_up(FASTA_dir)
                    distance_matrix = mash(FASTA_dir)
                    stats = generate_fasta_stats(FASTA_dir, distance_matrix)
                    for num in [2, 2.5, 3, 3.5]:
                        filter_med_ad(FASTA_dir, stats, multiplier=num, max_n_count=100)
                  # passed_dir = make_passed_dir(FASTA_dir)
                  # generate_links_to_passed(passed_df, passed_dir, FASTA_dir)

Main()
