#!/usr/bin/env python

import os, argparse
from subprocess import Popen
from shutil import rmtree
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

    return distance_matrix

def clean_up_matrix(distance_matrix):

    """
    Sets indices and headers to the accession ID's
    """

    distances_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    new_index = []
    for i in distances_df.index:
        name = i.split("/")[-1]
        name = name.split("_")[0:2]
        new_index.append("_".join(name))

    distances_df.index = new_index
    distances_df.columns = new_index
    distances_df.to_csv(distance_matrix, sep="\t")

def dist_matrix_quantiles():

  # distances_pass_fail = os.path.join(fasta_dir, "distance_matrix_pass_fail.csv")
  # distances_upper_percentile = distances_df.quantile(upper_percentile)
  # pass_fail_distances = pd.DataFrame(index=distances_df.index)

  # for accession in distances_df.columns:
  #     for distance in distances_df.iloc[:][accession]:
  #         if distance > distances_upper_percentile.loc[accession]:
  #             pass_fail_distances[accession] = "F"
  #         else:
  #             pass_fail_distances[accession] = "P"

  # distances_df.to_csv(distances_pass_fail, sep="\t")
  pass

def generate_fasta_stats(fasta_dir, distance_matrix):

    distances_df = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    avg_distances = ["{:05.4f}".format(i) for i in distances_df.mean()]

    for root, dirs, files, in os.walk(fasta_dir):
        file_names = []
        contig_totals = []
        lengths = []
        n_counts = []

        for f in files:
            if f.endswith(".fasta"):
                print("Getting stats for {}".format(f))
                fasta = (os.path.join(root, f))
                file_names.append(fasta)

                # Read all contigs for current fasta into list
                try:
                    contigs = [ seq.seq for seq in SeqIO.parse(fasta, "fasta") ]
                except UnicodeDecodeError:
                    print("{} threw UnicodeDecodeError".format(f))
                    with open(os.path.join(fasta_dir, "quality_filters_log.txt"), "a") as log:
                        log.write("{} threw UnicodeDecodeError".format(f))

                # Append the total number of contigs to contig_totals
                contig_totals.append(len(contigs))

                # Read the length of each contig into a list
                contig_lengths = [ len(str(seq)) for seq in contigs ]
                # Append the sum of all contig lengths to lengths
                lengths.append(sum(contig_lengths))

                # Read the N_Count for each contig into a list
                N_Count = [ int(str(seq.upper()).count("N")) for seq in contigs ]
                # Append the total N_Count to n_counts
                n_counts.append(sum(N_Count))

        SeqDataSet = list(zip(contig_totals, lengths, n_counts, avg_distances))
        stats_df = pd.DataFrame(data = SeqDataSet, index=file_names, columns=["Contigs", "Total_Length", "N_Count", "Avg_Distances"], dtype="float64")
        stats_df.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")

        return stats_df

def assess_stats_df(fasta_dir, stats_df, max_n_count, max_contigs, lower_percentile, upper_percentile, stds):


    print("Stats:", end="\n\n")
    print(stats_df)
    quantiles_stats = stats_df.quantile([lower_percentile, upper_percentile])
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]

    # Filter based on N_Count and Contigs first
    passed_df = stats_df[stats_df["N_Count"] <= max_n_count ]

    passed_df = passed_df[passed_df["Contigs"] <= max_contigs ]

    def calculate_homogeneity():
        
        """
        Determine whether or not the dataset is homogeneous, 
        i.e. the mean, median, max, and min values
        are within a 1.5 stds of each other.
        """

        # Get the absolute difference between mean/median, mean/max, and mean/min.
        abs_dif_mean_median = abs(passed_df.mean()["Total_Length"] - passed_df.median()["Total_Length"]) 
        abs_dif_mean_max = abs(passed_df.mean()["Total_Length"] - passed_df.max()["Total_Length"])
        abs_dif_mean_min = abs(passed_df.mean()["Total_Length"] - passed_df.min()["Total_Length"])
        # Calculate the standard deviation for Total_Length; multiple by stds
        length_stds = (passed_df.std()["Total_Length"]) * stds

        # Check if absolute differences are <= std * stds
        check_median = int(abs_dif_mean_median <= (length_stds))
        check_max = int(abs_dif_mean_max <= (length_stds))
        check_min = int(abs_dif_mean_max <= (length_stds))
        homogeneity_rating = check_median + check_max + check_min

        print("abs_dif_mean_median:  {}".format(abs_dif_mean_median))
        print("abs_dif_mean_max:  {}".format(abs_dif_mean_max))
        print("abs_dif_mean_min:  {}".format(abs_dif_mean_min))
        print("length_stds:  {}".format(length_stds))
        print("Median results:  {}".format(check_median))
        print("Max results:  {}".format(check_max))
        print("Min results:  {}".format(check_min))

        return int(homogeneity_rating)

    homogeneity_rating = calculate_homogeneity()

    if homogeneity_rating < 3:
        print("Homogeneity rating:  {}".format(homogeneity_rating))
        print("Dataset not homogeneous.  Filtering to commence based on percentiles.")
        passed_df = passed_df[(passed_df["Total_Length"] >= lower_percentiles["Total_Length"]) &
                (passed_df["Total_Length"] <= upper_percentiles["Total_Length"])]

        passed_df = passed_df[(passed_df["Avg_Distances"] >= lower_percentiles["Avg_Distances"]) &
                (passed_df["Avg_Distances"] <= upper_percentiles["Avg_Distances"])]
    else:
        print("Data set is homoegenous")

    print("Passed:", end="\n\n")
    print(passed_df)
    return passed_df

def generate_links_to_passed(passed_df, passed_dir):

    for source in passed_df.index:
        print("{} passed.".format(source.split("/")[-1]))
        dst = os.path.join(passed_dir, source.split("/")[-1])
        print("Linking to {}".format(dst))
        os.link(source, dst)

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

def make_passed_dir(fasta_dir):
    
    passed_dir = os.path.join(fasta_dir, "passed")
    if os.path.isdir(passed_dir):
        rmtree(passed_dir)
    if not os.path.isdir(passed_dir):
        os.mkdir(passed_dir)

    return passed_dir
    

def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("fasta_dir", help = "directory containing your FASTA files")
    parser.add_argument("-d", "--directory", help = "Specifies that `fasta_dir` is a parent directory with \
            subdirectories for each species", action="store_true")
    parser.add_argument("--from_list", help = "Specify a list of one more more directories to fun filters on.", nargs="+")
    parser.add_argument("-f", "--file", help = "Specify a file containing one species directory per line.", action="store_true")
    parser.add_argument("-n", "--max_n_count", help = "Maximum number of N's acceptable", type=int, default=2000)
    parser.add_argument("-c", "--max_contigs", help = "Maximum number of contigs acceptable", type=int, default=500)
    parser.add_argument("-l", "--lower", help = "Enter the lower percentile range to filter fastas.", type=float, default=.05)
    parser.add_argument("-u", "--upper", help = "Enter the upper percentile range to filter fastas.", type=float, default=.9)
    parser.add_argument("-s", "--stds", help = "The number of standard deviations used to define the acceptable distance \
            between the mean and the median, min, and max value, when judging the homogeneity of the dataset", type=float, default=1.5)
    parser.add_argument("-m", "--mash", help = "Create a sketch file and distance matrix", action="store_true")
    args = parser.parse_args()

    fasta_dir = args.fasta_dir

    if args.from_list:
        for name in args.from_list:
            species_dir = os.path.join(fasta_dir, name)
            distance_matrix = os.path.join(species_dir, "distance_matrix.csv")
            mash(species_dir)
            stats_df = generate_fasta_stats(species_dir, distance_matrix)
            passed_df = assess_stats_df(species_dir, stats_df, args.max_n_count, args.max_contigs, args.lower, args.upper, args.stds)
            passed_dir = make_passed_dir(species_dir)
            generate_links_to_passed(passed_df, passed_dir)
    else:
        for name in os.listdir(fasta_dir):
            if os.path.join(fasta_dir, name):
                species_dir = os.path.join(fasta_dir, name)
                distance_matrix = os.path.join(species_dir, "distance_matrix.csv")
                mash(species_dir)
                stats_df = generate_fasta_stats(species_dir, distance_matrix)
                passed_df = assess_stats_df(species_dir, stats_df, args.max_n_count, args.max_contigs, args.lower, args.upper, args.stds)
                passed_dir = make_passed_dir(species_dir)
                generate_links_to_passed(passed_df, passed_dir)

Main()
