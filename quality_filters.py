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
    summary = os.path.join(FASTA_dir, "summary.txt")

    files = [sketch_file, distance_matrix, filter_log, summary]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

def mash(FASTA_dir, mash_exe):
    sketch_file = os.path.join(FASTA_dir, "all.msh")
    all_fastas = os.path.join(FASTA_dir, "*.fasta")
    distance_matrix = os.path.join(FASTA_dir, "distance_matrix.csv")
    sketch_command = "{} sketch -o {} {}".format(mash_exe, sketch_file, all_fastas)
    distance_command = "{} dist -p 4 -t {} {} > {}".format(mash_exe, sketch_file, all_fastas, distance_matrix)
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
                N_Count = [len(findall("[^ATCG]", str(seq))) for seq in contigs]
                # Append the total N_Count to n_counts
                n_counts.append(sum(N_Count))

        SeqDataSet = list(zip(n_counts, contig_totals, assembly_sizes))
        stats = pd.DataFrame(data=SeqDataSet, index=file_names, columns=["N_Count", "Contigs", "Assembly_Size"], dtype="float64")
        mean_distances = distance_matrix.mean()
        stats["MASH"] = mean_distances
        stats.to_csv(os.path.join(root, "stats.csv"), index_label="Accession")

        return stats

def filter_med_ad(FASTA_dir, stats, multiplier, max_n_count):

    # Create new log files and remove old ones
    passed_log = os.path.join(FASTA_dir, "passed_{}.csv".format(multiplier))
    failed_log = os.path.join(FASTA_dir, "failed_{}.csv".format(multiplier))

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

    # Filter using special function for contigs
    #if not len(passed_I) <=5:
    passed_II = filter_contigs(stats, passed_I, multiplier, failed_df)

    Assembly_Size_med_ad = abs(passed_II["Assembly_Size"] - passed_II["Assembly_Size"].median()).mean()# Median absolute deviation
    deviation_reference = Assembly_Size_med_ad * multiplier
    passed_III = passed_II[abs(passed_II["Assembly_Size"] - passed_II["Assembly_Size"].median()) <= deviation_reference]
    for i in passed_II.index:
        if i not in passed_III.index:
            failed_df["Assembly_Size"][i] = stats["Assembly_Size"][i]

    MASH_med_ad = abs(passed_III["MASH"] - passed_III["MASH"].median()).mean()# Median absolute deviation
    deviation_reference = MASH_med_ad * multiplier
    passed_IV = passed_III[abs(passed_III["MASH"] - passed_III["MASH"].median()) <= deviation_reference]
    for i in passed_III.index:
        if i not in passed_IV.index:
            failed_df["MASH"][i] = stats["MASH"][i]

    failed_df.dropna(axis=0, how="all", inplace=True)
    summarize_results(FASTA_dir, failed_df, multiplier)
    failed_df.fillna(value="N/A", inplace=True)
    failed_df.to_csv(failed_log, header=True)
    passed_IV.to_csv(passed_log, header=True)

def filter_contigs(stats, passed_I, multiplier, failed_df):

    contigs = passed_I["Contigs"]
    contigs_above_median = contigs[contigs >= contigs.median()]
    contigs_below_median = contigs[contigs <= contigs.median()]
    contigs_lower = contigs[contigs <= 10] # Save genomes with < 10 contigs to add them back in later.
    contigs = contigs[contigs > 10] # Only look at genomes with > 10 contigs to avoid throwing off the Median AD
  # Should this value be obtained from contigs > 10 or contigs >contigs.median()?
    contigs_med_ad = abs(contigs - contigs.median()).mean() # Median absolute deviation
    contigs_deviation_reference = contigs_med_ad * multiplier
  # failed_df.insert(2, "Med_AD", contigs_med_ad)
  # ["Med_AD"] = contigs_med_ad
    contigs = contigs[abs(contigs - contigs.median()) <= contigs_deviation_reference]
    contigs = pd.concat([contigs, contigs_lower])
  # contigs = pd.concat([contigs, contigs_lower, contigs_below_median])

    # Avoid returning empty DataFrame when no genomes are removed above
    if len(contigs) == len(passed_I):
        passed_II = passed_I
    else:
        failed_contigs = [i for i in passed_I.index if i not in contigs.index]
        passed_II = passed_I.drop(failed_contigs)

        for i in failed_contigs:
            failed_df["Contigs"][i] = stats["Contigs"][i]

    return passed_II

def summarize_results(FASTA_dir, failed_df, assembly_med_ad, mash_med_ad, multiplier):
    summary = os.path.join(FASTA_dir, "summary.txt")
    failed_N = len(failed_df[failed_df.N_Count.notnull()])
    failed_contigs = len(failed_df[failed_df.Contigs.notnull()])
    failed_assembly_size = len(failed_df[failed_df.Assembly_Size.notnull()])
    failed_mash = len(failed_df[failed_df.MASH.notnull()])
    total_failed = failed_N+failed_contigs+failed_assembly_size+failed_mash
    with open (summary, "a") as f:
        f.write("Results from multiplying the deviation reference value by {}\n".format(multiplier))
        f.write("Total filtered genomes: {}\n".format(total_failed))
        f.write("N Count: {} genomes filtered\n".format(failed_N))
        f.write("Contigs: {} genomes filtered\n".format(failed_contigs))
        f.write("Assembly Size: {} genomes filtered\n".format(failed_assembly_size))
        f.write("MASH: {} genomes filtered\n".format(failed_mash))
        f.write("#"*61+"\n")

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
    parser.add_argument("-x", "--mash_exe", help = "Path to MASH")
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

    mash_exe = "/common/contrib/bin/mash-Linux64-v1.1/mash"
    if args.mash_exe:
        mash_exe = args.mash_exe

    if args.directories:
        for name in args.FASTA_dir:
            FASTA_dir = name
            if len(os.listdir(FASTA_dir)) <= 5: # pass if there are <= 5 FASTA's
                continue
            clean_up(FASTA_dir)
            distance_matrix = mash(FASTA_dir, mash_exe)
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
                    distance_matrix = mash(FASTA_dir, mash_exe)
                    stats = generate_fasta_stats(FASTA_dir, distance_matrix)
                    for num in [2, 2.5, 3, 3.5]:
                        filter_med_ad(FASTA_dir, stats, multiplier=num, max_n_count=100)
                  # passed_dir = make_passed_dir(FASTA_dir)
                  # generate_links_to_passed(passed_df, passed_dir, FASTA_dir)

Main()
