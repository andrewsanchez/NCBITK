#!/usr/bin/env python

import os, argparse
from subprocess import Popen
from shutil import rmtree
import pandas as pd
from Bio import SeqIO
from re import findall

def clean_up(species_dir):
    sketch_file = os.path.join(species_dir, "all.msh")
    distance_matrix = os.path.join(species_dir, "distance_matrix.csv")
    filter_log = os.path.join(species_dir, "filter_log.txt")

    files = [sketch_file, distance_matrix, filter_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

def mash(species_dir):
    mash="/home/asanchez/local/bin/mash-Linux64-v1.1/mash"
    sketch_file = os.path.join(species_dir, "all.msh")
    all_fastas = os.path.join(species_dir, "*.fasta")
    distance_matrix = os.path.join(species_dir, "distance_matrix.csv")
    sketch_command = "{} sketch -o {} {}".format(mash, sketch_file, all_fastas)
    distance_command = "{} dist -p 4 -t {} {} > {}".format(mash, sketch_file, all_fastas, distance_matrix)
    Popen(sketch_command, shell="True").wait()
    Popen(distance_command, shell="True").wait()
    distance_matrix = pd.read_csv(distance_matrix, index_col=0, delimiter="\t")
    clean_up_matrix(species_dir, distance_matrix)

    return distance_matrix

def clean_up_matrix(species_dir, distance_matrix):

    """
    Set indices and headers to the accession ID's
    """

    new_index = []
    for i in distance_matrix.index:
        name = i.split("/")[-1].strip(".fasta")
        new_index.append(name)

    distance_matrix.index = new_index
    distance_matrix.columns = new_index
    distance_matrix.to_csv(os.path.join(species_dir, "distance_matrix.csv"), sep="\t")

def generate_fasta_stats(species_dir, distance_matrix):

    med_ad_distances = {}
    for i in distance_matrix:
        med_ad_distances[i] = abs(distance_matrix[i] - distance_matrix[i].median()).mean()

    med_ad_distances = pd.Series(med_ad_distances)

    for root, dirs, files, in os.walk(species_dir):
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
                    with open(os.path.join(species_dir, "filter_log.txt"), "a") as log:
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

def filter_med_ad(species_dir, stats, multiplier, max_n_count):

    passed_log = "passed_{}.txt".format(multiplier)
    passed_log = os.path.join(species_dir, passed_log)
    failed_log = "failed_{}.txt".format(multiplier)
    failed_log = os.path.join(species_dir, failed_log)
    files = [passed_log, failed_log]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)

    failed_df = pd.DataFrame(index=stats.index, columns=stats.columns)
    passed = stats[stats["N_Count"] <= max_n_count ] # Filter based on N's first
    failed_n = [i for i in stats.index if i not in passed.index]
    for i in failed_n:
        failed_df["N_Count"][i] = stats["N_Count"][i]
    filter_contigs(stats, passed, multiplier, failed_df, failed_log) # Deal with contigs separately

    axes = ["Assembly_Size", "MASH"]
    for axis in axes:
        before = passed.index
        med_ad = abs(passed[axis] - passed[axis].median()).mean()# Median absolute deviation
        deviation_reference = med_ad * multiplier
        passed = passed[abs(passed[axis] - passed[axis].median()) <= deviation_reference]
        failed = [i for i in before if i not in passed.index]
        for i in failed:
            failed_df[axis][i] = stats[axis][i]
    failed_df.dropna(axis=0, how="all", inplace=True)
    failed_df.fillna(value="N/A", inplace=True)
    failed_df.to_csv(failed_log, header=True)
    passed.to_csv(passed_log, header=True)

def filter_contigs(stats, passed, multiplier, failed_df, failed_log):

    contigs = passed["Contigs"]
    contigs = contigs[contigs > 10] # Genomes with < 10 contigs automatically pass.
    contigs_lower = contigs[contigs < 10] # Save genomes with < 10 contigs to add them back in later.
    contigs_med_ad = abs(contigs - contigs.median()).mean()# Median absolute deviation
    contigs_deviation_reference = contigs_med_ad * multiplier
    # should this be <= contigs_deviation_reference] OR < contigs.median() ?
    contigs = contigs[abs(contigs - contigs.median()) <= contigs_deviation_reference]
    contigs = pd.concat([contigs, contigs_lower])
    failed_contigs = [i for i in passed.index if i not in contigs.index]

    for i in failed_contigs:
        passed.drop(i, inplace=True, errors="ignore") # Update the passed DataFrame
        failed_df["Contigs"][i] = stats["Contigs"][i]

def assess_failed(stats, failed_ids, failed_log, passed, axes, multiplier):

    with open(failed_log, "a") as log:
        log.write("Failed:  {} out of {} = {:4.2f}%\n".format(len(failed_ids), len(stats), len(failed_ids)/len(stats)*100))
        for axis in axes:
            med_ad = abs(passed[axis] - passed[axis].median()).mean()# Median absolute deviation
            deviation_reference = med_ad * multiplier
            for i in failed:
                value = passed.ix[i][axis]
                if not abs(value - passed[axis].median()) <= deviation_reference:
                    accession = "_".join(i.split("_")[:2])
                    log.write("{} {}:  {} > {}\n".format(accession, axis, value, deviation_reference))

def assess_stats(species_dir, stats, max_n_count, max_contigs, lower_percentile, upper_percentile, multiplier=1.5, deviation="med_ads"):

    log_file = "{}_filter_log.txt".format(deviation)
    filter_log = os.path.join(species_dir, log_file)
    if os.path.isfile(filter_log):
        os.remove(filter_log)

    axes = ["Assembly_Size", "Contigs"]
    for axis in axes:
        if axis == "Contigs":
            passed_df = stats[stats["Contigs"] >= 10 ]
            passed_df = stats[stats["N_Count"] <= max_n_count ] # Filter based on N_Count first
        else:
            passed_df = stats[stats["N_Count"] <= max_n_count ] # Filter based on N_Count first

    # Get the absolute difference between mean/median, mean/max, and mean/min.
    # Maybe get these values, not just for median, max, and min, but for the range of values
    # that fall within 5% of the median, max, and min?  Then take the average of those?
    median_mean_ad = abs(passed_df.median()[axis] - passed_df.mean()[axis]) 
    median_max_ad = abs(passed_df.median()[axis] - passed_df.max()[axis])
    median_min_ad = abs(passed_df.median()[axis] - passed_df.min()[axis])
    # The average or STD between the above three measures might be helpful.

    def calculate_homogeneity(axis, passed_df, deviation="med_ads"):
        
        """
        Determine whether or not the dataset is homogeneous, 
        i.e. the mean, median, max, and min values
        are within `std * multiplier` of each other.
        Axis is the column of the DataFrame you wish to check.
        """

        # Check if absolute differences are <= deviation * multiplier
        # Calculate reference points for "normal" deviation checks
        # Peform these checks for the range of values referenced above?
        def create_log(check_median, check_max, check_min, homogeneity_rating):
            with open(filter_log, "a") as log:
                log.write("Results for {} using {} as the deviation reference point:\n".format(axis, deviation))
                log.write("abs(median - mean) = {}\n".format(median_mean_ad))
                log.write("abs(median - max) = {}\n".format(median_max_ad))
                log.write("abs(median - min) = {}\n".format(median_min_ad))
                log.write("{} * multiplier:\n".format(deviation))
                if deviation == "std":
                    log.write("{} * {} = {}\n".format(std, multiplier, stds))
                elif deviation == "med_ad":
                    log.write("{} * {} = {}\n".format(med_ad, multiplier, med_ads))
                elif deviation == "mad":
                    log.write("{} * {} = {}\n".format(mad, multiplier, mads))
                log.write("1 if abs(x-y) <= deviation * multiplier, 0 zero if not.\n")
                log.write("Median results:  {}\n".format(check_median))
                log.write("Max results:  {}\n".format(check_max))
                log.write("Min results:  {}\n".format(check_min))
                log.write("Homogeneity rating:  {}\n".format(homogeneity_rating))
                log.write("{} total FASTA's".format(len(stats)))

        def med_max_min_checks():
            check_median = int(median_mean_ad <= (deviation_reference))
            check_max = int(median_max_ad <= (deviation_reference))
            check_min = int(median_min_ad <= (deviation_reference))
            homogeneity_rating = check_median + check_max + check_min
            create_log(check_median, check_max, check_min, homogeneity_rating)

        if deviation == "stds":
            std = passed_df.std()[axis] #  Standard deviation
            deviation_reference = std * multiplier
            med_max_min_checks()
        elif deviation == "med_ads":
            med_ad = abs(passed_df - passed_df.median()).mean()# Median absolute deviation
            deviation_reference = med_ad[axis] * multiplier
            med_max_min_checks()
        elif deviation == "mads":
            mad = abs(stats - stats.mean()).mean()# Mean absolute deviation
            deviation_reference = mad[axis] * multiplier
            med_max_min_checks()

            # If this is used: `passed = passed_df[abs(passed_df[axis] - med_ad) <= deviation_reference]`
            # The code below becomes completely unecessary

        if axis == "Contigs":
            if not check_max:
                log.write("Filtering FASTA's above the upper percentile")
                return "upper"
            elif check_max:
                log.write("Dataset is homogenous")
                return True
        else:
            if check_max and not check_min:
                if homogeneity_rating == 2:
                    log.write("Filtering FASTA's below the lower percentile")
                    return "lower"
                else:
                    log.write("Filtering FASTA's outside the upper and lower percentiles")
                    return False
            elif check_min and not check_max:
                if homogeneity_rating == 2:
                    log.write("Filtering FASTA's above the upper")
                    return "upper"
                else:
                    log.write("Filtering FASTA's outside the upper and lower percentiles")
                    return False
            elif homogeneity_rating == 3:
                log.write("Dataset is homogenous")
                return True

        log.write("################################################################################")

    failed_N_count_df = stats[stats["N_Count"] >= max_n_count ] # Initialize the failed DataFrame

    # Calculate percentiles
    quantiles_stats = passed_df.quantile([lower_percentile, upper_percentile])
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]

    rating = calculate_homogeneity(axis, passed_df)
    if not rating:
        passed_df = passed_df[(passed_df[axis] > lower_percentiles[axis]) & (passed_df[axis] < upper_percentiles[axis])]
        failed_lower_df = passed_df[(passed_df[axis] <= lower_percentiles[axis])]
        failed_upper_df = passed_df[(passed_df[axis] >= upper_percentiles[axis])]
    elif rating == "lower":
        # the data are homogenous at the lower end but not the upper end
        # filter out the upper end
        passed_df = passed_df[(passed_df[axis] < upper_percentiles[axis])]
        failed_upper_df = passed_df[(passed_df[axis] >= upper_percentiles[axis])]
    elif rating == "upper":
        # the data are homogenous at the upper end but not the lower end
        # filter out the lower end
        passed_df = passed_df[(passed_df[axis] > lower_percentiles[axis])]
        failed_lower_df = passed_df[(passed_df[axis] <= lower_percentiles[axis])]
    elif rating:
        """
        if rating is True, filtering will only occur based on N Count
        """
        pass

    frames = [failed_N_count_df, failed_lower_df, failed_upper_df]
    failed_df = pd.concat(frames)
    failed_df.to_csv(os.path.join(species_dir, "failed.csv"))

    return passed_df

def generate_links_to_passed(passed_df, passed_dir, species_dir):

    filter_log = os.path.join(species_dir, "filter_log.txt")
    with open(filter_log, "a") as log:
        log.write("{} FASTA's passed".format(len(passed_df)))

    for source in passed_df.index:
        print("{} passed.".format(source.split("/")[-1]))
        dst = os.path.join(passed_dir, source.split("/")[-1])
        print("Linking to {}".format(dst))
        os.link(source, dst)

def bool_df_for_failed(species_dir, stats, max_n_count, max_contigs, lower_percentile, upper_percentile):

    quantiles_stats = stats.quantile([lower_percentile, upper_percentile])
    lower_percentiles = quantiles_stats.loc[lower_percentile]
    upper_percentiles = quantiles_stats.loc[upper_percentile]

    n_count_bool = (stats.iloc[:]["N_Count"]) <= max_n_count
    contigs_bool = (stats.iloc[:]["Contigs"]) <= max_contigs

    distances_bool = (stats.iloc[:]["Avg_Distances"] >= lower_percentiles["Avg_Distances"]) \
                   & (stats.iloc[:]["Avg_Distances"] <= upper_percentiles["Avg_Distances"])

    lengths_bool = (stats.iloc[:]["Assembly_Size"] >= lower_percentiles["Total_Length"]) \
                   & (stats.iloc[:]["Assembly_Size"] <= upper_percentiles["Total_Length"])

    passed = os.listdir(os.path.join(species_dir, passed))
    total = os.listdir(species_dir)
    filter_log = os.path.join(species_dir, "filter_log.txt")
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
    bool_df.to_csv(os.path.join(species_dir, "failed_tf.csv"))
    stats.to_csv(os.path.join(species_dir, "failed.csv"))

def make_passed_dir(species_dir):
    
    passed_dir = os.path.join(species_dir, "passed")
    if os.path.isdir(passed_dir):
        rmtree(passed_dir)
    if not os.path.isdir(passed_dir):
        os.mkdir(passed_dir)

    return passed_dir
    
def Main():
    parser = argparse.ArgumentParser(description = "Assess the integrity of your FASTA collection")
    parser.add_argument("species_dir", help = "directory containing your FASTA files")
    parser.add_argument("-d", "--directory", help = "Specifies that `species_dir` is a parent directory with \
            subdirectories for each species", action="store_true")
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

    species_dir = args.species_dir

    max_ns = args.max_n_count

    if args.from_list:
        for name in args.from_list:
            species_dir = os.path.join(species_dir, name)
            distance_matrix = os.path.join(species_dir, "distance_matrix.csv")
            mash(species_dir)
            stats = generate_fasta_stats(species_dir, distance_matrix)
            passed_df = assess_stats(species_dir, stats, args.max_n_count, args.max_contigs, args.lower, args.upper, args.multiplier, args.deviation_type)
            passed_dir = make_passed_dir(species_dir)
            generate_links_to_passed(passed_df, passed_dir, species_dir)
    else:
        for name in os.listdir(species_dir):
            species_dir = os.path.join(species_dir, name)
            if os.path.isdir(species_dir):
                contents = os.listdir(species_dir) # pass if there are <= 5 FASTA's
                if len(contents) <= 5:
                    continue
                else:
                    clean_up(species_dir)
                    distance_matrix = mash(species_dir)
                    stats = generate_fasta_stats(species_dir, distance_matrix)
                    for num in [2.5]:
                        filter_med_ad(species_dir, stats, multiplier=num, max_n_count=100)
                  # passed_dir = make_passed_dir(species_dir)
                  # passed_df = assess_stats(species_dir, stats, args.max_n_count, args.max_contigs, args.lower, args.upper, args.multiplier, args.deviation_type)
                  # generate_links_to_passed(passed_df, passed_dir, species_dir)

Main()
