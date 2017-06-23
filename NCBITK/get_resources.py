#!/usr/bin/env python

import os
import logging
import subprocess
import pandas as pd
import tarfile
from urllib.request import urlretrieve

bacteria_assembly_summary = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
taxdump_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


def get_assembly_summary(update,
                         assembly_summary_url=bacteria_assembly_summary):
    """Get current version of assembly_summary.txt and load into DataFrame"""

    if update:
        assembly_summary = pd.read_csv(
            bacteria_assembly_summary, sep="\t", index_col=0, skiprows=1)
    else:
        assembly_summary = pd.read_csv(
            path_assembly_summary, sep="\t", index_col=0)

    return assembly_summary


def get_scientific_names(genbank_mirror, assembly_summary, update=True):
    """
    Get names.dmp from the taxonomy dump
    """

    info_dir = os.path.join(genbank_mirror, ".info")
    names_dmp = os.path.join(genbank_mirror, ".info", 'names.dmp')

    # TODO: Create else statement to read local names.dmp
    if update:
        taxdump = urlretrieve(taxdump_url)
        taxdump_tar = tarfile.open(taxdump[0])
        taxdump_tar.extract('names.dmp', info_dir)

        # TODO: Use fileinput instead of sed
        sed_cmd = "sed -i '/scientific name/!d' {}".format(
            names_dmp)  # we only want rows with the scientific name
        subprocess.Popen(sed_cmd, shell='True').wait()
        names = pd.read_csv(
            names_dmp, sep='\t', index_col=0, header=None, usecols=[0, 2])
        names = names.loc[set(assembly_summary.species_taxid.tolist())]
        names.index.name = 'species_taxid'
        names.columns = ['scientific_name']
        names.scientific_name.replace({' ': '_'}, regex=True, inplace=True)
        names.scientific_name.replace({'/': '_'}, regex=True, inplace=True)
        names.to_csv(names_dmp)
    else:
        names = pd.read_csv(names_dmp, index_col=0)

    return names


def update_assembly_summary(assembly_summary, names):

    # TODO: Maybe keep separate from assembly_summary?
    # TODO: Very slow
    for taxid in names.index:
        scientific_name = names.scientific_name.loc[taxid]
        # get the list of indices that share the same species_taxid
        ixs = assembly_summary.index[assembly_summary.species_taxid ==
                                     taxid].tolist()
        assembly_summary.loc[ixs, 'scientific_name'] = scientific_name

    return assembly_summary


def clean_up_assembly_summary(genbank_mirror, assembly_summary):

    non_word_chars = [
        '?', '=', '"', '!', '*', '>', '+', ':', '\x93', '´', '@', '“', ']',
        '-', '’', ',', ')', '”', '&', '\x94', '\u200d', '/', '.', ';', '#',
        '(', '%', '[', "'"
    ]

    # If both infraspecific_name and isolate are empty, fill infraspecific_name with "NA"
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].isnull())].fillna('NA'))

    # If infraspecific_name column is empty and isolate column is not empty, fill infraspecific_name with the value of isolate.
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].notnull())].fillna(assembly_summary['isolate']))

    # TODO: Rewrite using fileinput
    assembly_summary.assembly_level.replace(
        {
            ' ': '_'
        }, regex=True, inplace=True)
    assembly_summary.organism_name.replace(
        {
            ' ': '_'
        }, regex=True, inplace=True)
    assembly_summary.organism_name.replace(
        {
            '[\W]': '_'
        }, regex=True, inplace=True)
    assembly_summary.infraspecific_name.replace(
        {
            '[\W]': '_'
        }, regex=True, inplace=True)

    # Extract organism_name column

    updated_assembly_summary = os.path.join(genbank_mirror, '.info',
                                            'assembly_summary.txt')
    assembly_summary.to_csv(updated_assembly_summary, sep='\t')
    return assembly_summary


def get_resources(genbank_mirror, update):
    """
    Get assembly summary and taxonomy dump file for bacteria.
    Parse and load into Pandas DataFrames.
    """

    if update:
        assembly_summary = get_assembly_summary(genbank_mirror, update)
        names = get_scientific_names(genbank_mirror, assembly_summary)
        assembly_summary = update_assembly_summary(genbank_mirror,
                                                   assembly_summary, names)
        assembly_summary = clean_up_assembly_summary(genbank_mirror,
                                                     assembly_summary)
    else:
        assembly_summary = get_assembly_summary(genbank_mirror, update)

    return assembly_summary
