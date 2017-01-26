#!/usr/bin/env python

import os, re, argparse
import pandas as pd
from urllib.request import urlretrieve

def get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"):

    """Get current version of assembly_summary.txt and load into DataFrame"""

    assembly_summary_dst = os.path.join(genbank_mirror, ".info", "assembly_summary.txt")
    urlretrieve(assembly_summary_url, assembly_summary_dst)
    assembly_summary = pd.read_csv(assembly_summary_dst, sep="\t", index_col=0, skiprows=1)

    return assembly_summary

def unzip_genome(root, f, genome_id):

    """
    Decompress genome and remove the compressed genome.
    """

    zipped_src = os.path.join(root, f)
    zipped = gzip.open(zipped_src)
    decoded = zipped.read()
    unzipped = "{}.fasta".format(genome_id)
    unzipped = os.path.join(root, unzipped)
    unzipped = open(unzipped, "wb")
    unzipped.write(decoded)
    zipped.close()
    unzipped.close()
    os.remove(zipped_src)

def unzip_genbank_mirror(genbank_mirror):

    for root, files, dirs, in os.walk(genbank_mirror):
        for f in files:
            if f.endswith("gz"):
                genome_id = "_".join(f.split("_")[:2])
                unzip_genome(root, f, genome_id)

def rm_duplicates(seq):

    """
    remove duplicate strings during renaming
    """

    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def rename(target_dir, assembly_summary):

    """
    Clean up assembly_summary.txt and renamed FASTA's.
    """

    # If infraspecific_name and isolate columns are empty, fill infraspecific_name with "NA"
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].isnull())].fillna('NA'))

    # If infraspecific_name column is empty and isolate column is not empty, fill infraspecific_name with the value of isolate.
    assembly_summary.update(assembly_summary['infraspecific_name'][(assembly_summary['infraspecific_name'].isnull()) &\
            (assembly_summary['isolate'].notnull())].fillna(assembly_summary['isolate']))

    assembly_summary.assembly_level.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({' ': '_'}, regex=True, inplace=True)
    assembly_summary.organism_name.replace({'[\W]': '_'}, regex=True, inplace=True)
    assembly_summary.infraspecific_name.replace({'[\W]': '_'}, regex=True, inplace=True)

    for root, dirs, files in os.walk(target_dir):
        for f in files:
            if f.endswith("fasta"):
                assembly_summary_id = "_".join(f.split('_')[0:2])
                if assembly_summary_id in assembly_summary.index:
                    org_name = assembly_summary.get_value(assembly_summary_id, 'organism_name')
                    strain = assembly_summary.get_value(assembly_summary_id, 'infraspecific_name')
                    assembly_level  = assembly_summary.get_value(assembly_summary_id, 'assembly_level')
                    new_name = '{}_{}_{}_{}.fasta'.format(assembly_summary_id, org_name, strain, assembly_level)
                    rm_words = re.compile( r'((?<=_)(sp|sub|substr|subsp|str|strain)(?=_))' )
                    new_name = rm_words.sub('_', new_name)
                    new_name = re.sub(r'_+', '_', new_name )
                    new_name = new_name.split('_')
                    new_name = rm_duplicates(new_name)
                    new_name = '_'.join(new_name)
                    old = os.path.join(root, f)
                    new = os.path.join(root, new_name)
                    os.rename(old, new)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dir', help = 'The folder whose contents will be renamed', type=str)
    parser.add_argument('-s', '--source', help = 'Specify a directory to rename.', action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.target_dir

    assembly_summary = get_assembly_summary(genbank_mirror, assembly_summary_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt")
    unzip_genbank_mirror(genbank_mirror)
   #rename(genbank_mirror, assembly_summary)

if __name__ == '__main__':
    main()
