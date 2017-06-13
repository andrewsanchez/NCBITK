#!/usr/bin/env python

import os
import re
import argparse
import subprocess

from urllib.request import urlretrieve
from urllib.error import URLError
from ftplib import error_temp
from time import strftime, sleep


def grab_zipped_genome(genbank_mirror,
                       species,
                       genome_id,
                       genome_url,
                       ext=".fna.gz"):
    """
    Download compressed genome from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/
    """

    zipped_path = "{}_genomic{}".format(genome_id, ext)
    zipped_url = "{}/{}".format(genome_url, zipped_path)
    zipped_dst = os.path.join(genbank_mirror, species, zipped_path)
    urlretrieve(zipped_url, zipped_dst)


def get_genome_id_and_url(assembly_summary, accession):

    genome_id = assembly_summary.ftp_path[accession].split('/')[-1]
    genome_url = assembly_summary.ftp_path[accession]

    return genome_id, genome_url


def sync_latest_genomes(genbank_mirror, assembly_summary, new_genomes, logger):

    for accession in new_genomes:
        genome_id, genome_url = get_genome_id_and_url(assembly_summary,
                                                      accession)
        species = assembly_summary.scientific_name.loc[accession]
        try:
            grab_zipped_genome(genbank_mirror, species, genome_id, genome_url)
            logger.info("Downloaded {}".format(genome_id))
        except error_temp as e:
            logger.info('error_temp for {}\n{}'.format(genome_id, e))
            sleep(2)
            grab_zipped_genome(genbank_mirror, species, genome_id, genome_url)
            logger.info("Downloaded {}".format(genome_id))
        except URLError as e:
            logger.info('URLError[1] for {}\n{}'.format(genome_id, e))
            grab_zipped_genome(
                genbank_mirror,
                species,
                genome_id,
                genome_url,
                ext=".fasta.gz")
            logger.info("Downloaded {}".format(genome_id))
        except URLError as e:
            logger.info('URLError[2] for {}\n{}'.format(genome_id, e))
            continue


def write_ftp_paths(genbank_mirror, assembly_summary, new_genomes):

    ftp_paths_file = os.path.join(genbank_mirror, '.info', 'ftp_paths.txt')

    if os.path.isfile(ftp_paths_file):
        os.remove(ftp_paths_file)

    with open(ftp_paths_file, 'a') as f:
        for accession in new_genomes:
            genome_url = assembly_summary.ftp_path[accession]
            genome_url = re.sub(r'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/', '',
                                genome_url)
            genome_parent_dir = genome_url.split('/')[-1]
            genome_url = '{}/{}_genomic.fna.gz'.format(genome_url,
                                                       genome_parent_dir)
            f.write(genome_url)
            f.write('\n')

    return ftp_paths_file


def rsync_latest_genomes(genbank_mirror, assembly_summary, new_genomes):

    ftp_paths_file = write_ftp_paths(genbank_mirror, assembly_summary,
                                     new_genomes)
    rsync_log = os.path.join(genbank_mirror, '.info',
                             'rsync_{}.out'.format(strftime('%Y.%m.%d.%H:%M')))
    incoming = os.path.join(genbank_mirror, 'incoming')
    if not os.path.isdir(incoming):
        os.mkdir(incoming)

    cmd = 'rsync --chmod=ugo=rwX --times --progress --itemize-changes --stats --files-from={}\
    --log-file={} --prune-empty-dirs ftp.ncbi.nlm.nih.gov::genomes/all/ {}'.format(
        ftp_paths_file, rsync_log, incoming)

    subprocess.Popen(cmd, shell=True).wait()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "genbank_mirror", help="Directory to save fastas", type=str)
    parser.add_argument("-u", "--unzip", action="store_true")
    args = parser.parse_args()
    genbank_mirror = args.genbank_mirror

    assembly_summary = get_assembly_summary(
        genbank_mirror,
        assembly_summary_url=
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
    )
    unzip_genbank_mirror(genbank_mirror)
    rename(genbank_mirror, assembly_summary)


if __name__ == "__main__":
    main()
