import glob
import os
import shutil
import tempfile
import unittest

import pandas as pd

from NCBITK import config, curate, get_resources, sync


class TestCurate(unittest.TestCase):
    def setUp(self):

        self.genbank_mirror = tempfile.mkdtemp(prefix='Genbank_')
        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        self.incoming = os.path.join(self.genbank_mirror, 'incoming')
        self.info_dir, self.slurm, self.out, self.logger = self.path_vars
        self.original_assembly_summary = pd.read_csv(
            'NCBITK/test/resources/original_assembly_summary.txt',
            sep="\t",
            index_col=0)
        self.updated_assembly_summary = pd.read_csv(
            'NCBITK/test/resources/updated_assembly_summary.txt',
            sep="\t",
            index_col=0)
        get_resources.clean_up_assembly_summary(
            self.updated_assembly_summary)
        self.updated_assembly_summary_len = len(
            self.updated_assembly_summary.index)
        self.test_species = 'Acinetobacter_nosocomialis'
        self.test_genomes = self.updated_assembly_summary.index[
            self.updated_assembly_summary.scientific_name == self.test_species]
        self.species_list = curate.get_species_list(
            self.updated_assembly_summary, self.test_species)
        self.species_dir = os.path.join(self.genbank_mirror, self.test_species)
        self.all_species_from_assembly_summary = self.updated_assembly_summary.scientific_name[
            self.updated_assembly_summary.scientific_name.notnull()]
        self.all_species_from_assembly_summary = set(
            self.all_species_from_assembly_summary.tolist())
        # excludes genomes with NaN value for scientific name
        self.all_genomes_from_assembly_summary = self.updated_assembly_summary.index[
            self.updated_assembly_summary.scientific_name.notnull()]

        os.mkdir(self.incoming)

    def test_get_species_list(self):

        species_list = ['Bacillus_anthracis', 'Escherichia_coli']
        species_from_string = curate.get_species_list(
            self.updated_assembly_summary, self.test_species)
        species_from_list = curate.get_species_list(
            self.updated_assembly_summary, species_list)
        all_species = curate.get_species_list(self.updated_assembly_summary,
                                              'all')

        self.assertTrue(len(species_from_string), 1)
        self.assertTrue(len(species_from_list), len(species_list))
        self.assertTrue(
            len(all_species), len(self.all_species_from_assembly_summary))

    def test_create_species_dirs_all(self):

        species_list = curate.get_species_list(self.updated_assembly_summary,
                                               'all')
        curate.create_species_dirs(self.genbank_mirror, self.logger,
                                   species_list)
        local_species = os.listdir(self.genbank_mirror)
        local_species.remove('.info')
        local_species.remove('incoming')

        self.assertEqual(len(local_species), len(species_list))

    def test_assess_fresh(self):

        genbank_assessment = curate.assess_genbank_mirror(
            self.genbank_mirror, self.updated_assembly_summary,
            self.all_species_from_assembly_summary, self.logger)

        local_genomes, new_genomes, old_genomes = genbank_assessment

        self.assertTrue(
            len(new_genomes) == len(self.all_genomes_from_assembly_summary))
        self.assertTrue(len(local_genomes) == 0)
        self.assertTrue(len(old_genomes) == 0)

    def test_assess_changes(self):

        curate.create_species_dirs(self.genbank_mirror, self.logger,
                                   self.species_list)

        genbank_assessment = curate.assess_genbank_mirror(
            self.genbank_mirror, self.updated_assembly_summary,
            self.species_list, self.logger)

        local_genomes, new_genomes, old_genomes = genbank_assessment
        before_sync_new_genomes = new_genomes[:10]
        after_sync_new_genomes = new_genomes[10:]

        for genome in before_sync_new_genomes:
            dst = os.path.join(self.species_dir, genome)
            tempfile.mkstemp(
                prefix='{}.fasta'.format(genome), dir=self.species_dir)

        genbank_assessment = curate.assess_genbank_mirror(
            self.genbank_mirror,
            self.updated_assembly_summary.drop(before_sync_new_genomes[0]),
            self.species_list, self.logger)

        local_genomes, new_genomes, old_genomes = genbank_assessment

        self.assertTrue(len(local_genomes) == len(before_sync_new_genomes))
        self.assertFalse(len(new_genomes) == 0)
        self.assertTrue(len(new_genomes) == len(after_sync_new_genomes))
        self.assertTrue(len(old_genomes) == 1)

    # def test_get_old_genomes(self):

    #     local_genomes = self.test_genomes
    #     not_in_assembly_summary = self.test_genomes[:5].tolist()
    #     self.updated_assembly_summary.drop(not_in_assembly_summary, inplace=True)

    #     old_genomes = curate.get_old_genomes(self.genbank_mirror, local_genomes)

    #     self.assertTrue(sorted(not_in_assembly_summary) == sorted(old_genomes))

    # def test_get_local_genomes(self):

    #     curate.create_species_dirs(self.genbank_mirror,
    #                                self.logger,
    #                                self.species_list)

    #     for genome in self.test_genomes:
    #         dst = os.path.join(self.species_dir, genome)
    #         tempfile.mkstemp(prefix=genome, dir=self.species_dir)

    #     local_genomes = curate.get_local_genomes(self.genbank_mirror)

    #     self.assertTrue(len(local_genomes) == len(self.test_genomes))

    # def test_post_rsync_cleanup(self):

    #     for genome in self.test_genomes:
    #         dst = os.path.join(self.species_dir, genome)
    #         tempfile.mkstemp(prefix='{}_'.format(genome), dir=self.incoming)

    #     curate.post_rsync_cleanup(self.genbank_mirror, self.updated_assembly_summary, self.logger)
    #     self.assertTrue(len(self.test_genomes) == len(os.listdir(self.species_dir)))
    #     self.assertFalse(os.listdir(self.incoming))

    # def test_sync_latest_genomes(self):

    #     curate.create_species_dirs(self.genbank_mirror, self.updated_assembly_summary, self.logger, self.species_list)

    #     genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.updated_assembly_summary, self.species_list)
    #     local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment
    #     before_sync_new_genomes = new_genomes

    #     sync.sync_latest_genomes(self.genbank_mirror, self.updated_assembly_summary, new_genomes, self.logger)

    #     species_dir = os.path.join(self.genbank_mirror, self.test_species)

    #     for f in glob.glob('{}/GCA*'.format(species_dir)):
    #         self.assertTrue(os.path.isfile(f))

    #     genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.updated_assembly_summary, self.species_list)
    #     local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

    #     self.assertTrue(len(before_sync_new_genomes) == len(local_genomes))
    #     self.assertTrue(len(local_genomes) == len(self.test_genomes))
    #     self.assertTrue(len(missing_sketch_files) == len(self.test_genomes))
    #     self.assertTrue(len(new_genomes) == 0)
    #     self.assertTrue(len(old_genomes) == 0)
    #     self.assertTrue(len(sketch_files) == 0)

    def test_rename(self):
        species_list = curate.get_species_list(self.assembly_summary, 'all')
        correct_name = 'GCA_000007365.1_Buchnera_aphidicola_Sg_Schizaphis_graminum'
        for genome in self.all_genomes_from_assembly_summary:
            tempfile.mkstemp(
                prefix=genome, suffix='.fasta', dir=self.genbank_mirror)
        curate.rename(self.genbank_mirror, self.assembly_summary)
        renamed_genomes = len([
            x for x in os.listdir(self.genbank_mirror) if x.startswith('GCA')
        ])
        self.assertEqual(self.assembly_summary_len, renamed_genomes)

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)


if __name__ == '__main__':
    unittest.main()
