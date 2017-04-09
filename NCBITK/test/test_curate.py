from NCBITK import config
from NCBITK import curate
from NCBITK import sync
from NCBITK import get_resources

import unittest
import os
import tempfile
import shutil
import pandas as pd


class TestCurate(unittest.TestCase):

    def setUp(self):

        # self.genbank_mirror = tempfile.mkdtemp()
        self.genbank_mirror = '/Users/andrew/scratch/genbank'
        self.assembly_summary = pd.read_csv('NCBITK/test/resources/assembly_summary.txt', sep="\t", index_col=0)
        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        self.info_dir, self.slurm, self.out, self.logger = self.path_vars

        self.test_species = 'Acinetobacter_nosocomialis'
        self.test_genomes = self.assembly_summary.index[self.assembly_summary.scientific_name == self.test_species]
        self.species_list = curate.get_species_list(self.assembly_summary, [self.test_species])

        self.genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        self.local_genomes, self.new_genomes,\
        self.old_genomes, self.sketch_files,\
        self.missing_sketch_files = self.genbank_assessment

    def test_create_species_dirs_all(self):

        species_list = curate.get_species_list(self.assembly_summary, 'all')
        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, species_list)
        local_species = os.listdir(self.genbank_mirror)
        local_species.remove('.info')

        self.assertEqual(len(local_species), len(species_list))

    def test_create_species_dirs_list(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)
        local_species = os.listdir(self.genbank_mirror)
        local_species.remove('.info')

        self.assertEqual(len(local_species), len(self.species_list))

    def test_assess_fresh(self):

        self.assertTrue(len(self.new_genomes) == len(self.test_genomes))
        self.assertTrue(len(self.missing_sketch_files) == len(self.test_genomes))
        self.assertTrue(len(self.sketch_files) == 0)
        self.assertTrue(len(self.local_genomes) == 0)
        self.assertTrue(len(self.old_genomes) == 0)

    def test_sync_latest_genomes(self):

        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, self.species_list)

        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

        print(self.genbank_mirror)
        sync.sync_latest_genomes(self.genbank_mirror, self.assembly_summary, new_genomes, self.logger)
        print(self.genbank_mirror)
        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, self.species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment

        self.assertTrue(len(local_genomes) == len(self.test_genomes))
        self.assertTrue(len(new_genomes) == 0)
        self.assertTrue(len(old_genomes) == 0)
        self.assertTrue(len(missing_sketch_files) == len(self.test_genomes))
        self.assertTrue(len(sketch_files) == 0)

    # def test_assess_after(self):

    #     self.assertTrue(len(self.new_genomes) == len(self.test_genomes))
    #     self.assertTrue(len(self.missing_sketch_files) == len(self.test_genomes))
    #     self.assertTrue(len(self.sketch_files) == 0)
    #     self.assertTrue(len(self.local_genomes) == 0)

    # def test_update_genbank(self):
    #     None

    # def test_remove_old_genomes(self):
    #     None

    # def test_get_resources(self):
    #     None

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)

class TestArrays(unittest.TestCase):
    None

if __name__ == '__main__':
    unittest.main()
