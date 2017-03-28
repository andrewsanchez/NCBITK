from NCBITK import config
from NCBITK import curate
from NCBITK import get_resources

import unittest
import os
import tempfile
import shutil
import pandas as pd


class TestCurate(unittest.TestCase):

    def setUp(self):

        self.genbank_mirror = tempfile.mkdtemp()
        self.assembly_summary = pd.read_csv('NCBITK/test/resources/assembly_summary.txt', sep="\t", index_col=0)
        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        self.info_dir, self.slurm, self.out, self.logger = self.path_vars

        self.species_list_slice = ['Acinetobacter_nosocomialis',
                        'Escherichia_coli',
                        'Bacillus_anthracis']

    def test_create_species_dirs_all(self):

        species_list = curate.get_species_list(self.assembly_summary, 'all')
        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, species_list)
        total_species_assembly_summary = len(set(self.assembly_summary.scientific_name.tolist()))
        total_species_local = len(os.listdir(self.genbank_mirror))
        self.assertEqual(total_species_local, total_species_assembly_summary)

    def test_create_species_dirs_list(self):

        species_list = self.species_list_slice
        species_list = curate.get_species_list(self.assembly_summary, species_list)
        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, self.logger, species_list)
        total_species_local = len([i for i in os.listdir(self.genbank_mirror) if not i.startswith('.')])
        self.assertEqual(total_species_local, len(species_list))

    def test_assess_genbank_mirror(self):

        species_list = self.species_list_slice
        # local_genomes = curate.get_local_genomes(self.genbank_mirror)
        # new_genomes = curate.get_new_genome_list(self.genbank_mirror, self.assembly_summary, local_genomes, species_list)
        # sketch_files = curate.get_sketch_files(self.genbank_mirror)
        # missing_sketch_files = curate.get_missing_sketch_files(local_genomes, new_genomes, sketch_files)
        # old_genomes = curate.get_old_genomes(self.genbank_mirror, self.assembly_summary, local_genomes)
        genbank_assessment = curate.assess_genbank_mirror(self.genbank_mirror, self.assembly_summary, species_list)
        local_genomes, new_genomes, old_genomes, sketch_files, missing_sketch_files = genbank_assessment
        self.assertTrue(len(new_genomes) > len(local_genomes))

        self.assertTrue(len(local_genomes) == len(sketch_files))
        self.assertTrue(len(missing_sketch_files) > len(sketch_files))

    def test_sync_latest_genomes(self):
        None

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)

class TestArrays(unittest.TestCase):
    None

#     def test_assess_genbank_mirror(self):
#         None

#     def test_remove_old_genomes(self):
#         None

#     def test_sync_latest_genomes(self):
#         None

# class TestGetResources(unittest.TestCase):

#     def test_get_resources(self):
#         None

if __name__ == '__main__':
    unittest.main()
