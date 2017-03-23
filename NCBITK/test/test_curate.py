from NCBITK import config
from NCBITK import curate
from NCBITK import get_resources
import unittest
import os
import shutil
import pandas as pd


class TestCurate(unittest.TestCase):

    def setUp(self):

        self.genbank_mirror = "NCBITK/test/resources/genbank_mirror"

        if os.path.isdir(self.genbank_mirror):
            shutil.rmtree(self.genbank_mirror)

        self.path_vars = config.instantiate_path_vars(self.genbank_mirror)
        info_dir, slurm, out, logger = self.path_vars
        shutil.copy('NCBITK/test/resources/assembly_summary.txt', info_dir)
        shutil.copy('NCBITK/test/resources/names.dmp', info_dir)

        self.assembly_summary = get_resources.get_assembly_summary(self.genbank_mirror, False)
        self.names = get_resources.get_scientific_names(self.genbank_mirror, self.assembly_summary, fetch_new=False)
        self.assembly_summary = get_resources.update_assembly_summary(self.genbank_mirror, self.assembly_summary, self.names)
        # assembly_summary = get_resources.get_resources(genbank_mirror, logger, False)

    def test_create_species_dirs(self):

        info_dir, slurm, out, logger = self.path_vars
        curate.create_species_dirs(self.genbank_mirror, self.assembly_summary, logger)
        self.assertEqual(len(os.listdir(self.genbank_mirror)), len(set(self.assembly_summary.scientific_name.tolist())))

    def tearDown(self):
        shutil.rmtree(self.genbank_mirror)


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
